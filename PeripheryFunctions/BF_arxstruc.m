function V = BF_arxstruc(ze,zv,NN)
% BF_arxstruc   Loss functions for a set of model structures (arxstruc()).
%
% For each candidate order, an AR model is estimated on ZE by least squares and
% its one-step prediction error is evaluated on ZV. hctsa only uses the
% output-only (AR) case, so NN is a column of na values.
%
%---OUTPUT: V, in MATLAB's layout
%   V(1,1:m)   loss function for each structure, on the validation data
%   V(1,m+1)   number of validation samples
%   V(2,1:m)   the orders themselves
%
% The loss is sum over a COMMON validation range, shared by every structure so
% the values are comparable across orders, divided by the full validation length
% Nv. The common range starts at max(NN)+2, i.e. the first max(NN)+1 prediction
% errors are discarded.
%
% Both details were determined against reference output from a licensed MATLAB:
% using a per-order warm-up instead leaves maxv ~1e-2 off, and using max(NN)
% rather than max(NN)+1 leaves it ~1e-3 off. With this convention the endpoint
% losses agree to ~1e-5.

ze = ze(:);
zv = zv(:);
NN = NN(:);
m = numel(NN);
Nv = numel(zv);
w = max(NN) + 1;

V = zeros(2,m+1);
for j = 1:m
    na = NN(j);

    % Least-squares AR fit on the estimation data (forward regression, which is
    % what arxstruc uses -- not the forward-backward default of ar()).
    % Every structure is fit over the SAME range so the candidates are
    % estimated on identical data, and that range starts at max(NN)+1 -- the
    % same offset used on the validation side. Determined against reference
    % output: a per-order range leaves the losses ~1.3e-5 off, max(NN) leaves
    % ~4e-6, and max(NN)+1 gives ~1e-14.
    Ne = numel(ze);
    s = w;   % same offset as the validation range: max(NN)+1
    rows = Ne - s;
    X = zeros(rows,na);
    for k = 1:na
        X(:,k) = ze(s+1-k : Ne-k);
    end
    a = -(X\ze(s+1:Ne));

    % One-step prediction error on the validation data, common range
    A = [1, a(:)'];
    e = filter(A,1,zv);
    V(1,j) = sum(e(w+1:end).^2)/Nv;
    V(2,j) = na;
end
V(1,m+1) = Nv;
V(2,m+1) = 0;
end
