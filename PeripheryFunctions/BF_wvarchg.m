function [pts_Opt,kopt,t_est] = BF_wvarchg(y,K,d)
% BF_wvarchg   Find variance change points (wvarchg()).
%
% Toolbox-free replacement, written against MathWorks' wvarchg.m source.
%
% The change points come from dynamic-programming minimisation of the
% Gaussian log-likelihood contrast, which is standard. The part that is NOT
% standard -- and which no threshold rule reproduces -- is how the number of
% change points is chosen. MathWorks does not threshold the contrast curve at
% all. Instead it builds the intervals of penalty values over which each
% candidate k would be selected, and takes the k owning the **widest** such
% interval, subject to that width exceeding the last upper bound.

if nargin < 2 || isempty(K), K = 6;  end
if nargin < 3 || isempty(d), d = 10; end

y = y(:) - mean(y);
K = K + 1;

%-------------------------------------------------------------------------------
%% Contrast matrix and dynamic programming
%-------------------------------------------------------------------------------
N  = length(y);
y2 = y.^2;
matD = NaN(N,N);
for i = 1:N-d
    vi = 1:N-i+1;
    dummy = vi.*log(cumsum(y2(i:N))'./vi);
    matD(i,i+d-1:N) = dummy(d:end);
end
% isinf(-M) is identical to isinf(M) for every input including NaN, so the
% negation only builds a full N-by-N temporary for nothing. The sign flip on the
% selected entries is kept exactly as MathWorks has it: it maps +Inf to -Inf and
% -Inf to +Inf, which is not the same as forcing everything non-finite to +Inf.
ind = isinf(matD);   matD(ind) = -matD(ind);
ind = isnan(matD);   matD(ind) = Inf;

I = zeros(K,N);
I(1,:) = matD(1,:);
t = zeros(K,N);
if K > 2
    for k = 2:K-1
        for L = k:N
            [I(k,L),t(k-1,L)] = min(I(k-1,1:L-1) + matD(2:L,L)');
        end
    end
end
t_est = diag(ones(1,K)*N);
[I(K,N),t(K-1,N)] = min(I(K-1,1:N-1) + matD(2:N,N)');
for j = 2:K
    for k = j-1:-1:1
        col = t_est(j,k+1);
        if col > 0
            t_est(j,k) = t(k,col);
        end
    end
end

%-------------------------------------------------------------------------------
%% Number of change points, by widest penalty interval
%-------------------------------------------------------------------------------
V  = I(:,N);
g2 = zeros(1,K);
for j = 2:K
    g2(j) = min((V(1:j-1)-V(j))./(j-1:-1:1)');
end
k = 0;
G2 = []; M = [];
for j = 2:K
    if g2(j) > max(g2(j+1:K))
        k = k+1; G2(k) = g2(j); M(k) = j;
    end
end
M(k+1)  = K;
G2(k+1) = g2(K);
M  = M(:);
G2 = G2(:);

G1 = [G2(1:k+1); 0];
G2 = [Inf; G2];
M  = [1; M] - 1;

% G1 and G2 hold the lower and upper bounds of the penalty intervals; M(i) is
% the number of change points selected for a penalty inside interval i.
if length(G1) == 2
    kopt = 0;
    pts_Opt = [];
else
    [lmax,indopt] = max(G2(2:end-1)-G1(2:end-1));
    kopt = M(indopt+1)*(lmax > G2(end));
    pts_Opt = t_est(kopt+1,1:kopt);
end
end
