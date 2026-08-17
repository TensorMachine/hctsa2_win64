function m = BF_ar(y,n,approach)
% BF_ar   Estimate a scalar AR model (toolbox-free System Identification ar()).
%
% Fits  A(q) y(t) = e(t),  A(q) = 1 + a_1 q^-1 + ... + a_n q^-n.
%
% approach: 'fb' forward-backward least squares (MATLAB's default), 'ls'
% forward-only least squares, 'yw' Yule-Walker. Only 'fb' is exercised by the
% hctsa library; the others are provided because they are cheap and verified.
% Burg is deliberately NOT provided rather than shipped unverified.
%
% Returns a struct mirroring the parts of an idpoly object that hctsa reads:
% A (and its lowercase alias a), NoiseVariance, np, N, Ts and EstimationInfo
% with LossFcn and FPE.

if nargin < 3 || isempty(approach), approach = 'fb'; end
y = y(:);
y = y(isfinite(y));
N = numel(y);
if n >= N
    error('BF_ar:orderTooLarge','AR order %u too large for %u samples.',n,N);
end

switch lower(approach)
case {'fb','ls'}
    % Forward equations: y(t) = -[y(t-1)...y(t-n)] a
    rows = N - n;
    Xf = zeros(rows,n);
    for k = 1:n
        Xf(:,k) = y(n+1-k : N-k);
    end
    yf = y(n+1:N);
    if strcmpi(approach,'fb')
        % Backward equations: y(t) = -[y(t+1)...y(t+n)] a
        Xb = zeros(rows,n);
        for k = 1:n
            Xb(:,k) = y(1+k : N-n+k);
        end
        yb = y(1:N-n);
        X = [Xf; Xb];
        b = [yf; yb];
    else
        X = Xf; b = yf;
    end
    a = -(X\b);

case 'yw'
    r = zeros(n+1,1);
    yc = y - mean(y);
    for k = 0:n
        r(k+1) = sum(yc(1+k:end).*yc(1:end-k))/N;
    end
    R = toeplitz(r(1:n));
    a = -(R\r(2:n+1));

otherwise
    error('BF_ar:unknownApproach','Unknown AR approach ''%s''.',approach);
end

m.A = [1, a(:)'];
% MATLAB's idpoly exposes the polynomials in both cases; hctsa reads `m.a` in
% MF_FitSubsegments and `m.A` elsewhere, so provide both.
m.a = m.A;
e = filter(m.A,1,y);

% Loss function and Akaike's Final Prediction Error.
%
% The loss is the mean squared prediction error averaged over ALL N samples with
% the startup transient zeroed -- i.e. sum(e(n+1:end).^2)/N, not /(N-n). This is
% the same convention as pe()'s estimated initial conditions; dividing by (N-n)
% instead is wrong by a factor (N-n)/N, which showed up as a 4.9e-3 error in
% MF_FitSubsegments' fpe fields on ~410-sample segments.
e(1:min(n,numel(e))) = 0;
m.NoiseVariance = sum(e.^2)/N;
m.EstimationInfo.LossFcn = m.NoiseVariance;
m.EstimationInfo.FPE = m.NoiseVariance*(1 + n/N)/(1 - n/N);
m.order = n;
m.np = n;      % number of estimated parameters, for BF_aic
m.N = N;       % number of data samples, for BF_aic
m.Ts = 1;
end
