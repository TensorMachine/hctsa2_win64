function [h,pValue,stat,cValue] = BF_archtest(res,varargin)
% BF_archtest   Engle's ARCH test for conditional heteroskedasticity.
%
% Regresses squared residuals on their own lags; under the null of no ARCH
% effects, T*R^2 ~ chi2(lags).
%
% Name/value pairs: 'lags' (default 1), 'alpha' (0.05).

res = res(:); res = res(isfinite(res));
lags = 1; alpha = 0.05;
for k = 1:2:numel(varargin)-1
    switch lower(varargin{k})
    case 'lags',  lags  = varargin{k+1};
    case 'alpha', alpha = varargin{k+1};
    end
end
lags = lags(:)'; alpha = alpha(:)';
n = max(numel(lags),numel(alpha));
if numel(lags)==1,  lags  = repmat(lags,1,n);  end
if numel(alpha)==1, alpha = repmat(alpha,1,n); end

% NOT demeaned: archtest assumes the input is already a residual series, and
% demeaning here is visibly wrong (30x larger error against the reference).
e2 = res.^2;
T = numel(e2);

stat = zeros(1,n); pValue = zeros(1,n); cValue = zeros(1,n);
for i = 1:n
    L = lags(i);
    yv = e2(L+1:end);
    % Preallocate rather than growing Xv one column per lag: same values, but
    % without reallocating and copying the whole matrix L times.
    Xv = ones(numel(yv),L+1);
    for j = 1:L
        Xv(:,j+1) = e2(L+1-j:end-j);
    end
    b = Xv\yv;
    r = yv - Xv*b;
    sst = sum((yv-mean(yv)).^2);
    R2 = 1 - sum(r.^2)/sst;
    stat(i) = numel(yv)*R2;
    pValue(i) = gammainc(stat(i)/2,L/2,'upper');  % accurate chi2 upper tail
    cValue(i) = chi2inv(1-alpha(i),L);
end
h = pValue < alpha;
end
