function [h,pValue,stat,cValue] = BF_lbqtest(res,varargin)
% BF_lbqtest   Ljung-Box Q-test (toolbox-free lbqtest()).
%
%   Q = T(T+2) * sum_{k=1}^{L} r_k^2/(T-k),  Q ~ chi2(L - dof)
%
% Name/value pairs: 'lags' (default min(20,T-1)), 'alpha' (0.05), 'dof' (lags).

res = res(:);
res = res(isfinite(res));
T = numel(res);

lags = min(20,T-1);
alpha = 0.05;
dof = [];
for k = 1:2:numel(varargin)-1
    switch lower(varargin{k})
    case 'lags',  lags  = varargin{k+1};
    case 'alpha', alpha = varargin{k+1};
    case 'dof',   dof   = varargin{k+1};
    end
end
if isempty(dof), dof = lags; end

lags = lags(:)'; alpha = alpha(:)'; dof = dof(:)';
n = max([numel(lags) numel(alpha) numel(dof)]);
if numel(lags)==1,  lags  = repmat(lags,1,n);  end
if numel(alpha)==1, alpha = repmat(alpha,1,n); end
if numel(dof)==1,   dof   = repmat(dof,1,n);   end

% Sample autocorrelations of the (mean-removed) series
r = res - mean(res);
maxLag = max(lags);
c0 = sum(r.^2);
ac = zeros(1,maxLag);
for k = 1:maxLag
    ac(k) = sum(r(1+k:end).*r(1:end-k))/c0;
end

stat = zeros(1,n); pValue = zeros(1,n); cValue = zeros(1,n); h = false(1,n);
for i = 1:n
    L = lags(i);
    stat(i) = T*(T+2)*sum(ac(1:L).^2./(T-(1:L)));
    pValue(i) = gammainc(stat(i)/2,dof(i)/2,'upper');  % accurate chi2 upper tail
    cValue(i) = chi2inv(1-alpha(i),dof(i));
    h(i) = pValue(i) < alpha(i);
end
end
