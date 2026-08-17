function [h,pValue,stat,cValue] = BF_kpsstest(y,varargin)
% BF_kpsstest   KPSS stationarity test (toolbox-free kpsstest()).
%
% Null hypothesis: the series is trend stationary.
%
%   e_t   residuals from regressing y on [1, t]   (or just 1 if trend is false)
%   S_t   cumsum(e)
%   s2    Newey-West long-run variance with Bartlett weights over `lags`
%   stat  sum(S_t^2) / (T^2 * s2)
%
% Name/value pairs: 'lags' (default 0), 'trend' (default true), 'alpha' (0.05).
%
% p-values are interpolated from the Kwiatkowski, Phillips, Schmidt & Shin (1992)
% Table 1 critical values and clipped to [0.01, 0.10], matching MATLAB (which
% warns econ:kpsstest:StatTooSmall / StatTooBig at the boundaries).

y = y(:);
y = y(isfinite(y));
T = numel(y);

lags = 0; trend = true; alpha = 0.05;
for k = 1:2:numel(varargin)-1
    switch lower(varargin{k})
    case 'lags',  lags  = varargin{k+1};
    case 'trend', trend = varargin{k+1};
    case 'alpha', alpha = varargin{k+1};
    end
end
lags = lags(:)'; alpha = alpha(:)'; trend = logical(trend(:))';
n = max([numel(lags) numel(alpha) numel(trend)]);
if numel(lags)==1,  lags  = repmat(lags,1,n);  end
if numel(alpha)==1, alpha = repmat(alpha,1,n); end
if numel(trend)==1, trend = repmat(trend,1,n); end

% KPSS (1992) Table 1, for alpha = [0.10 0.05 0.025 0.01]
alphaTable = [0.100 0.050 0.025 0.010];
cvLevel    = [0.347 0.463 0.574 0.739];   % trend = false
cvTrend    = [0.119 0.146 0.176 0.216];   % trend = true

stat = zeros(1,n); pValue = zeros(1,n); cValue = zeros(1,n);
for i = 1:n
    if trend(i)
        X = [ones(T,1), (1:T)'];
        cv = cvTrend;
    else
        X = ones(T,1);
        cv = cvLevel;
    end
    e = y - X*(X\y);
    S = cumsum(e);

    % Newey-West long-run variance, Bartlett kernel
    s2 = sum(e.^2)/T;
    for j = 1:lags(i)
        w = 1 - j/(lags(i)+1);
        s2 = s2 + 2*w*sum(e(j+1:end).*e(1:end-j))/T;
    end

    stat(i) = sum(S.^2)/(T^2*s2);

    % Interpolate the p-value, clipping at the table edges
    pValue(i) = interp1(cv,alphaTable,stat(i),'linear');
    if stat(i) <= cv(1),   pValue(i) = alphaTable(1); end
    if stat(i) >= cv(end), pValue(i) = alphaTable(end); end

    cValue(i) = interp1(alphaTable,cv,alpha(i),'linear','extrap');
end
h = pValue <= alpha;
end
