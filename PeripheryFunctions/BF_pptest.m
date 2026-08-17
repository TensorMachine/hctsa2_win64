function [h,pValue,stat,cValue,reg] = BF_pptest(y,varargin)
% BF_pptest   Phillips-Perron unit root test (toolbox-free pptest()).
%
% Null hypothesis: a unit root. The Dickey-Fuller t statistic is corrected
% non-parametrically for serial correlation in the residuals:
%
%   Z_t = sqrt(g0/lam2)*t_a - (lam2-g0)/(2*lam) * (n*se(a)/s)
%
% where g0 is the residual variance, lam2 the Newey-West long-run variance with
% Bartlett weights, and t_a the ordinary t statistic for a = 1.
%
% Name/value pairs:
%   'lags'  Newey-West lags (default 0)
%   'model' 'ar' (no deterministic terms), 'ard' (drift), 'ts' (trend). Default 'ar'.
%   'test'  't1' (default) or 't2'
%   'alpha' significance level (default 0.05)
%
% p-values are interpolated from Dickey-Fuller tables and clipped to
% [0.001, 0.999], matching MATLAB (which warns econ:pptest:LeftTailStatTooSmall
% at the lower boundary).

y = y(:);
y = y(isfinite(y));

lags = 0; model = 'ar'; testStat = 't1'; alpha = 0.05;
for k = 1:2:numel(varargin)-1
    switch lower(varargin{k})
    case 'lags',  lags     = varargin{k+1};
    case 'model', model    = varargin{k+1};
    case 'test',  testStat = varargin{k+1};
    case 'alpha', alpha    = varargin{k+1};
    end
end
lags = lags(:)'; alpha = alpha(:)';
n = max(numel(lags),numel(alpha));
if numel(lags)==1,  lags  = repmat(lags,1,n);  end
if numel(alpha)==1, alpha = repmat(alpha,1,n); end
if ischar(model),    model    = {model};    end
if ischar(testStat), testStat = {testStat}; end
if numel(model)==1,    model    = repmat(model,1,n);    end
if numel(testStat)==1, testStat = repmat(testStat,1,n); end

% Dickey-Fuller critical values (asymptotic), by model
alphaTable = [0.001 0.010 0.025 0.050 0.100 0.900 0.950 0.975 0.990 0.999];
switch lower(model{1})
case 'ar',  cvTable = [-3.43 -2.58 -2.23 -1.95 -1.62 0.89 1.28 1.62 2.00 2.77];
case 'ard', cvTable = [-4.38 -3.43 -3.12 -2.86 -2.57 -0.44 -0.07 0.23 0.60 1.30];
otherwise,  cvTable = [-4.90 -3.96 -3.66 -3.41 -3.13 -1.25 -0.94 -0.66 -0.33 0.34];
end

yLag = y(1:end-1);
dy   = y(2:end);
nObs = numel(dy);

switch lower(model{1})
case 'ar',  X = yLag;
case 'ard', X = [ones(nObs,1), yLag];
otherwise,  X = [ones(nObs,1), (1:nObs)', yLag];
end
k = size(X,2);

b   = X\dy;
e   = dy - X*b;
SSE = sum(e.^2);
s2  = SSE/(nObs-k);
XtXi = inv(X'*X);
seB = sqrt(s2*diag(XtXi));

aIdx = k;                       % the lagged-level coefficient is last
aHat = b(aIdx);
seA  = seB(aIdx);
tA   = (aHat-1)/seA;

% Regression diagnostics, matching MATLAB's econometrics conventions
LL  = -0.5*nObs*(1 + log(2*pi) + log(SSE/nObs));
reg.coeff = b(:)';
reg.LL    = LL;
reg.AIC   = -2*LL + 2*k;
reg.BIC   = -2*LL + k*log(nObs);
reg.HQC   = -2*LL + 2*k*log(log(nObs));
reg.RMSE  = sqrt(s2);

g0 = SSE/nObs;
s  = sqrt(s2);

stat = zeros(1,n); pValue = zeros(1,n); cValue = zeros(1,n);
for i = 1:n
    lam2 = g0;
    for j = 1:lags(i)
        w = 1 - j/(lags(i)+1);
        lam2 = lam2 + 2*w*sum(e(j+1:end).*e(1:end-j))/nObs;
    end
    lam = sqrt(lam2);

    if strcmpi(testStat{i},'t2')
        stat(i) = nObs*(aHat-1)*sqrt(g0/lam2) - 0.5*(lam2-g0)*nObs^2*seA^2/s2;
    else
        stat(i) = sqrt(g0/lam2)*tA - (lam2-g0)/(2*lam)*(nObs*seA/s);
    end

    pValue(i) = interp1(cvTable,alphaTable,stat(i),'linear');
    if stat(i) <= cvTable(1),   pValue(i) = alphaTable(1);   end
    if stat(i) >= cvTable(end), pValue(i) = alphaTable(end); end
    cValue(i) = interp1(alphaTable,cvTable,alpha(i),'linear','extrap');
end
h = pValue <= alpha;
end
