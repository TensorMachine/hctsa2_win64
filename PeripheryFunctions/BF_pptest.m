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
% The 't2' statistic is the normalised bias, corrected the same way:
%
%   Z_rho = n*(a-1) - (n*se(a))^2/(2*s^2) * (lam2-g0)
%
% Name/value pairs:
%   'lags'  Newey-West lags (default 0), integers in [0, T-1]; may be a vector
%   'model' 'ar' (no deterministic terms), 'ard' (drift), 'ts' (trend). Default 'ar'.
%   'test'  't1' (default) or 't2'
%   'alpha' significance level (default 0.05), in [0.001, 0.999]; may be a vector
% 'lags' and 'alpha' may be vectors of the same length (one test per element), but
% unlike pptest() this function takes a single 'model' and a single 'test': a
% vector of either raises BF_pptest:vectorModel rather than being collapsed.
%
% Critical values and p-values come from finite-sample tables of the statistics'
% null distributions (BF_pptestTable: 15 effective sample sizes T from 10 to 10000,
% 50 significance levels), interpolated linearly in T and then in the statistic,
% the method pptest documents. The tables are HCTSA2's own Monte Carlo simulation
% of pptest's documented design (docs/hctsa2_2.1.1/tools/make_pptest_tables.py), so
% p-values agree with MATLAB's to within the two simulations' sampling error
% (about 1e-4). p-values are clipped to [0.001, 0.999], as MATLAB does (it warns
% econ:pptest:LeftTailStatTooSmall at the lower boundary). As in pptest, T must be
% at least 10, and the log-likelihood is Gaussian with the residual variance
% SSE/(T-k).

y = y(:);
y = y(isfinite(y));

lags = 0; model = 'ar'; testStat = 't1'; alpha = 0.05;
if mod(numel(varargin),2) ~= 0
    error('BF_pptest:badOptions','Options must be name/value pairs.');
end
for k = 1:2:numel(varargin)-1
    switch lower(varargin{k})
    case 'lags',  lags     = varargin{k+1};
    case 'model', model    = varargin{k+1};
    case 'test',  testStat = varargin{k+1};
    case 'alpha', alpha    = varargin{k+1};
    otherwise
        error('BF_pptest:unknownOption','Unknown option ''%s''.',char(varargin{k}));
    end
end
lags = lags(:)'; alpha = alpha(:)';
n = max(numel(lags),numel(alpha));
if (numel(lags) ~= 1 && numel(lags) ~= n) || (numel(alpha) ~= 1 && numel(alpha) ~= n)
    error('BF_pptest:badOptions','''lags'' and ''alpha'' must be scalars or have the same length.');
end
if numel(lags)==1,  lags  = repmat(lags,1,n);  end
if numel(alpha)==1, alpha = repmat(alpha,1,n); end
if ischar(model),    model    = {model};    end
if ischar(testStat), testStat = {testStat}; end
% (HCTSA2 2.2.1: one regression and one critical-value table serve every test,
% so a vector of models or tests, which pptest() runs one per element, is refused
% rather than silently run with the first model.)
if numel(model) ~= 1 || numel(testStat) ~= 1
    error('BF_pptest:vectorModel','''model'' and ''test'' must each be a single name.');
end
model    = repmat(model,1,n);
testStat = repmat(testStat,1,n);

yLag = y(1:end-1);
dy   = y(2:end);
nObs = numel(dy);                % the effective sample size T
if nObs < 10
    error('econ:pptest:EffectiveSampleSizeLessThanTabulatedValues', ...
          'Effective sample size must be at least 10 (the smallest tabulated).');
end
if any(alpha < 0.001 | alpha > 0.999)
    error('econ:pptest:AlphaOutOfRange','Alpha must be between 0.001 and 0.999.');
end
if any(lags < 0 | lags ~= fix(lags) | lags >= nObs)
    error('econ:pptest:LagsOutOfRange', ...
          'Lags must be integers between 0 and the effective sample size minus 1.');
end

switch lower(model{1})
case 'ar',  X = yLag;
case 'ard', X = [ones(nObs,1), yLag];
otherwise,  X = [ones(nObs,1), (1:nObs)', yLag];
end
k = size(X,2);
% Constant data (or any rank-deficient regressors) have no unit-root test
if all(y == y(1)) || rcond(X'*X) < eps
    error('econ:pptest:RankDeficient', ...
          'The regressors are rank deficient (is the series constant?).');
end

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

% Regression diagnostics, with pptest's conventions: the Gaussian log-likelihood
% of the residuals at the residual variance s2 = SSE/(T-k)
LL  = -0.5*nObs*log(2*pi*s2) - 0.5*(nObs-k);
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
        stat(i) = nObs*(aHat-1) - 0.5*(lam2-g0)*nObs^2*seA^2/s2;
    else
        stat(i) = sqrt(g0/lam2)*tA - (lam2-g0)/(2*lam)*(nObs*seA/s);
    end

    [cvRow,sigLevels] = criticalValues(model{1},testStat{i},nObs);
    if stat(i) <= cvRow(1)
        pValue(i) = sigLevels(1);
    elseif stat(i) >= cvRow(end)
        pValue(i) = sigLevels(end);
    else
        pValue(i) = interp1(cvRow,sigLevels,stat(i),'linear');
    end
    cValue(i) = interp1(sigLevels,cvRow,alpha(i),'linear');
end
h = stat < cValue;               % a left-tailed test
end

%-------------------------------------------------------------------------------
function [row,sigLevels] = criticalValues(model,test,T)
% The table's critical values at effective sample size T, linear in T between
% the tabulated sizes; the largest tabulated size (10000) stands for all T above.
[cv,sizes,sigLevels] = BF_pptestTable(model,test);
if T >= sizes(end)
    row = cv(end,:);
    return
end
i = find(sizes <= T,1,'last');
w = (T - sizes(i))/(sizes(i+1) - sizes(i));
row = (1 - w)*cv(i,:) + w*cv(i+1,:);
end
