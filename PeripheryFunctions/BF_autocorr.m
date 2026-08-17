function [acf,lags,bounds] = BF_autocorr(y,varargin)
% BF_autocorr   Sample autocorrelation function (toolbox-free autocorr()).
%
% Accepts NumLags as a name/value pair, in either the 'NumLags',L or NumLags=L
% form. acf(1) is lag 0 and equals 1. bounds are the +/-2/sqrt(N) significance
% lines for a white-noise null.

y = y(:); y = y(isfinite(y));
N = numel(y);
numLags = min(20,N-1);
for k = 1:2:numel(varargin)-1
    if strcmpi(varargin{k},'NumLags'), numLags = varargin{k+1}; end
end
r = y - mean(y);
c0 = sum(r.^2);
acf = zeros(numLags+1,1);
acf(1) = 1;
for k = 1:numLags
    acf(k+1) = sum(r(1+k:end).*r(1:end-k))/c0;
end
lags = (0:numLags)';
bounds = [2;-2]/sqrt(N);
end
