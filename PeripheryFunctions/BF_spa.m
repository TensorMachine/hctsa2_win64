function g = BF_spa(y,M,nFreq)
% BF_spa   Smoothed spectral estimate (toolbox-free System Identification spa()).
%
% Blackman-Tukey estimate: the sample covariance function is tapered with a
% symmetric Hann lag window of length 2M-1 and Fourier transformed.
%
%   R(k)   = (1/N) sum_t (y_t - ybar)(y_{t-k} - ybar)
%   w(k)   = 0.5*(1 + cos(pi*k/(M-1))),   |k| <= M-1
%   Phi(f) = sum_k w(k) R(k) cos(f*k)
%
% Defaults verified against reference output from a licensed MATLAB:
%   M = min(30, floor(N/10)) and 128 frequencies at (1:128)/128*pi.
% Note the window denominator is M-1, not M -- using M is visibly wrong (~0.6%
% error in the normalised band powers).
%
% Returns a struct with the fields hctsa reads from an idfrd object:
% .frequency and .Spectrumdata.

y = y(:);
y = y(isfinite(y));
N = numel(y);
if nargin < 2 || isempty(M),     M = min(30,floor(N/10)); end
if nargin < 3 || isempty(nFreq), nFreq = 128; end
M = max(M,2);

e = y - mean(y);
L = M - 1;
R = zeros(2*L+1,1);
for k = 0:L
    c = sum(e(1+k:end).*e(1:end-k))/N;
    R(L+1+k) = c;
    R(L+1-k) = c;
end

tau = (-L:L)';
w = 0.5*(1 + cos(pi*tau/L));

f = ((1:nFreq)/nFreq*pi)';
S = zeros(nFreq,1);
for i = 1:nFreq
    S(i) = sum(w.*R.*cos(f(i)*tau));
end

g.frequency = f;
g.Spectrumdata = S;
end
