function f = BF_centfrq(wname,prec)
% BF_centfrq   Centre frequency of a wavelet (toolbox-free centfrq()).
%
% The dominant frequency of the wavelet function, found from the peak of its
% FFT magnitude. Reference values: db3 -> 0.8, sym2 -> 0.6667.
if nargin < 2 || isempty(prec), prec = 8; end
[~,psi,xval] = BF_wavefun(wname,prec);
T = xval(end) - xval(1);
n = numel(psi);
psi = psi - mean(psi);
sp = abs(fft(psi));
sp = sp(1:floor(n/2));
[~,idx] = max(sp);
f = (idx-1)/T;
end
