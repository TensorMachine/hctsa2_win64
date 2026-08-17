function coefs = BF_cwt(x,scales,wname,prec)
% BF_cwt   Continuous wavelet transform (toolbox-free legacy cwt()).
%
% Reproduces MATLAB's cwt(x,scales,wname) for real, non-analytic wavelets: the
% integrated wavelet is sampled at each scale, convolved with the signal, and
% differenced.

if nargin < 4 || isempty(prec), prec = 10; end
x = x(:)';
n = numel(x);
[psiInt,xval] = BF_intwave(wname,prec);
xval = xval - xval(1);
dx = xval(2) - xval(1);
xmax = xval(end);

scales = scales(:)';
coefs = zeros(numel(scales),n);
for k = 1:numel(scales)
    a = scales(k);
    j = 1 + floor((0:a*xmax)/(a*dx));
    if numel(j) == 1, j = [1 1]; end
    f = psiInt(j);
    f = f(end:-1:1);
    z = diff(conv(x,f));
    coefs(k,:) = -sqrt(a)*wkeepCentre(z,n);
end
end

function y = wkeepCentre(z,n)
s = floor((numel(z)-n)/2);
y = z(s+1:s+n);
end
