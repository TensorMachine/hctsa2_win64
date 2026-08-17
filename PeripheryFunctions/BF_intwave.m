function [psiInt,xval] = BF_intwave(wname,prec)
% BF_intwave   Integral of the wavelet function (toolbox-free intwave()).
if nargin < 2 || isempty(prec), prec = 10; end
[~,psi,xval] = BF_wavefun(wname,prec);
dx = xval(2) - xval(1);
psiInt = cumsum(psi)*dx;
end
