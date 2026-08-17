function [V,logL] = BF_garchInfer(Mdl,y)
% BF_garchInfer   Infer conditional variances (toolbox-free infer()).
%
% Returns the conditional VARIANCES V (not standard deviations), matching
% MATLAB, and the Gaussian log-likelihood.

K = Mdl.Constant;
G = cell2mat(Mdl.GARCH);
A = cell2mat(Mdl.ARCH);
if isempty(G), G = []; end
if isempty(A), A = []; end
[V,e] = BF_garchVar(y,K,G,A,Mdl.Offset);
logL = -0.5*sum(log(2*pi) + log(V) + e.^2./V);
end
