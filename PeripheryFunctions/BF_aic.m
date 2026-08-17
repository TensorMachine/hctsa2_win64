function v = BF_aic(m)
% BF_aic   Akaike's Information Criterion for an identified model.
%
% Toolbox-free replacement for the System Identification aic().
%
%   AIC = log(V) + 2*d/N
%
% with V the loss function (mean squared prediction error), d the number of
% estimated parameters and N the number of data samples.
%
% Verified against reference output from a licensed MATLAB: this is numerically
% equal to log(FPE) to ~3e-9, consistent with FPE = V*(N+d)/(N-d), which is how
% hctsa's own comment describes it ("~ log(fpe)").

if ~isfield(m,'EstimationInfo') || ~isfield(m.EstimationInfo,'LossFcn')
    error('BF_aic:noLossFcn','Model has no EstimationInfo.LossFcn.');
end
V = m.EstimationInfo.LossFcn;
d = m.np;
N = m.N;
v = log(V) + 2*d/N;
end
