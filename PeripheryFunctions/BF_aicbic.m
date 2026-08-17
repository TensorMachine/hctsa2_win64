function [aic,bic] = BF_aicbic(logL,numParam,numObs)
% BF_aicbic   Information criteria (toolbox-free aicbic()).
%   AIC = -2*logL + 2*numParam
%   BIC = -2*logL + numParam*log(numObs)
logL = logL(:); numParam = numParam(:);
aic = -2*logL + 2*numParam;
if nargin > 2
    bic = -2*logL + numParam.*log(numObs(:));
else
    bic = [];
end
end
