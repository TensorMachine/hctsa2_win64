function [transdat,lambda] = BF_boxcox(y)
% BF_boxcox   Toolbox-free replacement for the Financial Toolbox's boxcox().
%
% Box-Cox transform with the exponent chosen by maximum likelihood:
%   lambda ~= 0 :  (y^lambda - 1)/lambda
%   lambda == 0 :  log(y)
%
% Requires strictly positive data, as boxcox() does.
%
% NOTE: no master operation in the default hctsa library requests the 'boxcox'
% detrending in PP_Compare, so this exists to remove the last hard-error path
% rather than to reproduce any library feature. The lambda search is a plain
% profile-likelihood maximisation over [-2,2]; MathWorks does not document their
% search, so lambda may differ marginally from theirs.

y = y(:);
if any(y <= 0)
    error('finance:boxcox:nonPositiveData','Box-Cox requires positive data.');
end

n = numel(y);
sumLogY = sum(log(y));

% Profile log-likelihood, negated for minimisation
negLL = @(lam) 0.5*n*log(var(bcTransform(y,lam),1)) - (lam-1)*sumLogY;

lambda = fminbnd(negLL,-2,2,optimset('TolX',1e-6));
transdat = bcTransform(y,lambda);

end

function z = bcTransform(y,lambda)
if abs(lambda) < eps
    z = log(y);
else
    z = (y.^lambda - 1)/lambda;
end
end
