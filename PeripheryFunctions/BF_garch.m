function Mdl = BF_garch(P,Q)
% BF_garch   GARCH(P,Q) conditional variance model (toolbox-free garch()).
%
% Variance equation:
%   sigma2_t = Constant + sum_i GARCH{i}*sigma2_{t-i} + sum_j ARCH{j}*e_{t-j}^2
%   e_t      = y_t - Offset
%
% P is the GARCH degree (lagged variances), Q the ARCH degree (lagged squared
% innovations), matching MATLAB's ordering.
%
% Unestimated coefficients are NaN, as in the Econometrics Toolbox. Fit with
% BF_garchEstimate and filter with BF_garchInfer.

if nargin < 1 || isempty(P), P = 0; end
if nargin < 2 || isempty(Q), Q = 0; end

Mdl.P = P;
Mdl.Q = Q;
Mdl.Constant = NaN;
Mdl.GARCH = repmat({NaN},1,P);
Mdl.ARCH  = repmat({NaN},1,Q);
Mdl.Offset = 0;
Mdl.Distribution = 'Gaussian';
end
