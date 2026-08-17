function e = BF_pe(m,y)
% BF_pe   Prediction errors of an AR model (toolbox-free pe()).
%
%   e(t) = A(q) y(t)
%
% MATLAB's pe() defaults to InitialCondition = 'e', i.e. it estimates the
% presample values rather than assuming zeros. For an AR(n) model those n free
% presample values enter only the first n prediction errors, and choosing them
% to minimise the squared error sets those errors exactly to zero. So the
% estimated-initial-condition result is simply the filtered series with its
% startup transient zeroed.
%
% This matters: leaving the transient in inflates mean(e.^2) by roughly n/T and
% put PP_ModelFit's ratios 4.5e-4 off the reference. Zeroing it brings them to
% ~1e-7.

y = y(:);
e = filter(m.A,1,y);
n = numel(m.A) - 1;
if n > 0 && numel(e) > n
    e(1:n) = 0;
end
end
