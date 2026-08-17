function s = BF_wnoisest(C,L,S)
% BF_wnoisest   Toolbox-free replacement for wnoisest().
%
% Robust estimate of the noise standard deviation from the level-S detail
% coefficients: median(|cD|)/0.6745, the standard MAD estimator for Gaussian
% noise that MATLAB uses.

if nargin < 3 || isempty(S), S = 1; end
s = zeros(size(S));
for k = 1:numel(S)
    d = BF_detcoef(C,L,S(k));
    s(k) = median(abs(d(:)))/0.6745;
end
end
