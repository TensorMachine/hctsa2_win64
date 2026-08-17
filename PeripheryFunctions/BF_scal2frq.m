function f = BF_scal2frq(a,wname,delta)
% BF_scal2frq   Scale to pseudo-frequency (toolbox-free scal2frq()).
%   F = centfrq(wname) / (a * delta)
if nargin < 3 || isempty(delta), delta = 1; end
f = BF_centfrq(wname)./(a(:)'*delta);
end
