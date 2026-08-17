function a = BF_appcoef(C,L,wname,N)
% BF_appcoef   Extract/reconstruct level-N approximation (toolbox-free appcoef()).
J = numel(L) - 2;
if nargin < 4 || isempty(N), N = J; end
a = C(1:L(1));
for k = J:-1:N+1                       % step back down to level N
    d = BF_detcoef(C,L,k);
    a = BF_idwt(a(:),d(:),wname,L(numel(L)-k));
end
end
