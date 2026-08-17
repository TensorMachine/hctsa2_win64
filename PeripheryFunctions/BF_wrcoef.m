function y = BF_wrcoef(type,C,L,wname,N)
% BF_wrcoef   Reconstruct a single approximation or detail branch (wrcoef()).
%
% type: 'a' for approximation, 'd' for detail; N is the level.
J = numel(L) - 2;
lx = L(end);
if strcmpi(type,'a')
    a = BF_appcoef(C,L,wname,N);
    d = zeros(size(a));
else
    d = BF_detcoef(C,L,N);
    a = zeros(size(d));
end
% Cascade back up to full length, zeroing all other branches
for k = N:-1:1
    lenNext = L(numel(L)-k+1);         % length of the signal one level up
    a = BF_idwt(a(:),d(:),wname,lenNext);
    d = zeros(size(a));
end
y = a(1:lx);
y = y(:)';
end
