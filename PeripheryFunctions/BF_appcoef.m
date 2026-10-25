function a = BF_appcoef(C,L,wname,N)
% BF_appcoef   Extract/reconstruct level-N approximation (toolbox-free appcoef()).
% L = [len(cA_J), len(cD_J), ..., len(cD_1), len(x)]; N = 0..J (N = 0 gives the signal).
J = numel(L) - 2;
if nargin < 4 || isempty(N), N = J; end
if ~isnumeric(N) || ~isscalar(N) || N ~= fix(N) || N < 0 || N > J
    error('BF_appcoef:badLevel','Level %g outside 0..%d.',N,J);
end
a = C(1:L(1));
for k = J:-1:N+1                       % step back down to level N
    d = BF_detcoef(C,L,k);
    % (HCTSA2 2.2.1: the target is the length of cA_(k-1), as in BF_wrcoef)
    a = BF_idwt(a(:),d(:),wname,L(numel(L)-k+1));
end
a = a(:)';                             % a row, as C (and BF_detcoef)
end
