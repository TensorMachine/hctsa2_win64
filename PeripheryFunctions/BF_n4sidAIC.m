function n = BF_n4sidAIC(y,nmax,maxsize)
% BF_n4sidAIC   AIC order selection for the n4sid past horizon.
%
% Direct port of the local function `localAIC` in MATLAB's
% toolbox/ident/ident/+idpack/@ssdata/n4sid_time.m, for the single-output,
% no-input case (ny = 1, nu = 0, so nz = 1).
%
% n4sid picks its past horizon (`sy` in N4Horizon = [r sy su]) as the AIC-optimal
% ARX order, capped at 4*modelOrder:
%
%   radef = ceil(1.5*order);
%   maxo  = ceil(min([4*order, ...]));
%   auxord = localAIC(z1,maxo,ny,maxsize);
%   n4h = [radef, auxord, auxord];
%
% The QR is accumulated blockwise so the regression matrix never exceeds
% MaxSize entries; the diagonal of the triangular factor then gives the
% residual variance at every order at once, since column k+1 of the lagged
% data matrix regressed on columns 1..k is exactly the AR(k) fit.

if nargin < 3 || isempty(maxsize), maxsize = 250e3; end

y = y(:);
N = numel(y);
nz = 1;                       % ny + nu
p  = 1;                       % ny

M = max(floor(maxsize/nz/nmax),nmax+1);
R1 = zeros(0,nmax*nz);
nr = 1;
for k = nmax:M:N-1
    jj = (k : min(N,k+M-1))';
    phi = zeros(numel(jj),nmax*nz);
    for kz = 1:nmax
        phi(:,kz) = y(jj-nmax+kz);
    end
    R1 = triu(qr([R1;phi]));
    [nr,nrc] = size(R1);
    R1 = R1(1:min(nr,nrc),:);
end

Neff = N - nmax + 1;
m = min(nmax-1,floor((Neff-p-1)/nz));
V = zeros(1,m+2);
V(1) = Inf;
for k = 0:m
    blk = R1(k*nz+1:k*nz+p, k*nz+1:k*nz+p)/nr;
    V(k+1) = sum(log(diag(blk).^2)) + 2*nz*p*k/Neff;
end

[~,n] = min(V);
n = n - 1;
end
