function [phi,psi,xval] = BF_wavefun(wname,iter)
% BF_wavefun   Scaling and wavelet functions (wavefun()).
%
% Toolbox-free replacement, written against MathWorks' wavefun.m source for the
% orthogonal case (wtype 1), which is the only case hctsa uses.
%
% The construction is a delta-seeded cascade via `upcoef`, NOT an eigenvector of
% the refinement operator:
%
%   coef = sqrt(2)^iter
%   phi  = coef * upcoef('a',1,Lo_R,'dummy',iter)
%   psi  = coef * upcoef('d',1,Lo_R,Hi_R,iter)
%
% then padded as [0, wkeep1(.,nb), zeros(1,1+dn)] with nb and dn from getNBpts.
% `upcoef` reconstructs by repeatedly upsampling by 2 and convolving: the 'd'
% form uses Hi_R on the first step and Lo_R thereafter.
%
% Symlets, coiflets and dmeyer additionally have psi negated -- a sign
% convention with no mathematical content but which the reference output
% follows.

if nargin < 2 || isempty(iter), iter = 8; end

[~,~,Lo_R,Hi_R] = BF_wfilters(wname);
Lo_R = Lo_R(:)';
Hi_R = Hi_R(:)';

coef = sqrt(2)^iter;
pas  = 1/(2^iter);
long = numel(Lo_R);
nbpts = (long-1)/pas + 1;

phi = coef*upcoefLocal(Lo_R,Lo_R,iter);
psi = coef*upcoefLocal(Hi_R,Lo_R,iter);

[nbpts,nb,dn] = getNBpts(nbpts,iter,long);
phi = [0, wkeep1Local(phi,nb), zeros(1,1+dn)];
psi = [0, wkeep1Local(psi,nb), zeros(1,1+dn)];

% Sign convention: coiflet, symlet, dmeyer
debut = lower(wname(1:min(2,end)));
if any(strcmp(debut,{'co','sy','dm'}))
    psi = -psi;
end

xval = linspace(0,(nbpts-1)*pas,nbpts);
end

%-------------------------------------------------------------------------------
function y = upcoefLocal(fFirst,fRest,iter)
% Reconstruct from a single unit coefficient: upsample by 2 and convolve, using
% fFirst on the first step and fRest on the rest.
y = 1;
for k = 1:iter
    n = numel(y);
    yu = zeros(1,2*n-1);
    yu(1:2:end) = y;              % dyadup
    if k == 1, f = fFirst; else f = fRest; end
    y = conv(yu,f);
end
end

%-------------------------------------------------------------------------------
function y = wkeep1Local(x,len)
% Keep the central len samples.
n = numel(x);
if len >= n
    y = x;
    return
end
d = (n-len)/2;
first = floor(d) + 1;
y = x(first:first+len-1);
end

%-------------------------------------------------------------------------------
function [nbpts,nb,dn] = getNBpts(nbpts,iter,long)
lplus = long - 2;
nb = 1;
for kk = 1:iter, nb = 2*nb + lplus; end
dn = nbpts - nb - 2;
if dn < 0
    nbpts = nbpts - dn;
    dn = 0;
end
end
