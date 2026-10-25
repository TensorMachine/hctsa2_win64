function lev = BF_wmaxlev(n,wname)
% BF_wmaxlev   Toolbox-free replacement for wmaxlev().
% MATLAB: lev = fix(log2(lx/(lf-1))), lf the filter length.
% n is the signal length or its size ([1 N] or [N 1]; [R C] gives min(R,C)).
if numel(n) == 2
    if min(n) == 1, n = max(n); else, n = min(n); end
elseif numel(n) ~= 1
    error('BF_wmaxlev:badSize','The first argument must be a length or a 2-element size.');
end
Lo_D = BF_wfilters(wname);
lev = max(0,fix(log2(n/(numel(Lo_D)-1))));
end
