function lev = BF_wmaxlev(n,wname)
% BF_wmaxlev   Toolbox-free replacement for wmaxlev().
% MATLAB: lev = fix(log2(lx/(lf-1))), lf the filter length.
if numel(n) > 1, n = numel(n); end
Lo_D = BF_wfilters(wname);
lev = max(0,fix(log2(n/(numel(Lo_D)-1))));
end
