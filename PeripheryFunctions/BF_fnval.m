function v = BF_fnval(sp,x)
% BF_fnval   Toolbox-free replacement for fnval(), for splines from BF_spap2.
%
% Returns values in the same orientation as X, matching fnval's behaviour.

isRow = (size(x,1) == 1 && size(x,2) > 1);
v = BF_bsplineBasis(sp.knots,sp.order,x(:)) * sp.coefs(:);
if isRow
    v = v';
end
end
