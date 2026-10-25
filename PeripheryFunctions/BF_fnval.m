function v = BF_fnval(sp,x)
% BF_fnval   Toolbox-free replacement for fnval(), for splines from BF_spap2.
%
% Returns values in the same shape as X, matching fnval's behaviour. (HCTSA2 2.2.1:
% a matrix X used to give a column.) Outside the knot range the end pieces are
% extrapolated, on both sides (BF_bsplineBasis).

v = reshape(BF_bsplineBasis(sp.knots,sp.order,x(:)) * sp.coefs(:),size(x));
end
