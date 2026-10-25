function v = BF_fitEval(cfun,x)
% BF_fitEval   Evaluate a BF_fit result, replacing feval() on a cfit object.

m = BF_fitModel(cfun.BF_modelName);
p = cfun.BF_coeffs;
if m.linear
    if isfield(cfun,'BF_scale')      % polyN, N >= 5, or 'Normalize','on': fitted in (x - mean)/std
        x = (x - cfun.BF_scale(1))/cfun.BF_scale(2);
    end
    v = m.basis(x)*p(:);
else
    v = m.fun(p,x);
end
end
