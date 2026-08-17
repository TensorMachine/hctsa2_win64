function v = BF_fitEval(cfun,x)
% BF_fitEval   Evaluate a BF_fit result, replacing feval() on a cfit object.

m = BF_fitModel(cfun.BF_modelName);
p = cfun.BF_coeffs;
if m.linear
    v = m.basis(x)*p(:);
else
    v = m.fun(p,x);
end
end
