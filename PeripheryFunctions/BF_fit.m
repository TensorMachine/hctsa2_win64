function [cfun,gof,output] = BF_fit(x,y,model,varargin)
% BF_fit   Toolbox-free replacement for the Curve Fitting Toolbox's fit().
%
% Covers every model form used in hctsa. Linear-in-parameter models are solved in
% closed form; the rest use Levenberg-Marquardt with analytic Jacobians.
%
%---INPUTS:
% x, y    data (vectors, any orientation)
% model   either a built-in name ('poly2','sin1','gauss1','exp1','power1',...)
%         or a struct from BF_fittype
%
%---OUTPUTS:
% cfun    struct of fitted coefficients, accessed as cfun.a, cfun.a1, cfun.p1, ...
%         (matches how the Curve Fitting Toolbox's cfit object is indexed).
%         Evaluate it with BF_fitEval(cfun,x), not feval.
% gof     struct: sse, rsquare, dfe, adjrsquare, rmse
% output  struct: numobs, numparam, residuals, exitflag, iterations
%
% Errors use the same identifiers the Curve Fitting Toolbox raises
% ('curvefit:fit:nanComputed', 'curvefit:fit:powerFcnsRequirePositiveData', ...)
% so existing hctsa try/catch blocks keep working unchanged.

x = x(:);
y = y(:);

startPoint = [];
if isstruct(model)
    if isfield(model,'options') && isfield(model.options,'StartPoint')
        startPoint = model.options.StartPoint;
    end
    modelName = model.model;
else
    modelName = model;
end

% Late StartPoint override, e.g. BF_fit(x,y,'exp1','StartPoint',[1 -1])
for k = 1:2:numel(varargin)-1
    if strcmpi(varargin{k},'StartPoint')
        startPoint = varargin{k+1};
    end
end

% Drop non-finite pairs, as fit() does
ok = isfinite(x) & isfinite(y);
x = x(ok);
y = y(ok);

m = BF_fitModel(modelName);
nP = numel(m.coeffs);
N  = numel(x);
if N < nP
    error('curvefit:fit:notEnoughData', ...
          'Not enough data (%u points) to fit %u coefficients.',N,nP);
end
if strcmpi(modelName,'power1') && any(x <= 0)
    error('curvefit:fit:powerFcnsRequirePositiveData', ...
          'Power functions cannot be fit to non-positive xdata.');
end

%-------------------------------------------------------------------------------
%% Solve
%-------------------------------------------------------------------------------
if m.linear
    A = m.basis(x);
    if ~all(isfinite(A(:)))
        error('curvefit:fit:nanComputed','NaN computed by model function.');
    end
    p = A\y;
    r = y - A*p;
    iters = 0;
    exitflag = 1;
else
    if isempty(startPoint)
        startPoint = m.start(x,y);
    elseif numel(startPoint) ~= nP
        error('curvefit:fit:startPointSize', ...
              'StartPoint has %u elements but the model has %u coefficients.', ...
              numel(startPoint),nP);
    end
    [p,r,iters,exitflag] = levenbergMarquardt(m.fun,m.jac,x,y,startPoint(:));
end

if any(~isfinite(r))
    error('curvefit:fit:nanComputed','NaN computed by model function.');
end

%-------------------------------------------------------------------------------
%% Package outputs
%-------------------------------------------------------------------------------
cfun = struct();
for k = 1:nP
    cfun.(m.coeffs{k}) = p(k);
end
% Retained so BF_fitEval can reconstruct the model; leading underscore-free names
% would collide with coefficient names like 'a' or 'p1'.
cfun.BF_modelName = modelName;
cfun.BF_coeffs    = p(:)';

sse = sum(r.^2);
sst = sum((y - mean(y)).^2);
dfe = N - nP;

gof.sse        = sse;
gof.rsquare    = 1 - sse/sst;
gof.dfe        = dfe;
gof.adjrsquare = 1 - (sse/dfe)/(sst/(N-1));
gof.rmse       = sqrt(sse/dfe);

output.numobs     = N;
output.numparam   = nP;
output.residuals  = r;
output.exitflag   = exitflag;
output.iterations = iters;

end

%-------------------------------------------------------------------------------
function [p,r,iter,exitflag] = levenbergMarquardt(fun,jac,x,y,p)
% Damped Gauss-Newton. Defaults mirror the Curve Fitting Toolbox
% (TolFun/TolX 1e-6, MaxIter 400).

tolFun  = 1e-6;
tolX    = 1e-6;
maxIter = 400;
lambda  = 1e-3;

r   = y - fun(p,x);
sse = sum(r.^2);
exitflag = 0;

for iter = 1:maxIter
    J = jac(p,x);
    if ~all(isfinite(J(:))) || ~isfinite(sse)
        exitflag = -1;
        break
    end

    JtJ = J'*J;
    Jtr = J'*r;
    scale = diag(JtJ);
    scale(scale == 0) = 1; % keep the damping term well-defined

    stepTaken = false;
    for attempt = 1:12 % raise damping until the step actually reduces sse
        A = JtJ + lambda*diag(scale);
        warnState = warning('off','all');
        delta = A\Jtr;
        warning(warnState);

        if ~all(isfinite(delta))
            lambda = lambda*10;
            continue
        end

        pNew = p + delta;
        rNew = y - fun(pNew,x);
        sseNew = sum(rNew.^2);

        if isfinite(sseNew) && sseNew < sse
            stepTaken = true;
            break
        end
        lambda = lambda*10;
    end

    if ~stepTaken
        exitflag = 2; % no downhill step found: at a local minimum
        break
    end

    relSSE = (sse - sseNew)/max(sse,realmin);
    relX   = norm(delta)/(tolX + norm(p));

    p   = pNew;
    r   = rNew;
    sse = sseNew;
    lambda = max(lambda/10,1e-12);

    if relSSE < tolFun || relX < tolX
        exitflag = 1;
        break
    end
end

if exitflag == 0
    exitflag = 0; % hit maxIter without converging
end

end
