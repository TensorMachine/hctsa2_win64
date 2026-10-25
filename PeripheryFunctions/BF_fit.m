function [cfun,gof,output] = BF_fit(x,y,model,varargin)
% BF_fit   Toolbox-free replacement for the Curve Fitting Toolbox's fit().
%
% Covers every model form used in hctsa. Linear-in-parameter models are solved in
% closed form. Nonlinear models with Lower/Upper bounds (SP_Specparam) use a
% trust-region-reflective method from the given start, as fit() does.
%
% Unbounded nonlinear models report the best fit found by several solvers from
% one or two start points (HCTSA2 2.1.1):
%   starts   the caller's StartPoint; or else fit()'s default start rule for the
%            model and, for exp1, power1, gaussN and sinN, HCTSA2 2.0's rule too
%            (BF_fitModel's start and start2)
%   solvers  from each start: fit()'s own method -- variable projection for the
%            models fit() solves separably (exp1, power1, gaussN), a
%            trust-region-reflective method on all coefficients otherwise --
%            stopped by fit()'s tolerances in the form MATLAB documents them; the
%            same continued to convergence from there; the trust-region method on
%            all coefficients (separable models); and Levenberg-Marquardt
%            (HCTSA2 2.0's solver)
%   choice   the first candidate (fit()'s method from the first start) unless
%            another reaches a sum of squares lower by more than 1e-6 (relative)
% So a fit is never worse than HCTSA2 2.0's (beyond that 1e-6), and where fit() stops
% early or in a poorer local minimum, BF_fit reports the better fit. Replayed on
% all 3,137 fit() calls hctsa's library makes on the 20 reference series
% (licensed MATLAB R2025b; docs/hctsa2_2.1.1), BF_fit's sum of squares is never
% more than 1e-5 (relative) above fit()'s, and lower by more than 1e-6 on 775.
% gaussN and sinN coefficients follow fit()'s signs (c >= 0, b >= 0).
%
%---INPUTS:
% x, y    data (vectors, any orientation)
% model   either a built-in name ('poly2','sin1','gauss1','exp1','power1',...)
%         or a struct from BF_fittype
% Name/value pairs: 'StartPoint', 'Lower', 'Upper', 'Normalize'. (HCTSA2 2.2.1: any
% other name is an error. So is a fitoptions property BF_fit does not implement, such
% as Weights, Exclude, Robust or TolFun, and finite Lower/Upper bounds on a linear
% model; MaxIter and Display are accepted and ignored.)
%
% 'Normalize','on' (HCTSA2 3.0.0, for upstream's PP_Compare) solves a linear model
% (polyN, ...) in the variable (x - mean(x))/std(x), as fit() does, and reports the
% coefficients in that variable, as fit() does; BF_fitEval evaluates in it. It is
% implemented for linear models only: on a nonlinear model it is an error.
%
% polyN with N >= 5 is solved in the variable (x - mean(x))/std(x), as polyfit
% does: on the raw Vandermonde matrix the solve loses the optimum (the condition
% number is ~1e22 for degree 6 at N = 4097). The coefficients p1..p(N+1) are
% still reported for the raw x (with 'Normalize','on': for the normalised variable,
% above); cfun.BF_scale holds mean and std for BF_fitEval.
% poly1 to poly4 are solved on the raw x, as before.
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
lb = []; ub = [];
normalize = false;
if isstruct(model)
    if isfield(model,'options') && isfield(model.options,'StartPoint')
        startPoint = model.options.StartPoint;
    end
    if isfield(model,'options') && isfield(model.options,'Lower')
        lb = model.options.Lower;
    end
    if isfield(model,'options') && isfield(model.options,'Upper')
        ub = model.options.Upper;
    end
    if isfield(model,'options')
        checkUnsupported(model.options);
        if isfield(model.options,'Normalize')
            normalize = onOff(model.options.Normalize);
        end
    end
    modelName = model.model;
else
    modelName = model;
end
if isstring(modelName), modelName = char(modelName); end   % "poly2" as well as 'poly2'

% Late overrides, e.g. BF_fit(x,y,'exp1','StartPoint',[1 -1])
if mod(numel(varargin),2) ~= 0
    error('BF_fit:badOptions','Options must be name/value pairs.');
end
for k = 1:2:numel(varargin)-1
    switch lower(varargin{k})
        case 'startpoint', startPoint = varargin{k+1};
        case 'lower',      lb = varargin{k+1};
        case 'upper',      ub = varargin{k+1};
        case 'normalize',  normalize = onOff(varargin{k+1});
        case {'maxiter','display'}   % accepted and ignored, as in BF_fitoptions
        otherwise
            error('BF_fit:unknownOption','Unknown or unsupported option ''%s''.',char(varargin{k}));
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
if m.linear && (any(isfinite(lb(:))) || any(isfinite(ub(:))))
    error('BF_fit:unsupportedOption','Lower/Upper bounds are not supported for linear models.');
end
if normalize && ~m.linear
    error('BF_fit:unsupportedOption', ...
          'The fit option ''Normalize'' is implemented for linear models only.');
end

%-------------------------------------------------------------------------------
%% Solve
%-------------------------------------------------------------------------------
% polyN, N >= 5, or 'Normalize','on'
scaled = m.linear && ((isfield(m,'centre') && m.centre) || normalize);
if m.linear
    if scaled
        mu = mean(x); sg = std(x);
        if ~(sg > 0), sg = 1; end
        A = m.basis((x - mu)/sg);
    else
        A = m.basis(x);
    end
    if ~all(isfinite(A(:)))
        error('curvefit:fit:nanComputed','NaN computed by model function.');
    end
    p = A\y;
    r = y - A*p;
    iters = 0;
    exitflag = 1;
else
    userStart = ~isempty(startPoint);
    if ~userStart
        startPoint = m.start(x,y);
    elseif numel(startPoint) ~= nP
        error('curvefit:fit:startPointSize', ...
              'StartPoint has %u elements but the model has %u coefficients.', ...
              numel(startPoint),nP);
    end
    if isempty(lb), lb = -Inf(nP,1); end
    if isempty(ub), ub =  Inf(nP,1); end
    if any(isfinite(lb)) || any(isfinite(ub))
        [p,r,iters,exitflag] = trustRegionReflective(m.fun,m.jac,x,y, ...
            startPoint(:),lb(:),ub(:));
    else
        starts = {startPoint(:)};
        if ~userStart && isfield(m,'start2')
            s2 = m.start2(x,y);
            starts{2} = s2(:);
        end
        [p,r,iters,exitflag] = bestFit(m,x,y,starts);
        if isfield(m,'canonical')
            p = m.canonical(p);
        end
    end
end

if any(~isfinite(r))
    error('curvefit:fit:nanComputed','NaN computed by model function.');
end

%-------------------------------------------------------------------------------
%% Package outputs
%-------------------------------------------------------------------------------
cfun = struct();
pOut = p;
if scaled && ~normalize, pOut = unscalePoly(p,mu,sg); end   % fit() reports raw-x coefficients
for k = 1:nP
    cfun.(m.coeffs{k}) = pOut(k);
end
% Retained so BF_fitEval can reconstruct the model; leading underscore-free names
% would collide with coefficient names like 'a' or 'p1'.
cfun.BF_modelName = modelName;
cfun.BF_coeffs    = p(:)';       % (for a scaled polyN: in the scaled variable)
if scaled, cfun.BF_scale = [mu sg]; end

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
function checkUnsupported(opt)
% fitoptions properties BF_fit does not implement are an error when set (HCTSA2
% 2.2.1), not silently dropped. Robust, Weights and Exclude at their defaults
% ('off', [], []) are not "set". (Normalize is implemented: HCTSA2 3.0.0.)
f = fieldnames(opt);
for i = 1:numel(f)
    v = opt.(f{i});
    switch f{i}
    case 'Robust'
        bad = ~(ischar(v) && strcmpi(v,'off'));
    case 'Normalize'                     % implemented (HCTSA2 3.0.0); checked by onOff
        bad = false;
    case {'Weights','Exclude'}
        bad = ~isempty(v);
    case {'Algorithm','TolFun','TolX','MaxFunEvals','DiffMinChange','DiffMaxChange'}
        bad = true;
    otherwise
        bad = false;
    end
    if bad
        error('BF_fit:unsupportedOption','The fit option ''%s'' is not implemented.',f{i});
    end
end
end

%-------------------------------------------------------------------------------
function tf = onOff(v)
% fit()'s 'on'/'off' option values
if isstring(v), v = char(v); end
if ischar(v) && any(strcmpi(v,{'on','off'}))
    tf = strcmpi(v,'on');
else
    error('BF_fit:badOptions','Normalize must be ''on'' or ''off''.');
end
end

%-------------------------------------------------------------------------------
function pr = unscalePoly(ps,mu,sg)
% Coefficients (highest power first) of a polynomial in z = (x-mu)/sg, as the
% coefficients of the same polynomial in x
n = numel(ps) - 1;
pr = zeros(n+1,1);
zk = 1;                                  % z^k as a polynomial in x
for k = 0:n
    pr(end-k:end) = pr(end-k:end) + ps(end-k)*zk(:);
    zk = conv(zk,[1/sg, -mu/sg]);
end
end

%-------------------------------------------------------------------------------
function [p,r,iters,exitflag] = bestFit(m,x,y,starts)
% Unbounded fits: every solver from every start; the first candidate is kept
% unless another's sum of squares is lower by more than 1e-6 (relative), so that
% equally good fits resolve to fit()'s method from fit()'s start.
nP = numel(starts{1});
free = Inf(nP,1);
sep = isfield(m,'separable');
p = starts{1}; r = NaN(size(y)); iters = 0; exitflag = -1;
best = Inf;
for i = 1:numel(starts)
    s0 = starts{i};
    for solver = 1:4
        switch solver
        case 1      % fit()'s method, stopped by fit()'s tolerances as it applies them
            if sep
                [pc,rc,ic,ec] = varPro(m.separable,x,y,s0,true);
            else
                [pc,rc,ic,ec] = trustRegionReflective(m.fun,m.jac,x,y,s0,-free,free,true);
            end
            p1 = pc; i1 = ic;
        case 2      % ... continued to convergence from where it stopped
            if ~all(isfinite(p1)), continue; end
            if sep
                [pc,rc,ic,ec] = varPro(m.separable,x,y,p1,false);
            else
                [pc,rc,ic,ec] = trustRegionReflective(m.fun,m.jac,x,y,p1,-free,free);
            end
            ic = i1 + ic;
        case 3      % trust region on all coefficients (separable models)
            if ~sep, continue; end
            [pc,rc,ic,ec] = trustRegionReflective(m.fun,m.jac,x,y,s0,-free,free);
        case 4      % Levenberg-Marquardt (HCTSA2 2.0)
            [pc,rc,ic,ec] = levenbergMarquardt(m.fun,m.jac,x,y,s0);
        end
        sse = sum(rc.^2);
        if isfinite(sse) && (~isfinite(best) || sse < best*(1 - 1e-6))
            p = pc; r = rc; iters = ic; exitflag = ec;
            best = sse;
        end
    end
end
end

%-------------------------------------------------------------------------------
function [p,r,iters,exitflag] = varPro(sep,x,y,p0,fitStop)
% Separable least squares (variable projection, Golub & Pereyra): the coefficients
% the model is linear in are eliminated by solving for them at every step, and
% the trust-region method minimises over the others alone, with forward-difference
% derivatives -- as fit() solves its separable library models.
nl = setdiff(1:numel(p0),sep.lin);
fun = @(t,xx) projectedFit(sep.basis,t,xx,y);
jac = @(t,xx) fdJacobian(fun,t,xx);
free = Inf(numel(nl),1);
[t,~,iters,exitflag] = trustRegionReflective(fun,jac,x,y,p0(nl),-free,free,fitStop);
A = sep.basis(t,x);
w = warning('off','all');
a = A\y;
warning(w);
p = zeros(numel(p0),1);
p(sep.lin) = a;
p(nl) = t;
r = y - A*a;
end

function f = projectedFit(basis,t,x,y)
A = basis(t,x);
w = warning('off','all');
f = A*(A\y);
warning(w);
end

function J = fdJacobian(fun,t,x)
f0 = fun(t,x);
J = zeros(numel(f0),numel(t));
for k = 1:numel(t)
    h = sqrt(eps)*max(abs(t(k)),1);
    tk = t;
    tk(k) = tk(k) + h;
    J(:,k) = (fun(tk,x) - f0)/h;
end
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
exitflag = 0;    % stays 0 if maxIter is reached without converging

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

end

%-------------------------------------------------------------------------------
function [x,r,iter,exitflag] = trustRegionReflective(fun,jac,t,y,x,lb,ub,fitStop)
% Bound-constrained nonlinear least squares, for fits with Lower/Upper bounds.
%
% The trust-region-reflective method of Branch, Coleman & Li (SIAM J. Sci.
% Comput. 21(1), 1999) -- the family MATLAB's fit() uses when bounds are given
% (its 'Trust-Region' algorithm) -- structured as in SciPy's least_squares
% 'trf' (BSD-3-Clause): Coleman-Li scaling by the distance to the bound the
% gradient points at, an exact trust-region step in the scaled variables, and
% a choice between that step (kept strictly interior), its reflection off the
% bound it hits, and the scaled steepest-descent step.
%
% Why this and not a projected Levenberg-Marquardt: SP_Specparam's peak and knee
% models have several local minima and flat valleys, and where a solver ends up
% depends on its path. A projected LM reached a different local minimum from
% fit() on 3 of the 5 peak fits logged in the licensed run of 2026-09-23; this
% method reaches fit()'s on all 5.
%
% Stopping rules are fit()'s defaults, deliberately (TolFun = TolX = 1e-6,
% MaxIter 400, MaxFunEvals 600): in a flat valley the stopping point IS the
% answer, and stopping where fit() stops keeps the features comparable. So the
% coefficients agree with fit()'s to roughly the tolerance, not to machine
% precision (the two methods differ in detail, e.g. fit() differentiates
% numerically).
%
% FITSTOP (unbounded fits, bestFit's first candidate) applies the two tests in the
% form MATLAB documents for this solver, |change in SSE| < TolFun*(1 + SSE) and
% ||step|| < TolX*(1 + ||p||), which stop earlier than the relative forms below
% when the SSE is small or the coefficients large. Without it (and always on the
% bounded path) the tests are relative, as in SciPy.

if nargin < 8, fitStop = false; end
ftol = 1e-6; xtol = 1e-6; gtol = 1e-8;
maxIter = 400; maxFev = 600;
nP = numel(x);
m = numel(y);

x = makeStrictlyFeasible(x,lb,ub,1e-10);
f = fun(x,t) - y;                           % residuals, SciPy's sign convention
nfev = 1;
cost = 0.5*(f'*f);
J = jac(x,t);
g = J'*f;
[v,~] = clScaling(x,g,lb,ub);
Delta = norm(x./sqrt(v));
if Delta == 0, Delta = 1; end
alpha = 0;
exitflag = 0;
iter = 0;

while true
    [v,dv] = clScaling(x,g,lb,ub);
    gNorm = max(abs(g.*v));
    if gNorm < gtol
        exitflag = 1;
    end
    if exitflag ~= 0 || nfev >= maxFev || iter >= maxIter || ~all(isfinite(J(:)))
        break
    end
    d = sqrt(v);
    diagH = g.*dv;
    gH = d.*g;
    Jh = J.*d';
    [U,S,V] = svd([Jh; diag(sqrt(diagH))],'econ');
    s = diag(S);
    uf = U'*[f; zeros(nP,1)];
    theta = max(0.995, 1 - gNorm);

    actualReduction = -1;
    while actualReduction <= 0 && nfev < maxFev
        [pH,alpha] = lsqTrustRegion(nP,m,uf,s,V,Delta,alpha);
        p = d.*pH;
        [step,stepH,predictedReduction] = selectStep(x,Jh,diagH,gH,p,pH,d,Delta,lb,ub,theta);
        xNew = makeStrictlyFeasible(x + step,lb,ub,0);
        fNew = fun(xNew,t) - y;
        nfev = nfev + 1;
        stepHNorm = norm(stepH);
        if ~all(isfinite(fNew))
            Delta = 0.25*stepHNorm;
            continue
        end
        costNew = 0.5*(fNew'*fNew);
        actualReduction = cost - costNew;
        [DeltaNew,ratio] = updateRadius(Delta,actualReduction,predictedReduction, ...
            stepHNorm,stepHNorm > 0.95*Delta);
        if fitStop                           % cost is SSE/2
            ftolOK = actualReduction > 0 && 2*actualReduction < ftol*(1 + 2*cost);
            xtolOK = norm(step) < xtol*(1 + norm(x));
        else
            ftolOK = actualReduction < ftol*cost && ratio > 0.25;
            xtolOK = norm(step) < xtol*(xtol + norm(x));
        end
        if ftolOK || xtolOK
            exitflag = 2 + ftolOK;           % 3: change in residuals below TolFun
            break
        end
        alpha = alpha*Delta/DeltaNew;
        Delta = DeltaNew;
    end
    if actualReduction > 0
        x = xNew; f = fNew; cost = costNew;
        J = jac(x,t);
        g = J'*f;
    end
    iter = iter + 1;
end
r = -f;                                      % back to y - model
end

function x = makeStrictlyFeasible(x,lb,ub,rstep)
if rstep == 0
    x(x <= lb) = lb(x <= lb) + eps(max(1,abs(lb(x <= lb))));   % nudge off the bound
    x(x >= ub) = ub(x >= ub) - eps(max(1,abs(ub(x >= ub))));
else
    atLb = x <= lb; atUb = x >= ub;
    x(atLb) = lb(atLb) + rstep*max(1,abs(lb(atLb)));
    x(atUb) = ub(atUb) - rstep*max(1,abs(ub(atUb)));
end
tight = lb > x | x > ub;                    % bounds closer than the nudge
x(tight) = 0.5*(lb(tight) + ub(tight));
end

function [v,dv] = clScaling(x,g,lb,ub)
% Coleman-Li scaling: distance to the bound the (descent) direction points at.
v = ones(size(x)); dv = zeros(size(x));
k = g < 0 & isfinite(ub);  v(k) = ub(k) - x(k);  dv(k) = -1;
k = g > 0 & isfinite(lb);  v(k) = x(k) - lb(k);  dv(k) = 1;
end

function [p,alpha] = lsqTrustRegion(n,m,uf,s,V,Delta,alpha)
% Exact solution of min ||J p + f|| subject to ||p|| <= Delta, from the SVD of
% the (augmented) Jacobian, by the More-Hebden iteration on alpha.
suf = s.*uf;
fullRank = m >= n && s(end) > eps*m*s(1);
if fullRank
    p = -V*(uf./s);
    if norm(p) <= Delta
        alpha = 0;
        return
    end
end
alphaUpper = norm(suf)/Delta;
if fullRank
    [phi,dphi] = phiDeriv(0,suf,s,Delta);
    alphaLower = -phi/dphi;
else
    alphaLower = 0;
end
if ~fullRank && alpha == 0
    alpha = max(0.001*alphaUpper, sqrt(alphaLower*alphaUpper));
end
for it = 1:10
    if alpha < alphaLower || alpha > alphaUpper
        alpha = max(0.001*alphaUpper, sqrt(alphaLower*alphaUpper));
    end
    [phi,dphi] = phiDeriv(alpha,suf,s,Delta);
    if phi < 0, alphaUpper = alpha; end
    ratio = phi/dphi;
    alphaLower = max(alphaLower, alpha - ratio);
    alpha = alpha - (phi + Delta)*ratio/Delta;
    if abs(phi) < 0.01*Delta, break; end
end
p = -V*(suf./(s.^2 + alpha));
p = p*Delta/norm(p);
end

function [phi,dphi] = phiDeriv(alpha,suf,s,Delta)
den = s.^2 + alpha;
pNorm = norm(suf./den);
phi = pNorm - Delta;
dphi = -sum(suf.^2./den.^3)/pNorm;
end

function [step,stepH,predRed] = selectStep(x,Jh,diagH,gH,p,pH,d,Delta,lb,ub,theta)
if all(x + p >= lb & x + p <= ub)
    predRed = -quadValue(Jh,gH,pH,diagH);
    step = p; stepH = pH;
    return
end
[pStride,hits] = stepToBound(x,p,lb,ub);
rH = pH; rH(hits ~= 0) = -rH(hits ~= 0);    % reflect off the bound(s) hit
r = d.*rH;
p = p*pStride; pH = pH*pStride;
xOnBound = x + p;
[~,toTR] = intersectTR(pH,rH,Delta);
toBound = stepToBound(xOnBound,r,lb,ub);
rStride = min(toBound,toTR);
if rStride > 0
    rLo = (1 - theta)*pStride/rStride;
    if rStride == toBound, rHi = theta*toBound; else, rHi = toTR; end
else
    rLo = 0; rHi = -1;
end
if rLo <= rHi
    [a,b,c] = quad1d(Jh,gH,rH,diagH,pH);
    [rStride,rValue] = minQuad1d(a,b,rLo,rHi,c);
    rH = pH + rH*rStride;
    r = d.*rH;
else
    rValue = Inf;
end
p = p*theta; pH = pH*theta;                  % keep the truncated step interior
pValue = quadValue(Jh,gH,pH,diagH);
agH = -gH; ag = d.*agH;
toTR = Delta/norm(agH);
toBound = stepToBound(x,ag,lb,ub);
if toBound < toTR, agStride = theta*toBound; else, agStride = toTR; end
[a,b] = quad1d(Jh,gH,agH,diagH,[]);
[agStride,agValue] = minQuad1d(a,b,0,agStride,0);
agH = agH*agStride; ag = ag*agStride;
if pValue < rValue && pValue < agValue
    step = p; stepH = pH; predRed = -pValue;
elseif rValue < pValue && rValue < agValue
    step = r; stepH = rH; predRed = -rValue;
else
    step = ag; stepH = agH; predRed = -agValue;
end
end

function [minStep,hits] = stepToBound(x,s,lb,ub)
steps = Inf(size(x));
nz = s ~= 0;
steps(nz) = max((lb(nz) - x(nz))./s(nz), (ub(nz) - x(nz))./s(nz));
minStep = min(steps);
hits = (steps == minStep).*sign(s);
end

function [t1,t2] = intersectTR(x,s,Delta)
a = s'*s; b = x'*s; c = x'*x - Delta^2;
dd = sqrt(b^2 - a*c);
q = -(b + sign(b + (b == 0))*dd);
t1 = q/a; t2 = c/q;
if t1 > t2, [t1,t2] = deal(t2,t1); end
end

function q = quadValue(J,g,s,diagH)
Js = J*s;
q = 0.5*(Js'*Js + s'*(diagH.*s)) + s'*g;
end

function [a,b,c] = quad1d(J,g,s,diagH,s0)
v = J*s;
a = 0.5*(v'*v + s'*(diagH.*s));
b = g'*s;
c = 0;
if ~isempty(s0)
    u = J*s0;
    b = b + u'*v + s0'*(diagH.*s);
    c = 0.5*(u'*u) + g'*s0 + 0.5*s0'*(diagH.*s0);
end
end

function [tBest,yBest] = minQuad1d(a,b,lo,hi,c)
tt = [lo, hi];
if a ~= 0
    ext = -0.5*b/a;
    if lo < ext && ext < hi, tt(end+1) = ext; end
end
yy = tt.*(a*tt + b) + c;
[yBest,k] = min(yy);
tBest = tt(k);
end

function [Delta,ratio] = updateRadius(Delta,actual,predicted,stepNorm,boundHit)
if predicted > 0
    ratio = actual/predicted;
elseif predicted == actual && actual == 0
    ratio = 1;
else
    ratio = 0;
end
if ratio < 0.25
    Delta = 0.25*stepNorm;
elseif ratio > 0.75 && boundHit
    Delta = 2*Delta;
end
end
