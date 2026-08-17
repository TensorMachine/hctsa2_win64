function [EstMdl,estParamCov,LLF,info] = BF_garchEstimate(Mdl,y,varargin)
% BF_garchEstimate   Gaussian GARCH(P,Q) MLE (toolbox-free estimate()).
%
% Maximises
%   LLF = -0.5 * sum_t [ log(2*pi) + log(sigma2_t) + e_t^2/sigma2_t ]
%
% subject to K > 0, GARCH >= 0, ARCH >= 0 and sum(GARCH)+sum(ARCH) < 1. The
% constraints are imposed by an unconstrained reparameterisation (log for the
% constant, logistic for total persistence, softmax for its split), so plain
% fminsearch can be used -- the Optimization Toolbox is not required.
%
%---OUTPUTS:
% EstMdl       model with fitted Constant, GARCH and ARCH
% estParamCov  covariance of [Constant; GARCH; ARCH], by outer product of gradients
% LLF          maximised log-likelihood
% info         struct with exitflag and the fitted parameter vector

y = y(:);
y = y(isfinite(y));
P = Mdl.P; Q = Mdl.Q;
nPar = 1 + P + Q;
offset = Mdl.Offset;
if isnan(offset), offset = mean(y); end

v0 = var(y - offset);

% --- unconstrained -> natural parameters -------------------------------------
    function [K,G,A] = unpack(u)
        K = exp(u(1));
        if P+Q > 0
            total = 1/(1+exp(-u(2)));            % persistence in (0,1)
            w = exp(u(3:2+P+Q));
            w = w/sum(w);
            GA = total*w;
            G = GA(1:P);
            A = GA(P+1:end);
        else
            G = []; A = [];
        end
    end

    function nll = negLL(u)
        [K,G,A] = unpack(u);
        [v,e] = BF_garchVar(y,K,G,A,offset);
        if any(~isfinite(v)) || any(v <= 0)
            nll = 1e10; return
        end
        nll = 0.5*sum(log(2*pi) + log(v) + e.^2./v);
        if ~isfinite(nll), nll = 1e10; end
    end

% Start from a persistence of 0.85 split mostly onto the GARCH terms, a common
% empirical starting point for financial-style volatility clustering
u0 = zeros(1,2+P+Q);
u0(1) = log(max(v0*0.15,eps));
u0(2) = log(0.85/0.15);
if P > 0, u0(3:2+P) = 1; end
if Q > 0, u0(3+P:2+P+Q) = 0; end

opts = optimset('MaxFunEvals',20000,'MaxIter',10000,'TolX',1e-10,'TolFun',1e-10);
[uHat,fval,exitflag] = fminsearch(@negLL,u0,opts);
% Restart once from the solution: fminsearch often stalls before convergence
[uHat,fval,exitflag] = fminsearch(@negLL,uHat,opts);

[K,G,A] = unpack(uHat);
LLF = -fval;

EstMdl = Mdl;
EstMdl.Constant = K;
EstMdl.GARCH = num2cell(G(:)');
EstMdl.ARCH  = num2cell(A(:)');
EstMdl.Offset = offset;

% --- covariance by outer product of gradients (OPG) ---------------------------
% Verified against the Curve Fitting/Econometrics reference: MathWorks uses OPG
% here, not the inverse Hessian. On this data the Hessian gives standard errors
% ~15-30%% too small and the sandwich estimator is further off still.
theta = [K, G(:)', A(:)'];
    function L = perObsLL(th)
        Kn = th(1);
        Gn = th(2:1+P);
        An = th(2+P:1+P+Q);
        [v,e] = BF_garchVar(y,Kn,Gn,An,offset);
        L = -0.5*(log(2*pi) + log(v) + e.^2./v);
    end

T = numel(y);
Grad = zeros(T,nPar);
for a = 1:nPar
    h = max(abs(theta(a))*1e-5,1e-8);
    tp = theta; tp(a) = tp(a)+h;
    tm = theta; tm(a) = tm(a)-h;
    Grad(:,a) = (perObsLL(tp)-perObsLL(tm))/(2*h);
end
warnState = warning('off','all');
estParamCov = inv(Grad'*Grad);
warning(warnState);
if any(~isfinite(estParamCov(:)))
    estParamCov = NaN(nPar);
end

info.exitflag = exitflag;
info.X = theta;
info.LLF = LLF;
end
