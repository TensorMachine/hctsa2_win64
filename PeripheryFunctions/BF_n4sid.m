function m = BF_n4sid(y,n,r,s1,W)
% BF_n4sid   Subspace state-space identification (toolbox-free n4sid()).
%
% Estimates the innovations-form model
%   x(t+1) = A x(t) + K e(t)
%   y(t)   = C x(t) + e(t)
% from output-only data, which is the only case hctsa uses.
%
%---ALGORITHM
% This follows the formulation in Ljung, "System Identification: Theory for the
% User" (1999), Section 10.6, which is what the System Identification Toolbox
% documentation cites for n4sid. Three details matter and are easy to get wrong:
%
%  1. Everything is built from the **LQ factorisation of [Phi; Y] without
%     normalising by N**. That is where the scale of C comes from: the reported
%     c_1 sits near sqrt(N) (about 63 for these 4097-sample series), which is
%     why it looks nearly constant across series and model orders.
%  2. **A and C come from the shift-invariance of the extended observability
%     matrix Or**, not from a least-squares fit of a state sequence.
%  3. The Kalman gain comes from residual covariances formed out of the L
%     factor (Ljung's find_PK construction), not from state-sequence residuals.
%
% Using a balanced SVD split, a state-sequence least squares for A and C, or
% covariance-based stochastic realization instead all give errors of 3-30%
% against reference output from a licensed MATLAB; this formulation gives ~1e-4.
%
%---INPUTS:
% y   time series
% n   model order
% r   forward prediction horizon        (default 2)
% s1  number of past outputs            (default 3)
% W   'CVA' (default) or 'MOESP'
%
% The defaults r = 2, s1 = 3 were determined against reference output: they are
% optimal for all 12 reference series tested, so the horizon appears fixed here
% rather than data-dependent.

% Order selection: n = 'best' picks the order from the singular-value spectrum,
% using the criterion in the reference implementation of Ljung's algorithm
% (count the singular values above the geometric mean of the largest and
% smallest).
autoOrder = false;
if ischar(n) && strcmpi(n,'best')
    autoOrder = true;
    n = 1;                   % provisional, refined below
end

% Default horizons, determined against reference output. The past horizon
% follows s1 = 4n-1 exactly (3, 7, 11 for orders 1, 2, 3); the forward horizon
% was found by search (2, 3, 5) and is extrapolated as 2n-1 beyond order 3,
% which is consistent with the observed values at n = 2 and 3.
% Horizon rule, read directly from MATLAB's own Report.N4Horizon on a licensed
% install (R2025b, System Identification Toolbox 25.2):
%
%   order:      1        2        3        4        10
%   N4Horizon:  [2 3 3]  [3 7 7]  [5 11 11] [6 15 15] [15 39 39]
%
% giving  r = ceil(3n/2)  and  sy = 4n-1  exactly at every observed order.
% `sy` is the AIC-optimal ARX order capped at 4n (BF_n4sidAIC is a direct port
% of localAIC in n4sid_time.m). It is usually 4n-1 but not always: on Bonn
% series 1 it is 34 at order 9 and 46 at order 12. Matches MATLAB at 12 of 12
% orders tested.
if nargin < 4 || isempty(s1), s1 = BF_n4sidAIC(y(isfinite(y)),4*n); end
if nargin < 3 || isempty(r),  r  = ceil(3*n/2);  end
if nargin < 5 || isempty(W),  W  = 'CVA'; end

y = y(:);
y = y(isfinite(y));
T = numel(y);
p = 1;                       % single output
t0 = s1 + 1;
s  = s1*p;
N  = T - r + 1 - t0;
if N < 50 || r < n+1
    error('BF_n4sid:tooShort','Series too short for order %u with r=%u, s1=%u.',n,r,s1);
end

% Block-Hankel of futures, and the past regressors (reverse time order)
Y = zeros(r,N);
for ri = 1:r
    Y(ri,:) = y(t0+ri-1 : t0+ri+N-2)';
end
Phi = zeros(s,N);
for k = 1:s1
    Phi(k,:) = y(t0-k : t0+N-1-k)';
end

if autoOrder
    % Order selection, following MATLAB's n4sid_time.m:
    %
    %  - the horizon is computed from the TOP of the order range (10), not from
    %    the order eventually chosen, and is kept for the final estimate. The
    %    probe confirms this: a synthetic AR(2) selects order 2 but still
    %    reports N4Horizon [15 10 10], i.e. r = ceil(1.5*10) = 15.
    %  - **the order is chosen from the UNWEIGHTED singular values** svd(L32).
    %    For CVA the source explicitly recomputes them:
    %        if strcmp(n4w,'CVA'), [~,Sn1,~] = svd(R(ind3,ind2)); end
    %    Applying the criterion to the CVA canonical correlations instead gave
    %    9, 8, 9 where MATLAB gives 10, 10, 10.
    %  - the count is capped at 10, the top of the default range.
    nMax = 10;
    rA = ceil(1.5*nMax);
    sA = BF_n4sidAIC(y,4*nMax);
    NA = T - rA + 1 - (sA+1);
    if NA > 50
        t0a = sA + 1;
        Ya = zeros(rA,NA); Pa = zeros(sA,NA);
        for ri = 1:rA, Ya(ri,:) = y(t0a+ri-1 : t0a+ri+NA-2)'; end
        for k = 1:sA,  Pa(k,:)  = y(t0a-k : t0a+NA-1-k)'; end
        [~,Ra] = qr([Pa;Ya]',0);
        La = Ra';
        sva = svd(La(sA+1:end,1:sA));
        n = min(nMax,max(1,sum(sva > sqrt(sva(1)*sva(end)))));
    else
        n = 1;
    end
    % the horizon stays at the top-of-range value
    r = rA; s1 = sA;
    t0 = s1 + 1; s = s1*p; N = T - r + 1 - t0;
    Y = zeros(r,N);
    for ri = 1:r, Y(ri,:) = y(t0+ri-1 : t0+ri+N-2)'; end
    Phi = zeros(s,N);
    for k = 1:s1, Phi(k,:) = y(t0-k : t0+N-1-k)'; end
end

[Qr,Rr] = qr([Phi;Y]',0);
L = Rr';
Q = Qr';
iPhi = 1:s;
iY   = s+1 : s+r*p;
L32 = L(iY,iPhi);
L2  = L(iY, 1:s+p);

switch upper(W)
case 'MOESP'
    G = L32*Q(iPhi,:);
    [U1,S1,~] = svd(G,'econ');
    sv = diag(S1);
    Or = U1(:,1:n)*diag(sqrt(sv(1:n)));
otherwise   % CVA
    W1 = L(iY,[iPhi iY]);
    [ul,sl,~] = svd(W1,'econ');
    sl = diag(sl);
    sl = diag(sl(1:r*p));
    [Or2,~,~] = svd(pinv(sl)*ul'*L32,'econ');
    Or = ul*sl*Or2;
    Or = Or(:,1:n);
end

C = Or(1:p,1:n);
A = Or(1:p*(r-1),1:n) \ Or(p+1:p*r,1:n);

% At high order the shift-invariance solve is ill-conditioned (the CVA whitening
% matrix reaches a condition number of ~2.5e5 at order 10) and can return a
% single spurious unstable mode. Comparing against MATLAB at order 10: nine of
% ten eigenvalues agree to ~4e-4, while a spurious 1.3785 replaces their 0.7427.
% Reflecting unstable eigenvalues inside the unit circle -- the standard
% minimum-phase operation, which preserves the magnitude response -- maps
% 1.3785 to 1/1.3785 = 0.7254 and keeps the predictor recursion bounded.
if any(abs(eig(A)) >= 1)
    [Vv,Dd] = eig(A);
    dd = diag(Dd);
    bad = abs(dd) >= 1;
    dd(bad) = conj(1./dd(bad));
    A = real(Vv*diag(dd)/Vv);
end

% Kalman gain: residual covariances from the L factor, then the Riccati
X1 = L2(p+1:r*p, 1:p*s1+p);
X2 = [L2(1:r*p,1:p*s1), zeros(r*p,p)];
vl = [Or(1:(r-1)*p,1:n)\X1 ; L2(1:p, 1:p*s1+p)];
hl = Or(:,1:n)\X2;
K0 = vl*pinv(hl);
Wm = (vl-K0*hl)*(vl-K0*hl)';
Wm = Wm/(numel(vl)-numel(K0));
Qk = Wm(1:n,1:n);
Sk = Wm(1:n,n+1:n+p);
Rk = Wm(n+1:n+p,n+1:n+p);

% The Riccati fixed-point diverges if the extracted A is unstable, which can
% happen at high order where the shift-invariance solve is ill-conditioned.
% Guard it: detect divergence and fall back to a stabilised A (eigenvalues
% projected just inside the unit circle) for the gain computation only.
[K,ok] = riccatiGain(A,C,Qk,Sk,Rk,n);
if ~ok
    [V,D] = eig(A);
    d = diag(D);
    bad = abs(d) >= 1;
    if any(bad)
        d(bad) = 0.99*d(bad)./abs(d(bad));
        As = real(V*diag(d)/V);
    else
        As = A;
    end
    [K,ok2] = riccatiGain(As,C,Qk,Sk,Rk,n);
    if ~ok2
        K = Sk*pinv(Rk);          % last resort: one-shot gain
    end
end

% Initial state and the prediction-error loss, with the initial state estimated
% (n4sid's InitialState default is 'estimate')
[V,x0,e] = BF_ssLoss(A,C,K,y);

% Parameter count used by FPE and AIC: 2n identifiable model parameters (A, K, C
% modulo the similarity transform) plus n estimated initial states, since
% n4sid's InitialState default is 'estimate'.
%
% Order 1 cannot distinguish this from n^2+2n (both give 3); order 2 can, and
% the reference settles it: d = 3n gives 1.75e-04 against 1.15e-03 for n^2+2n,
% and it removes a uniform ~1% offset across every fpe field of
% MF_FitSubsegments_ss_2.
d = 3*n;
m.A = A;  m.C = C;  m.K = K;
m.k = K;
m.X0 = x0;
m.NoiseVariance = V;
m.np = d;
m.N = T;
m.Ts = 1;
m.order = n;
m.EstimationInfo.LossFcn = V;
m.EstimationInfo.FPE = V*(T+d)/(T-d);
m.EstimationInfo.N4Horizon = [r s1 0];
m.ParameterVector = [A(:); K(:); C(:)];
m.residuals = e;
end

%-------------------------------------------------------------------------------
function [K,ok] = riccatiGain(A,C,Qk,Sk,Rk,n)
% Iterate the algebraic Riccati equation, reporting failure rather than
% returning Inf/NaN. The fixed point diverges when A is unstable, which the
% shift-invariance solve can produce at high order.
P = zeros(n); ok = false;
for it = 1:5000
    Pn = A*P*A' + Qk - (A*P*C'+Sk)*pinv(C*P*C'+Rk)*(A*P*C'+Sk)';
    if ~all(isfinite(Pn(:))) || norm(Pn,'fro') > 1e12
        K = zeros(n,1); return
    end
    if max(abs(Pn(:)-P(:))) < 1e-15, P = Pn; ok = true; break; end
    P = Pn;
end
K = (A*P*C'+Sk)*pinv(C*P*C'+Rk);
ok = ok && all(isfinite(K(:)));
end
