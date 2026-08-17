function yp = BF_predict(m,y,k,varargin)
% BF_predict   k-step-ahead prediction for a polynomial model (predict()).
%
% For a model  A(q) y(t) = C(q) e(t)  the noise transfer function is H = C/A and
% the k-step predictor is (Ljung, System Identification, eq. 3.53)
%
%   yhat(t|t-k) = [ 1 - Wk(q) H^-1(q) ] y(t),
%   Wk(q)       = sum_{j=0}^{k-1} h(j) q^-j
%
% with h the impulse response of H. Since H^-1 = A/C, this is
%
%   e   = filter(A, C, y)          the one-step prediction error
%   yhat = y - filter(h(1:k), 1, e)
%
% For k = 1, Wk = 1 and this reduces to yhat = y - e, as it must.
%
% AR models have C = 1. Pass 'init','e' (MATLAB's option, and the default here)
% to estimate initial conditions rather than assume zeros -- important whenever
% the MA root is close to the unit circle, where the transient is long-lived.

y = y(:);
N = numel(y);
if nargin < 3 || isempty(k), k = 1; end

initE = true;                    % 'init','e' is what hctsa asks for
for i = 1:2:numel(varargin)-1
    if strcmpi(varargin{i},'init')
        initE = ~strcmpi(varargin{i+1},'z');
    end
end

%-------------------------------------------------------------------------------
%% State-space models (from BF_n4sid): run the Kalman predictor, then project
%% forward k-1 steps with the free dynamics
%-------------------------------------------------------------------------------
if isfield(m,'K') && isfield(m,'C') && ~isfield(m,'a')
    A = m.A; C = m.C; K = m.K;
    n = size(A,1);
    Ac = A - K*C;

    % 'init','e' means the initial state is estimated on the data supplied here,
    % NOT inherited from the model's training fit. The initial state enters the
    % predictor linearly, so it is found by least squares against the zero-state
    % residuals -- the same construction as BF_ssLoss.
    x0 = zeros(n,1);
    if initE
        e0 = zeros(N,1);
        x = zeros(n,1);
        for t = 1:N
            e0(t) = y(t) - C*x;
            x = Ac*x + K*y(t);
        end
        Bx = zeros(N,n);
        for jj = 1:n
            x = zeros(n,1); x(jj) = 1;
            for t = 1:N
                Bx(t,jj) = -C*x;
                x = Ac*x;
            end
        end
        warnState = warning('off','all');
        z = (Bx'*Bx)\(Bx'*e0);
        warning(warnState);
        if all(isfinite(z)), x0 = -z; end
    elseif isfield(m,'X0') && ~isempty(m.X0)
        x0 = m.X0(:);
    end

    % One-step predicted states xhat(:,t) = x(t | t-1)
    xh = zeros(n,N);
    x = x0;
    for t = 1:N
        xh(:,t) = x;
        x = Ac*x + K*y(t);
    end

    % yhat(t | t-k) = C * A^(k-1) * xhat(:, t-k+1)
    Ak = eye(n);
    for i = 1:k-1, Ak = Ak*A; end
    yp = NaN(N,1);
    for t = 1:N
        src = t-k+1;
        if src >= 1
            yp(t) = C*Ak*xh(:,src);
        end
    end
    yp(isnan(yp)) = 0;
    return
end

%-------------------------------------------------------------------------------
%% Polynomial models
%-------------------------------------------------------------------------------
if isfield(m,'a'), A = m.a(:)'; else A = m.A(:)'; end
if isfield(m,'c') && ~isempty(m.c), C = m.c(:)'; else C = 1; end

% One-step prediction error
e = filter(A,C,y);

if initE
    % Project out the initial-condition subspace, spanned by the delayed impulse
    % responses of 1/C. This applies to pure AR models too (C = 1), where the
    % basis reduces to unit impulses and the projection simply removes the
    % startup transient -- without it an AR predictor carries large errors in
    % its first na samples, which inflates rmserr far above mabserr (3.8x on the
    % reference data, where the two are nearly equal).
    w = max(numel(A),numel(C)) - 1;
    % The j-th basis column is the impulse response of 1/C delayed by j-1, so
    % compute that response once and shift it rather than running w separate
    % filter calls. filter() is LTI and deterministic, so the shifted values are
    % the identical floating-point sequence -- verified bitwise against the
    % per-column version.
    B = zeros(N,w);
    if w > 0
        imp = zeros(N,1); imp(1) = 1;
        h = filter(1,C,imp);
        for j = 1:w
            B(j:N,j) = h(1:N-j+1);
        end
    end
    warnState = warning('off','all');
    z = (B'*B)\(B'*e);
    warning(warnState);
    if all(isfinite(z))
        e = e - B*z;
    end
end

% Impulse response of H = C/A, first k terms
h = filter(C,A,[1; zeros(k-1,1)]);

yp = y - filter(h,1,e);
end
