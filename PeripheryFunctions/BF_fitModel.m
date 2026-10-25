function m = BF_fitModel(name)
% BF_fitModel   Model registry for BF_fit (toolbox-free Curve Fitting replacement).
%
% Returns a struct describing model NAME:
%   coeffs  cell array of coefficient names, in parameter order
%   linear  true if the model is linear in its parameters (closed-form solve)
%   basis   @(x) design matrix                          [linear models only]
%   fun     @(p,x) model values                         [nonlinear models only]
%   jac     @(p,x) Jacobian, numel(x)-by-numel(p)       [nonlinear models only]
%   start   @(x,y) default start point                  [nonlinear models only]
% and, for some library models only:
%   start2     @(x,y) a second start point: HCTSA2 2.0's rule where fit()'s is now
%              the default (exp1, power1, gaussN, sinN)
%   separable  struct: lin, the indices of the coefficients the model is linear in,
%              and basis, @(t,x) the matrix they multiply, t being the other
%              coefficients in order. fit() solves exactly these models separably
%              (exp1, power1, gaussN), and so does BF_fit.
%   canonical  @(p) the same fit with fit()'s sign convention: gaussN widths c >= 0
%              and sinN frequencies b >= 0, the Lower bounds fit() gives them by
%              default. The model values do not change.
%
% Covers every model form used anywhere in hctsa. See BF_fit.

if isstring(name), name = char(name); end   % "poly2" as well as 'poly2' (numel("poly2") is 1)
name = strtrim(name);
x2c = @(x) x(:); % force column

% --- polyN: linear, coefficients p1..p(N+1), highest power first ---------------
if numel(name) == 5 && strcmp(name(1:4),'poly') && any(name(5) == '123456789')
    n = str2double(name(5));
    m.coeffs = arrayfun(@(k) sprintf('p%d',k), 1:n+1, 'UniformOutput', false);
    m.linear = true;
    m.basis  = @(x) bsxfun(@power, x2c(x), n:-1:0);
    m.centre = n >= 5;   % BF_fit solves in (x - mean(x))/std(x): raw powers are too ill-conditioned
    return
end

% --- sinN: a1*sin(b1*x+c1) + ... ----------------------------------------------
if numel(name) == 4 && strcmp(name(1:3),'sin') && any(name(4) == '12345678')
    n = str2double(name(4));
    m.coeffs = reshape(arrayfun(@(k) {sprintf('a%d',k),sprintf('b%d',k),sprintf('c%d',k)}, ...
                       1:n, 'UniformOutput', false), 1, []);
    m.coeffs = [m.coeffs{:}];
    m.linear = false;
    m.fun    = @(p,x) sinSum(p,x2c(x),n);
    m.jac    = @(p,x) sinJac(p,x2c(x),n);
    m.start  = @(x,y) sinFFTStart(x2c(x),x2c(y),n);
    m.start2 = @(x,y) sinStart(x2c(x),x2c(y),n);
    m.canonical = @sinCanonical;
    return
end

% --- gaussN: a1*exp(-((x-b1)/c1)^2) + ... -------------------------------------
if numel(name) == 6 && strcmp(name(1:5),'gauss') && any(name(6) == '12345678')
    n = str2double(name(6));
    m.coeffs = reshape(arrayfun(@(k) {sprintf('a%d',k),sprintf('b%d',k),sprintf('c%d',k)}, ...
                       1:n, 'UniformOutput', false), 1, []);
    m.coeffs = [m.coeffs{:}];
    m.linear = false;
    m.fun    = @(p,x) gaussSum(p,x2c(x),n);
    m.jac    = @(p,x) gaussJac(p,x2c(x),n);
    m.start  = @(x,y) gaussPeelStart(x2c(x),x2c(y),n);
    m.start2 = @(x,y) gaussStart(x2c(x),x2c(y),n);
    m.separable = struct('lin',1:3:3*n, ...
        'basis',@(t,x) exp(-((x2c(x) - t(1:2:end)')./t(2:2:end)').^2));
    m.canonical = @(p) gaussCanonical(p,n);
    return
end

% --- named / custom model strings ---------------------------------------------
% Whitespace is stripped so 'a*x +b' and 'a*x+b' are the same model.
switch regexprep(name,'\s','')
case {'a*x+b'}
    m.coeffs = {'a','b'};  m.linear = true;
    m.basis  = @(x) [x2c(x), ones(numel(x),1)];

case 'exp1'
    % The library model: fit() solves it separably from its own start rule. The
    % same formula written out ('a*exp(b*x)') is a custom equation to fit(), and
    % hctsa always gives that one a StartPoint.
    m.coeffs = {'a','b'};  m.linear = false;
    m.fun    = @(p,x) p(1)*exp(p(2)*x2c(x));
    m.jac    = @(p,x) [exp(p(2)*x2c(x)), p(1)*x2c(x).*exp(p(2)*x2c(x))];
    m.start  = @(x,y) expSumsStart(x2c(x),x2c(y));
    m.start2 = @(x,y) expStart(x2c(x),x2c(y));
    m.separable = struct('lin',1,'basis',@(t,x) exp(t*x2c(x)));

case 'a*exp(b*x)'
    m.coeffs = {'a','b'};  m.linear = false;
    m.fun    = @(p,x) p(1)*exp(p(2)*x2c(x));
    m.jac    = @(p,x) [exp(p(2)*x2c(x)), p(1)*x2c(x).*exp(p(2)*x2c(x))];
    m.start  = @(x,y) expStart(x2c(x),x2c(y));

case 'a*exp(b*x)+c'
    m.coeffs = {'a','b','c'};  m.linear = false;
    m.fun    = @(p,x) p(1)*exp(p(2)*x2c(x)) + p(3);
    m.jac    = @(p,x) [exp(p(2)*x2c(x)), p(1)*x2c(x).*exp(p(2)*x2c(x)), ones(numel(x),1)];
    m.start  = @(x,y) [expStart(x2c(x),x2c(y)-min(y)), min(y)];

case 'a*(1-exp(b*x))'
    m.coeffs = {'a','b'};  m.linear = false;
    m.fun    = @(p,x) p(1)*(1 - exp(p(2)*x2c(x)));
    m.jac    = @(p,x) [1-exp(p(2)*x2c(x)), -p(1)*x2c(x).*exp(p(2)*x2c(x))];
    m.start  = @(x,y) [max(y), -1/max(eps,mean(x2c(x)))];

case 'exp(-b*x)'
    m.coeffs = {'b'};  m.linear = false;
    m.fun    = @(p,x) exp(-p(1)*x2c(x));
    m.jac    = @(p,x) -x2c(x).*exp(-p(1)*x2c(x));
    m.start  = @(x,y) 1/max(eps,mean(x2c(x)));

case 'a*x^2/(b+x^2)'
    m.coeffs = {'a','b'};  m.linear = false;
    m.fun    = @(p,x) p(1)*x2c(x).^2 ./ (p(2) + x2c(x).^2);
    m.jac    = @(p,x) [x2c(x).^2./(p(2)+x2c(x).^2), ...
                      -p(1)*x2c(x).^2./(p(2)+x2c(x).^2).^2];
    m.start  = @(x,y) [max(y), median(x2c(x).^2)];

case 'power1'
    % a*x^b -- requires positive x, matching Curve Fitting Toolbox behaviour
    m.coeffs = {'a','b'};  m.linear = false;
    m.fun    = @(p,x) p(1)*x2c(x).^p(2);
    m.jac    = @(p,x) [x2c(x).^p(2), p(1)*x2c(x).^p(2).*log(max(x2c(x),realmin))];
    m.start  = @(x,y) powerRatioStart(x2c(x),x2c(y));
    m.start2 = @(x,y) powerStart(x2c(x),x2c(y));
    m.separable = struct('lin',1,'basis',@(t,x) x2c(x).^t);

case 'h*exp(-(x-m)^2/(2*w^2))'
    % SP_Specparam's spectral peak: a Gaussian in log-frequency (coefficient
    % order h, m, w as SP_Specparam declares it). No default start: the caller
    % always supplies one, as it must for fit() too.
    m.coeffs = {'h','m','w'};  m.linear = false;
    m.fun    = @(p,x) p(1)*exp(-(x2c(x)-p(2)).^2/(2*p(3)^2));
    m.jac    = @(p,x) gaussPeakJac(p,x2c(x));
    m.start  = @(x,y) error('BF_fit:noStart','This model needs a StartPoint.');

case 'a-log10(k+x^c)'
    % SP_Specparam's aperiodic component with a knee (coefficients a, k, c).
    m.coeffs = {'a','k','c'};  m.linear = false;
    m.fun    = @(p,x) p(1) - log10(p(2) + x2c(x).^p(3));
    m.jac    = @(p,x) kneeJac(p,x2c(x));
    m.start  = @(x,y) error('BF_fit:noStart','This model needs a StartPoint.');

otherwise
    error('BF_fit:unknownModel','Unknown model ''%s''.',name);
end

end

%-------------------------------------------------------------------------------
% Model evaluation helpers
%-------------------------------------------------------------------------------
function v = sinSum(p,x,n)
    v = zeros(numel(x),1);
    for k = 1:n
        v = v + p(3*k-2)*sin(p(3*k-1)*x + p(3*k));
    end
end

function J = sinJac(p,x,n)
    J = zeros(numel(x),3*n);
    for k = 1:n
        a = p(3*k-2); b = p(3*k-1); c = p(3*k);
        s = sin(b*x + c); co = cos(b*x + c);
        J(:,3*k-2) = s;
        J(:,3*k-1) = a*x.*co;
        J(:,3*k)   = a*co;
    end
end

function v = gaussSum(p,x,n)
    v = zeros(numel(x),1);
    for k = 1:n
        v = v + p(3*k-2)*exp(-((x - p(3*k-1))/p(3*k)).^2);
    end
end

function J = gaussJac(p,x,n)
    J = zeros(numel(x),3*n);
    for k = 1:n
        a = p(3*k-2); b = p(3*k-1); c = p(3*k);
        u = (x - b)/c; e = exp(-u.^2);
        J(:,3*k-2) = e;
        J(:,3*k-1) = a*e.*(2*u/c);
        J(:,3*k)   = a*e.*(2*u.^2/c);
    end
end

%-------------------------------------------------------------------------------
% Start-point heuristics (used only when the caller supplies no StartPoint)
%-------------------------------------------------------------------------------
function p0 = expStart(x,y)
    % Log-linear regression on the positive part, falling back to a crude guess
    ok = y > 0;
    if sum(ok) >= 2
        c = [x(ok), ones(sum(ok),1)] \ log(y(ok));
        p0 = [exp(c(2)), c(1)];
    else
        p0 = [max(abs(y)), -1/max(eps,mean(x))];
    end
end

function p0 = powerStart(x,y)
    ok = x > 0 & y > 0;
    if sum(ok) >= 2
        c = [log(x(ok)), ones(sum(ok),1)] \ log(y(ok));
        p0 = [exp(c(2)), c(1)];
    else
        p0 = [1, 1];
    end
end

function p0 = sinStart(x,y,n)
    % Greedy global frequency sweep. For a fixed frequency b, the model
    % a*sin(b*x+c) = A*sin(b*x) + B*cos(b*x) is linear in (A,B), so the best
    % amplitude and phase come from a 2-column least-squares solve. Sweeping b
    % and keeping the best therefore locates the globally optimal component,
    % which is then peeled off and the sweep repeated.
    %
    % HCTSA2 2.0's rule (defect D1), now BF_fit's second start for sinN: fit()
    % starts from an FFT peak (sinFFTStart), which on noisy data can lead its
    % optimiser to a near-zero frequency where the model degenerates into a ramp.
    % BF_fit solves from both and keeps the better fit (D1).
    N = numel(x);
    span = max(x) - min(x);
    if span <= 0, span = max(N-1,1); end
    nGrid = min(4000, max(400, 4*N));
    b = logspace(log10(2*pi/(4*span)), log10(pi*(N-1)/span), nGrid);
    r = y;
    p0 = zeros(1,3*n);
    for k = 1:n
        bestSSE = Inf; bestB = b(1); bestC = [0;0];
        for j = 1:numel(b)
            D = [sin(b(j)*x), cos(b(j)*x)];
            co = D\r;
            s = sum((r - D*co).^2);
            if s < bestSSE
                bestSSE = s; bestB = b(j); bestC = co;
            end
        end
        p0(3*k-2) = hypot(bestC(1),bestC(2));
        p0(3*k-1) = bestB;
        p0(3*k)   = atan2(bestC(2),bestC(1));
        r = r - [sin(bestB*x), cos(bestB*x)]*bestC;
    end
end

function p0 = gaussStart(x,y,n)
    % Successive peaks of the (positive) signal, with a fixed width (HCTSA2 2.0's
    % rule; now the second start)
    p0 = zeros(1,3*n);
    r = y;
    w = (max(x) - min(x))/(2*n);
    for k = 1:n
        [a,idx] = max(r);
        p0(3*k-2) = a;
        p0(3*k-1) = x(idx);
        p0(3*k)   = max(w,eps);
        r = r - a*exp(-((x - x(idx))/max(w,eps)).^2);
    end
end

%-------------------------------------------------------------------------------
% fit()'s start rules, implemented from their description (the default start
% points of the Curve Fitting Toolbox library models). Where fit() would fall back
% to random numbers, these return a fixed guess instead.
%-------------------------------------------------------------------------------
function p0 = expSumsStart(x,y)
    % Geometric sums (Petersson & Holmstrom, "Initial values for the exponential
    % sum least squares fitting problem"). On an even grid of step h, a*exp(b*x)
    % is a geometric sequence of ratio exp(b*h), so two consecutive blocks of q
    % points sum in the ratio exp(b*h*q). Then a = sum(y)/sum(exp(b*x)).
    [x,o] = sort(x); y = y(o);
    n = numel(x);
    if any(diff(diff(x)) > 1/n^2)
        % Uneven x: drop the first of any (near-)repeated x, and resample onto an
        % even grid by exponential interpolation between neighbouring points.
        % Signs are carried over by position, as fit() does.
        keep = [diff(x) >= eps^0.7; true];
        x = x(keep); y = y(keep);
        n = numel(x);
        if n < 2
            p0 = [mean(y), 0];
            return
        end
        xg = linspace(x(1),x(end),n)';
        xg(end) = xg(end) - eps(xg(end));
        k = zeros(n,1); j = 1;               % xg(i) lies in [x(k(i)), x(k(i)+1))
        for i = 1:n
            while j < n-1 && x(j+1) <= xg(i), j = j + 1; end
            k(i) = j;
        end
        ly = log(abs(y) + eps);
        slope = (ly(k+1) - ly(k))./(x(k+1) - x(k));
        y = sign(y).*exp(ly(k) + slope.*(xg - x(k)));
        x = xg;
    end
    q = floor(n/2);
    S1 = sum(y(1:q));
    while S1 == 0                            % leading zeros: shorter blocks
        q = q - 1;
        if q < 1
            p0 = [0, 0];
            return
        end
        S1 = sum(y(1:q));
    end
    S2 = sum(y(q+1:2*q));
    h = mean(diff(x));
    if h > 0
        b = log(abs(S2/S1))/(q*h);
    else
        b = 1;                               % constant x: no information
    end
    a = sum(y)/sum(exp(b*x));
    if ~isfinite(a), a = 1; end
    p0 = [a, b];
end

function p0 = powerRatioStart(x,y)
    % Every point on one curve a*x^b: b is the mean over the other points of
    % log(y0/y)/log(x0/x), (x0,y0) being the point of smallest x (its repeats
    % averaged), and a the mean of y./x.^b. Only points with x > 0 and y > 0.
    ok = x > 0 & y > 0;
    x = x(ok); y = y(ok);
    if isempty(x)
        p0 = [1, 1];
        return
    end
    x0 = min(x);
    at = x == x0;
    y0 = mean(y(at));
    xr = x(~at); yr = y(~at);
    if isempty(xr)
        p0 = [1, 1];
        return
    end
    b = mean(log(y0./yr)./log(x0./xr));
    a = mean([y0; yr]./[x0; xr].^b);
    p0 = [a, b];
end

function p0 = gaussPeelStart(x,y,n)
    % One peak at a time: the height a and centre b of the highest point, and the
    % width from inverting the Gaussian at every point below it,
    % c = mean(|x-b|/sqrt(log(a/y))), shrunk by (2n - k) for the k-th peak found
    % so far; the peak is subtracted before the next. Peaks it cannot find (no
    % point strictly between 0 and the maximum) are spread evenly over x, with the
    % remaining maximum as height and their spacing as width.
    [x,o] = sort(x); y = y(o);
    A = zeros(1,n); B = zeros(1,n); C = zeros(1,n);
    k = 0;
    while k < n
        j = find(y == max(y),1,'last');
        a = y(j); b = x(j);
        below = y > 0 & y < a;
        if ~any(below), break; end
        c = mean(abs(x(below) - b)./sqrt(log(a./y(below))))/(2*n - k);
        y = y - a*exp(-((x - b)/c).^2);
        k = k + 1;
        A(k) = a; B(k) = b; C(k) = c;
    end
    if k < n
        g = linspace(min(x),max(x),n - k + 2);
        B(k+1:n) = g(2:end-1);
        A(k+1:n) = max(y);
        C(k+1:n) = g(2) - g(1);
    end
    p0 = reshape([A; B; C],1,[]);
end

function p0 = sinFFTStart(x,y,n)
    % Frequencies from FFT peaks, one component at a time (each new peak searched
    % with the earlier ones zeroed, on the residual of the fit so far); then, for
    % those frequencies, amplitudes and phases by linear least squares on
    % sin/cos pairs. x is sorted and near-repeated x are merged (averaged).
    [x,o] = sort(x); y = y(o);
    dup = [diff(x) < eps^0.7; false];
    nxt = [false; dup(1:end-1)];
    x(dup) = (x(dup) + x(nxt))/2;  x(nxt) = [];
    y(dup) = (y(dup) + y(nxt))/2;  y(nxt) = [];
    N = numel(x);
    if N < 2
        p0 = [ones(1,n); ones(1,n); zeros(1,n)];
        p0 = p0(:)';
        return
    end
    w = zeros(n,1); used = zeros(0,1);
    res = y;
    for j = 1:n
        F = fft(res);                        % no mean removed: no constant term
        F(used) = 0;
        [~,k] = max(abs(F(1:floor(N/2))));
        used(end+1,1) = k; %#ok<AGROW>
        w(j) = 2*pi*max(0.5,k-1)/(x(end) - x(1));
        X = zeros(N,2*j);
        for q = 1:j
            X(:,2*q-1) = sin(w(q)*x);
            X(:,2*q)   = cos(w(q)*x);
        end
        ab = X\y;
        res = y - X*ab;
    end
    p0 = zeros(1,3*n);
    for q = 1:n
        p0(3*q-2:3*q) = [hypot(ab(2*q-1),ab(2*q)), w(q), atan2(ab(2*q),ab(2*q-1))];
    end
end

%-------------------------------------------------------------------------------
function p = gaussCanonical(p,n)
    k = 3*(1:n);
    p(k) = abs(p(k));
end

function p = sinCanonical(p)
    % a*sin(b*x + c) = -a*sin(-b*x - c)
    for k = 1:3:numel(p)
        if p(k+1) < 0
            p(k:k+2) = -p(k:k+2);
        end
    end
end

%-------------------------------------------------------------------------------
function J = gaussPeakJac(p,x)
d = x - p(2);
E = exp(-d.^2/(2*p(3)^2));
J = [E, p(1)*E.*d/p(3)^2, p(1)*E.*d.^2/p(3)^3];
end

%-------------------------------------------------------------------------------
function J = kneeJac(p,x)
xc = x.^p(3);
den = (p(2) + xc)*log(10);
J = [ones(size(x)), -1./den, -xc.*log(x)./den];
end
