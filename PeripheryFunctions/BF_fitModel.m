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
%
% Covers every model form used anywhere in hctsa. See BF_fit.

name = strtrim(name);
x2c = @(x) x(:); % force column

% --- polyN: linear, coefficients p1..p(N+1), highest power first ---------------
if numel(name) == 5 && strcmp(name(1:4),'poly') && any(name(5) == '123456789')
    n = str2double(name(5));
    m.coeffs = arrayfun(@(k) sprintf('p%d',k), 1:n+1, 'UniformOutput', false);
    m.linear = true;
    m.basis  = @(x) bsxfun(@power, x2c(x), n:-1:0);
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
    m.start  = @(x,y) sinStart(x2c(x),x2c(y),n);
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
    m.start  = @(x,y) gaussStart(x2c(x),x2c(y),n);
    return
end

% --- named / custom model strings ---------------------------------------------
% Whitespace is stripped so 'a*x +b' and 'a*x+b' are the same model.
switch regexprep(name,'\s','')
case {'a*x+b'}
    m.coeffs = {'a','b'};  m.linear = true;
    m.basis  = @(x) [x2c(x), ones(numel(x),1)];

case {'a*exp(b*x)','exp1'}
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
    m.start  = @(x,y) powerStart(x2c(x),x2c(y));

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
    % This is deliberately NOT what the Curve Fitting Toolbox does: its sinN fit
    % starts from an FFT peak and lets the optimiser wander, which on noisy data
    % often slides to a near-zero frequency where the model degenerates into a
    % ramp. That behaviour is path-dependent and not reproducible outside their
    % optimiser. See HCTSA2_Status.md (defect D1) for the measured comparison.
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
    % Successive peaks of the (positive) signal
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
