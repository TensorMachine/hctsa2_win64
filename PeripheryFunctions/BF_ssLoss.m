function [V,x0,e] = BF_ssLoss(A,C,K,y)
% BF_ssLoss   Prediction-error loss of an innovations-form state-space model,
% with the initial state estimated.
%
% The one-step predictor is
%   x(t+1) = (A - K*C) x(t) + K y(t),   e(t) = y(t) - C x(t)
%
% The initial state enters linearly, so the response to each unit initial state
% is computed and x0 is found by least squares -- matching n4sid's default
% InitialState = 'estimate'.

y = y(:);
T = numel(y);
n = size(A,1);
Ac = A - K*C;

% Response with zero initial state
x = zeros(n,1);
e0 = zeros(T,1);
for t = 1:T
    e0(t) = y(t) - C*x;
    x = Ac*x + K*y(t);
end

% Response to each unit initial state (input-free)
% One free-response recursion per state. The matrix-form alternative (propagate
% all n together) is 8.95x faster but NOT bitwise identical for n > 1, so it
% lives on the ''performance'' branch rather than here.
B = zeros(T,n);
for j = 1:n
    x = zeros(n,1); x(j) = 1;
    for t = 1:T
        B(t,j) = -C*x;
        x = Ac*x;
    end
end

warnState = warning('off','all');
x0 = (B'*B)\(B'*e0);
warning(warnState);
if ~all(isfinite(x0))
    x0 = zeros(n,1);
end
e = e0 - B*x0;
x0 = -x0;                 % sign: B holds -C*x response
V = sum(e.^2)/T;
end
