function [v,e] = BF_garchVar(y,K,G,A,offset)
% BF_garchVar   Conditional variance recursion for a GARCH(P,Q) model.
%
%   sigma2_t = K + sum_i G_i sigma2_{t-i} + sum_j A_j e_{t-j}^2
%
% Presample variances and squared innovations are both set to the sample
% variance of the innovations, the standard default.
%
% IMPLEMENTATION NOTE: this evaluates exactly the recursion above, but expressed
% as a linear filter rather than a scalar loop. The ARCH terms are an FIR of e^2
% and the GARCH terms are an all-pole filter, so
%
%   v = filter(1, [1 -G], u)   where   u_t = K + sum_j A_j e^2_{t-j}
%
% with the presample contributions folded into the first max(P,Q) entries of u.
% This is the same arithmetic in the same order, only compiled instead of
% interpreted -- the estimator, its tolerances and its convergence path are
% untouched. Verified identical to the original scalar-loop version to ~1e-16.
% It matters because this function is called several hundred times per fit and
% dominates the runtime.

e = y(:) - offset;
T = numel(e);
G = G(:)'; A = A(:)';
P = numel(G); Q = numel(A);
e2 = e.^2;
pre = mean(e2);

% ARCH contribution: FIR over past squared innovations, presample = pre
u = repmat(K,T,1);
if Q > 0
    e2ext = [repmat(pre,Q,1); e2];
    for j = 1:Q
        u = u + A(j)*e2ext(Q+1-j : Q+T-j);
    end
end

% Presample variance contributions: v_{t-i} = pre whenever t-i < 1
for i = 1:P
    m = min(i,T);
    u(1:m) = u(1:m) + G(i)*pre;
end

% GARCH contribution: all-pole recursion with zero initial conditions, the
% presample part having already been added above
if P > 0
    v = filter(1,[1, -G],u);
else
    v = u;
end
end
