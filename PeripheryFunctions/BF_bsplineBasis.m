function B = BF_bsplineBasis(t,k,x)
% BF_bsplineBasis   B-spline collocation matrix by the Cox-de Boor recursion.
%
% t  knot vector (nondecreasing), k  order (degree+1), x  evaluation points
% B  numel(x)-by-(numel(t)-k) matrix of basis function values

x = x(:);
m = numel(x);

% Order 1: indicator of each knot interval
N = zeros(m,numel(t)-1);
for i = 1:numel(t)-1
    N(:,i) = (x >= t(i)) & (x < t(i+1));
end
% Points at (or beyond) the final knot belong to the last interval of positive
% width, so the basis sums to 1 there rather than 0. (HCTSA2 2.2.1: points
% before the first knot likewise belong to the first one, k, so the spline
% extrapolates at both ends as fnval does, instead of evaluating to 0 on the left.)
N(x >= t(end),numel(t)-k) = 1;
N(x < t(1),k) = 1;

% Raise the order one at a time; terms over a zero-width knot span drop out
for d = 2:k
    Nd = zeros(m,numel(t)-d);
    for i = 1:numel(t)-d
        w1 = t(i+d-1) - t(i);
        w2 = t(i+d) - t(i+1);
        v = zeros(m,1);
        if w1 > 0, v = v + (x - t(i))/w1 .* N(:,i); end
        if w2 > 0, v = v + (t(i+d) - x)/w2 .* N(:,i+1); end
        Nd(:,i) = v;
    end
    N = Nd;
end
B = N;
end
