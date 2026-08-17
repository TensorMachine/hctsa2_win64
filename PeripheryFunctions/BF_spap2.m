function sp = BF_spap2(nPieces,k,x,y)
% BF_spap2   Toolbox-free replacement for the Curve Fitting Toolbox's spap2().
%
% Least-squares spline approximation of order K with NPIECES polynomial pieces.
% Matches the scalar-knots form of spap2, which is the only form hctsa uses.
%
% Knot placement follows spap2's aptknt rule, NOT uniform breaks:
%
%   sites          = linspace(min(x),max(x),nPieces+k-1)
%   interior knots = moving average of (k-1) consecutive sites
%   knots          = augknt([min(x), interior, max(x)], k)
%
% e.g. for k=4: 2 pieces -> [1/2];  4 pieces -> [1/3, 1/2, 2/3];
%               6 pieces -> [1/4, 3/8, 1/2, 5/8, 3/4]   (fractions of the span)
%
% Uniform breaks coincide with this only at nPieces=2, which is why that case
% agrees under either rule. Verified against the Bonn EEG reference values
% (real Curve Fitting Toolbox output) to ~1e-13 for 2, 4 and 6 pieces.
%
%---INPUTS:
% nPieces  number of polynomial pieces
% k        spline order (order = degree+1, so k = 4 is cubic)
% x, y     data, either orientation
%
%---OUTPUT:
% sp       struct with fields knots, coefs, order, number -- evaluate with BF_fnval

x = x(:);
y = y(:);

ok = isfinite(x) & isfinite(y);
x = x(ok);
y = y(ok);

if numel(x) < k + nPieces - 1
    error('curvefit:spap2:notEnoughData', ...
          'Not enough data (%u points) for %u B-spline coefficients.', ...
          numel(x),k+nPieces-1);
end

% aptknt rule: interior knots are moving averages of (k-1) consecutive sites
nSites = nPieces + k - 1;
sites  = linspace(min(x),max(x),nSites);
interior = zeros(1,nSites-k);
for j = 1:nSites-k
    interior(j) = mean(sites(j+1:j+k-1));
end
breaks = [min(x), interior, max(x)];
knots  = [repmat(breaks(1),1,k-1), breaks, repmat(breaks(end),1,k-1)];

B = BF_bsplineBasis(knots,k,x);

% Least squares; a rank-deficient basis (e.g. a gap in x spanning a whole piece)
% falls back to the minimum-norm solution rather than erroring
warnState = warning('off','all');
coefs = B\y;
warning(warnState);

if ~all(isfinite(coefs))
    error('curvefit:spap2:nanComputed','NaN computed in spline fit.');
end

sp.form   = 'B-';
sp.knots  = knots;
sp.coefs  = coefs(:)';
sp.number = numel(coefs);
sp.order  = k;
sp.dim    = 1;

end
