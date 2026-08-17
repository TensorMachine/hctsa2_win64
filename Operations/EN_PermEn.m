function out = EN_PermEn(y,m,tau)
% EN_PermEn     Permutation Entropy of a time series.
%
% "Permutation Entropy: A Natural Complexity Measure for Time Series"
% C. Bandt and B. Pompe, Phys. Rev. Lett. 88(17) 174102 (2002)
%
%---INPUTS:
% y, the input time series
% m, the embedding dimension (or order of the permutation entropy)
% tau, the time-delay for the embedding
%
%---OUTPUT:
% Outputs the permutation entropy and normalized version computed according to
% different implementations

% ------------------------------------------------------------------------------
% Copyright (C) 2020, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite the following two papers:
%
% (1) B.D. Fulcher and N.S. Jones, "hctsa: A Computational Framework for Automated
% Time-Series Phenotyping Using Massive Feature Extraction, Cell Systems 5: 527 (2017).
% DOI: 10.1016/j.cels.2017.10.001
%
% (2) B.D. Fulcher, M.A. Little, N.S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013).
% DOI: 10.1098/rsif.2013.0048
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program. If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Check inputs and set defaults:
%-------------------------------------------------------------------------------

if nargin < 2 || isempty(m)
    m = 2; % order 2
end

if nargin < 3 || isempty(tau)
    tau = 1;
end

% ------------------------------------------------------------------------------
% Embed the signal:
% ------------------------------------------------------------------------------
x = BF_Embed(y,tau,m,false);
Nx = size(x,1); % number of embedding vectors produced
if Nx < 5 % need at least 5 embedding vectors to actually do a computation
    error('Time series too short to embed');
end
% Generate permutations up to the embedding dimension, m:
permList = perms(1:m);
numPerms = length(permList);

% Initialize
countPerms = zeros(numPerms,1);

% Count each type of permutation through the time series
% Map each window's permutation to its index by table lookup instead of a
% linear scan. The original compared the sorted order against every row of
% permList until it matched -- up to factorial(m) comparisons per window, so
% ~490,000 for m = 5 on a 4097-point series.
%
% A permutation of 1..m is uniquely identified by sum(perm .* m.^(0:m-1)),
% which indexes a small dense table. The counts are integers and the mapping
% is exact, so countPerms is bit-for-bit what the scan produced.
powers = m.^(0:m-1);
permKey = permList*powers(:);
lookup = zeros(max(permKey),1);
lookup(permKey) = 1:numPerms;

[~,ixAll] = sort(x,2);
keys = ixAll*powers(:);
idx = lookup(keys);
for j = 1:Nx
    countPerms(idx(j)) = countPerms(idx(j)) + 1;
end

p = countPerms/Nx; %((Nx-(m-1))*tau);
p_0 = p(p>0); % makes log(0) = 0
out.permEn = -sum(p_0.*log2(p_0));

% Normalized permutation entropy (more comparable across m?)
mFact = factorial(m);
out.normPermEn = out.permEn/log2(mFact);

% ------------------------------------------------------------------------------
% Adapted implementation by Bruce Land and Damian Elias:
% cf.:
% http://people.ece.cornell.edu/land/PROJECTS/Complexity/
% http://people.ece.cornell.edu/land/PROJECTS/Complexity/logisticPE.m

% Not clear to me why you would make log(0) = log(1/Nx); (the minimum)
% rather than exclude it from the sum, as is done here:
p_LE = max(1/Nx,p);
out.permEnLE = -sum(p_LE.*log(p_LE))/(m-1);

end
