function T = BF_wpdec(x,depth,wname,entName)
% BF_wpdec   Full wavelet packet decomposition (toolbox-free wpdec()).
%
%   T = BF_wpdec(x,depth,wname)          % Shannon entropy
%   T = BF_wpdec(x,depth,wname,'shannon')
%
% Returns a plain struct in place of MATLAB's wptree object:
%   T.depth, T.wname
%   T.coefs{k+1}   coefficients of node k (node 0 = the signal; the children of
%                  node k are 2k+1 = approximation and 2k+2 = detail)
%   T.ent(k+1)     entropy of node k's coefficients
%   T.tn           terminal nodes, in left-to-right order (here: depth `depth`)
%
% Every node is split with the single-level DWT (BF_dwt, MATLAB's default 'sym'
% extension) down to `depth`, and each node's entropy is the legacy Shannon
% entropy -sum(c.^2 .* log(c.^2)) (BF_wentropy). Matches MATLAB's
% wpdec(x,depth,'db3','shannon') -- every node's coefficients and entropy --
% on the five series of the licensed run of 2026-09-23. See BF_besttree and
% BF_wenergy.

if nargin < 4 || isempty(entName), entName = 'shannon'; end
if ~strcmpi(entName,'shannon')
    error('BF_wpdec:entropy','Only Shannon entropy is supported (got ''%s'').',entName);
end

if isempty(x)
    error('BF_wpdec:emptyInput','The input signal is empty.');
end
if ~isnumeric(depth) || ~isscalar(depth) || ~isfinite(depth) || depth < 0 || depth ~= fix(depth)
    error('BF_wpdec:depth','The depth must be a non-negative integer.');
end
if depth > BF_wmaxlev(numel(x),wname)
    error('BF_wpdec:depth','Depth %d exceeds the maximum level %d for %d samples.', ...
          depth,BF_wmaxlev(numel(x),wname),numel(x));
end

nNodes = 2^(depth+1) - 1;
T.depth = depth;
T.wname = wname;
T.coefs = cell(nNodes,1);
T.coefs{1} = x(:);
for k = 0:(2^depth - 2)                      % every node above the bottom level
    [a,d] = BF_dwt(T.coefs{k+1},wname);
    T.coefs{2*k+2} = a;
    T.coefs{2*k+3} = d;
end
T.ent = cellfun(@(c) BF_wentropy(c,'shannon'),T.coefs);
T.tn = (2^depth - 1:nNodes - 1)';
end
