function E = BF_wenergy(T)
% BF_wenergy   Energy of each terminal node of a wavelet packet tree
% (toolbox-free wenergy() for wptree objects).
%
%   E = BF_wenergy(BF_besttree(BF_wpdec(x,depth,wname)))
%
% Percentage of the total energy of the terminal nodes' coefficients carried
% by each terminal node, in ascending node-index order (as MATLAB's wenergy,
% which reads the leaves sorted -- not in T.tn's left-to-right order). Matches
% MATLAB's wenergy on the five series of the licensed run of 2026-09-23.

tn = sort(T.tn);
e = cellfun(@(c) sum(c(:).^2),T.coefs(tn + 1));
E = 100*e(:)'/sum(e);
end
