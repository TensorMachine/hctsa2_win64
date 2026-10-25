function [T,ento] = BF_besttree(T)
% BF_besttree   Best-basis wavelet packet tree (toolbox-free besttree()).
%
%   Tbest = BF_besttree(BF_wpdec(x,depth,wname))
%
% The Coifman-Wickerhauser best basis: working up from the deepest level, a
% node keeps its split only when the children's optimal entropies sum to LESS
% than the node's own entropy; otherwise (ties included) its subtree is merged
% into it. Returns the pruned tree (T.tn = the best basis' terminal nodes, in
% left-to-right order) and each node's optimal entropy. Matches MATLAB's
% besttree on the five series of the licensed run of 2026-09-23.

nNodes = numel(T.coefs);
ent = T.ent(:);
ento = ent;
merged = false(nNodes,1);
for k = (2^T.depth - 2):-1:0                 % internal nodes, deepest first
    echild = ento(2*k+2) + ento(2*k+3);
    if echild < ent(k+1)
        ento(k+1) = echild;
    else
        merged(k+1) = true;                  % keep the node, drop its subtree
    end
end

% Terminal nodes of the pruned tree, left to right
tn = zeros(0,1);
stack = 0;
while ~isempty(stack)
    k = stack(end); stack(end) = [];
    if merged(k+1) || k >= 2^T.depth - 1
        tn(end+1,1) = k; %#ok<AGROW>
    else
        stack = [stack, 2*k+2, 2*k+1];       % visit the approximation child first
    end
end
T.tn = tn;
T.ento = ento;
end
