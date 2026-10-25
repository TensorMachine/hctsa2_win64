function d = BF_detcoef(C,L,N)
% BF_detcoef   Extract level-N detail coefficients (toolbox-free detcoef()).
% L = [len(cA_J), len(cD_J), ..., len(cD_1), len(x)]
J = numel(L) - 2;
if numel(N) > 1                        % several levels: a cell array, as detcoef
    d = cell(1,numel(N));
    for i = 1:numel(N), d{i} = BF_detcoef(C,L,N(i)); end
    return
end
if ~isnumeric(N) || ~isscalar(N) || N ~= fix(N) || N < 1 || N > J
    error('BF_detcoef:badLevel','Level %d outside 1..%d.',N,J);
end
idx = J - N + 2;
first = sum(L(1:idx-1)) + 1;
d = C(first:first+L(idx)-1);
end
