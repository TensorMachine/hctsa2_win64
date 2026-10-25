function fdr = BF_mafdr(p,varargin)
% BF_mafdr   Toolbox-free replacement for the Bioinformatics Toolbox's mafdr(p,'BHFDR',true).
%
%   fdr = BF_mafdr(p,'BHFDR',true)
%
% Benjamini-Hochberg false discovery rates (adjusted p-values) of the p-values p:
% with the m p-values sorted, p_(1) <= ... <= p_(m),
%
%   fdr_(i) = min over k >= i of  min(1, m*p_(k)/k),
%
% returned in the order and shape of p (Benjamini & Hochberg 1995, the linear
% step-up procedure that mafdr's 'BHFDR' option documents). Tied p-values get the
% same value. Only this option is implemented: mafdr's default, Storey's
% q-values (pFDR with an estimated pi0), is not, and asking for it is an error.
%
% (HCTSA2 3.0.0) Upstream's TS_TopFeatures computes its q-values with mafdr when
% 'numNulls' > 0; this keeps that option within Statistics and Signal Processing.

if nargin ~= 3 || ~(ischar(varargin{1}) || isstring(varargin{1})) || ...
        ~strcmpi(varargin{1},'BHFDR') || ~isequal(varargin{2},true)
    error('BF_mafdr:unsupported', ...
        'BF_mafdr implements only mafdr(p,''BHFDR'',true) (Benjamini-Hochberg).');
end
if ~isnumeric(p) || ~isreal(p) || any(~isfinite(p(:))) || any(p(:) < 0 | p(:) > 1)
    error('BF_mafdr:badPValues','p-values must be finite and in [0,1].');
end

fdr = zeros(size(p));
m = numel(p);
if m == 0
    return
end
[ps,order] = sort(p(:));
q = min(1, ps*m./(1:m)');
q = flipud(cummin(flipud(q)));          % the smallest value at or above each rank
fdr(order) = q;

end
