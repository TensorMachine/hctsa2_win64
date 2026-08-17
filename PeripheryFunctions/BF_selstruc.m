function [nn,vmod] = BF_selstruc(V,c)
% BF_selstruc   Select a model structure from an arxstruc loss table.
%
% Toolbox-free replacement for selstruc().
%
%---INPUTS:
% V   loss table from BF_arxstruc
% c   0      -> minimise the loss function itself (default)
%     'aic'  -> minimise log(V) + 2*d/N
%     'mdl'  -> minimise log(V) + d*log(N)/N
%     numeric -> minimise V*(1 + c*d/N)
%
%---OUTPUTS:
% nn    the selected structure (the order, for the AR case hctsa uses)
% vmod  the criterion values for every structure

if nargin < 2 || isempty(c), c = 0; end

m = size(V,2) - 1;
v = V(1,1:m);
d = V(2,1:m);          % number of estimated parameters per structure
N = V(1,m+1);          % number of validation samples

if ischar(c)
    switch lower(c)
    case 'aic'
        vmod = log(v) + 2*d/N;
    case 'mdl'
        vmod = log(v) + d.*log(N)/N;
    otherwise
        error('BF_selstruc:unknownCriterion','Unknown criterion ''%s''.',c);
    end
elseif c == 0
    vmod = v;
else
    vmod = v.*(1 + c*d/N);
end

[~,ix] = min(vmod);
nn = V(2:end,ix)';
end
