function e = BF_wentropy(x,tname,p)
% BF_wentropy   Toolbox-free replacement for wentropy().
%
% Entropy of a signal, using the conventions 0*log(0)=0 and log(0)=0.
x = x(:);
if nargin < 3, p = 0; end
switch lower(strrep(tname,' ',''))
case 'shannon'
    x2 = x.^2;
    t = x2.*log(x2);
    t(x2 == 0) = 0;
    e = -sum(t);
case {'logenergy','log_energy'}
    x2 = x.^2;
    t = log(x2);
    t(x2 == 0) = 0;
    e = sum(t);
case 'threshold'
    e = sum(abs(x) > p);
case 'sure'
    e = numel(x) - 2*sum(abs(x) <= p) + sum(min(x.^2,p^2));
case 'norm'
    e = sum(abs(x).^p);
otherwise
    error('BF_wentropy:unknownType','Unknown entropy type ''%s''.',tname);
end
end
