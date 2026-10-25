function e = BF_wentropy(x,varargin)
% BF_wentropy   Toolbox-free replacement for wentropy().
%
% Two interfaces, as in MATLAB:
%
%   e = BF_wentropy(x,'shannon')            legacy (R2022a and before): the
%   e = BF_wentropy(x,'logenergy')          entropy of the values themselves
%
%   e = BF_wentropy(x,'Wavelet','sym4','Distribution','global','Level',5)
%                                           current: wavelet entropy
%
%---LEGACY
% Entropy of a signal, using the conventions 0*log(0)=0 and log(0)=0. Only
% 'shannon' (used by BF_wpdec) and 'logenergy' are provided; MATLAB's
% 'threshold', 'sure' and 'norm' were never validated here and raise an error.
%
%---CURRENT (name-value pairs)
% The Shannon wavelet entropy of MATLAB's wentropy (R2022b onwards) with the
% options hctsa uses:
%   'Wavelet'          default 'sym4'
%   'Level'            default floor(log2(N)) - 1
%   'Distribution'     'scale' (default) or 'global'
%   'Scaled'           true (default): divide by log of the number of outcomes
%   'Transform'        'modwt' only (the default)
%   'Entropy'          'Shannon' only (the default)
%   'EnergyThreshold'  default 1e-8
% The series is decomposed with the MODWT (BF_modwt). 'global': the outcomes are
% the levels (wavelet levels plus the final scaling level), with probability
% each level's share of the total energy. 'scale': one entropy per level, the
% outcomes being the times, with probability each coefficient's share of that
% level's energy (returned as a column). A total or level energy below the
% threshold counts as zero, giving entropy 0. MATLAB's time alignment of the
% MODWT is a circular shift within each level, which neither distribution sees.
%
% Matches MATLAB's wentropy for 'global' and 'scale', scaled and unscaled, with
% a given and a default level, to 1e-12 on the eight series of the licensed run
% of 2026-09-23. A NaN or Inf in the input is an error ('BF_wentropy:nonFinite')
% rather than an entropy of 0.

if nargin < 2
    error('BF_wentropy:nargin','Specify an entropy type or name-value options.');
end
for i = 1:numel(varargin)              % string scalars as char (HCTSA2 2.2.1)
    if isstring(varargin{i}) && isscalar(varargin{i}), varargin{i} = char(varargin{i}); end
end
if ischar(varargin{1}) && any(strcmpi(varargin{1}, ...
        {'Wavelet','Transform','Entropy','Exponent','Scaled','Level', ...
         'Distribution','EnergyThreshold'}))
    e = waveletEntropy(x,varargin{:});
    return
end

%-------------------------------------------------------------------------------
%% Legacy interface
tname = varargin{1};
x = x(:);
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
case {'threshold','sure','norm'}
    error('BF_wentropy:unvalidated','Legacy entropy ''%s'' is not validated in HCTSA2.',tname);
otherwise
    error('BF_wentropy:unknownType','Unknown entropy type ''%s''.',tname);
end
end

%-------------------------------------------------------------------------------
function e = waveletEntropy(x,varargin)
if mod(numel(varargin),2) ~= 0
    error('BF_wentropy:option','Options must be name-value pairs.');
end
opts = struct('wavelet','sym4','transform','modwt','entropy','shannon', ...
    'scaled',true,'level',[],'distribution','scale','energythreshold',1e-8);
for i = 1:2:numel(varargin)-1
    name = lower(varargin{i});
    if ~isfield(opts,name) || strcmp(name,'exponent')
        error('BF_wentropy:option','Option ''%s'' is not supported.',varargin{i});
    end
    opts.(name) = varargin{i+1};
end
if ~strcmpi(opts.transform,'modwt')
    error('BF_wentropy:transform','Only the MODWT transform is supported (got ''%s'').',opts.transform);
end
if ~strcmpi(opts.entropy,'shannon')
    error('BF_wentropy:entropy','Only Shannon entropy is supported (got ''%s'').',opts.entropy);
end

x = x(:);
N = numel(x);
if N < 4
    error('BF_wentropy:tooShort','The input must have at least four samples.');
end
if any(~isfinite(x))
    error('BF_wentropy:nonFinite','The input contains NaN or Inf values.');
end
level = opts.level;
if isempty(level)
    level = floor(log2(N)) - 1;
end
w = BF_modwt(x,opts.wavelet,level);       % (level+1) x N
thresh = opts.energythreshold;

switch lower(opts.distribution)
    case 'global'
        E = sum(w.^2,2);
        total = sum(E);
        if total < thresh
            P = zeros(size(E));
        else
            P = E/total;
        end
        n = numel(P);
        e = -sumOmitNaN(P.*log(P));
        e(e < eps) = 0;
    case 'scale'
        E = sum(w.^2,2);
        zeroE = E < thresh;
        w(zeroE,:) = 0;
        E(zeroE) = 1;
        P = (w.^2)./E;
        n = N;
        e = -sumOmitNaN((P.*log(P))')';   % one value per level
        e(e < eps) = 0;
    otherwise
        error('BF_wentropy:distribution','Unknown distribution ''%s''.',opts.distribution);
end
e(~isfinite(e)) = 0;
if opts.scaled
    if n == 1
        e(:) = 0;
    else
        e = e/log(n);
    end
end
end

function s = sumOmitNaN(v)
% Column sums ignoring NaN (0*log(0) is NaN and counts as 0).
v(isnan(v)) = 0;
s = sum(v,1);
end
