function f = BF_fittype(model,varargin)
% BF_fittype   Toolbox-free stand-in for fittype().
%
% Returns a struct carrying the model string and any options, for BF_fit.
% The 'independent' argument is accepted and ignored: every hctsa model uses x.
% 'coefficients', when given, must name the model's coefficients in the order
% BF_fitModel has them -- the order StartPoint, Lower and Upper refer to.

f = struct('model',model,'options',BF_fitoptions());
coeffs = {};
if mod(numel(varargin),2) ~= 0
    error('BF_fittype:badOptions','Options must be name/value pairs.');
end
for k = 1:2:numel(varargin)-1
    switch lower(varargin{k})
        case 'options'
            f.options = varargin{k+1};
        case 'coefficients'
            coeffs = varargin{k+1};
        case 'independent'
            % accepted and ignored, as above
        otherwise
            error('BF_fittype:unknownOption','Unknown or unsupported option ''%s''.',char(varargin{k}));
    end
end
% Validate the model now, so a typo errors here rather than at fit time
m = BF_fitModel(model);
if ~isempty(coeffs) && ~isequal(cellstr(coeffs),m.coeffs)
    error('BF_fittype:coefficients', ...
        'Coefficients {%s} do not match the model''s order {%s}.', ...
        strjoin(cellstr(coeffs),','),strjoin(m.coeffs,','));
end
end
