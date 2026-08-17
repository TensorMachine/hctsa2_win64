function f = BF_fittype(model,varargin)
% BF_fittype   Toolbox-free stand-in for fittype().
%
% Returns a struct carrying the model string and any options, for BF_fit.
% The 'independent' argument is accepted and ignored: every hctsa model uses x.

f = struct('model',model,'options',BF_fitoptions());
for k = 1:2:numel(varargin)-1
    if strcmpi(varargin{k},'options')
        f.options = varargin{k+1};
    end
end
% Validate the model now, so a typo errors here rather than at fit time
BF_fitModel(model);
end
