function s = BF_fitoptions(varargin)
% BF_fitoptions   Toolbox-free stand-in for fitoptions().
%
% hctsa calls this as
%   fitoptions('Method','NonlinearLeastSquares','StartPoint',p0)
% and, in SP_Specparam, with 'Lower', 'Upper' (bounds on the coefficients, used
% by BF_fit), 'MaxIter' and 'Display'.
%
% As fitoptions does, names are matched case-insensitively and stored with
% fitoptions' own spelling; an unknown name is an error. (HCTSA2 2.2.1: names
% were stored as typed, so 'startpoint' was silently ignored by BF_fit.) BF_fit
% uses Method, StartPoint, Lower and Upper, ignores MaxIter and Display, and
% refuses the other fitoptions properties rather than dropping them.

s = struct('Method','NonlinearLeastSquares','StartPoint',[],'Lower',[],'Upper',[]);
known = {'Method','StartPoint','Lower','Upper','MaxIter','Display', ...
         'Weights','Exclude','Robust','Normalize','Algorithm','TolFun','TolX', ...
         'MaxFunEvals','DiffMinChange','DiffMaxChange'};
if mod(numel(varargin),2) ~= 0
    error('BF_fitoptions:badOptions','Options must be name/value pairs.');
end
for k = 1:2:numel(varargin)-1
    i = find(strcmpi(varargin{k},known),1);
    if isempty(i)
        error('BF_fitoptions:unknownOption','Unknown option ''%s''.',char(varargin{k}));
    end
    s.(known{i}) = varargin{k+1};
end
end
