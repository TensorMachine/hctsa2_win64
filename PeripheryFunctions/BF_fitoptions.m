function s = BF_fitoptions(varargin)
% BF_fitoptions   Toolbox-free stand-in for fitoptions().
%
% hctsa only ever calls this as
%   fitoptions('Method','NonlinearLeastSquares','StartPoint',p0)
% so only Method and StartPoint are carried. Anything else is stored but unused.

s = struct('Method','NonlinearLeastSquares','StartPoint',[]);
for k = 1:2:numel(varargin)-1
    s.(varargin{k}) = varargin{k+1};
end
end
