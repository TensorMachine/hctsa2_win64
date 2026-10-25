classdef BF_MutualInfoCalculator < handle
% BF_MutualInfoCalculator   Mutual information calculator with JIDT's interface.
%
% HCTSA2 2.4.0. IN_Initialize_MI returns one of these where it used to return a
% calculator from the bundled JIDT (infodynamics.jar), so its callers keep their
% calls:
%
%   miCalc = BF_MutualInfoCalculator('kraskov1',4);
%   miCalc.initialise(1,1);              % source and target dimensions
%   miCalc.setObservations(y1,y2);       % one observation per row
%   mi = miCalc.computeAverageLocalOfObservations();     % nats
%
% Methods:
%   'kraskov1', 'kraskov2'  the Kraskov-Stoegbauer-Grassberger estimator, algorithm
%                           1 or 2, with k nearest neighbours (KSG_mi MEX: build it
%                           with HCTSA2_BuildBinaries). Same conventions as JIDT 1.3
%                           as hctsa used it: columns normalised, max norm, no noise.
%   'gaussian'              assumes joint normality: 0.5*log(det(Rx)*det(Ry)/det(R)) from
%                           the correlation matrices (for two variables -0.5*log(1-r^2));
%                           det(R) is floored at eps, so identical or collinear
%                           variables give a large finite value, not Inf or a complex one
%
% JIDT's kernel estimator is not provided (no registered feature used it).

properties (SetAccess = private)
    method              % 'kraskov1', 'kraskov2' or 'gaussian'
    k = 4               % nearest neighbours (Kraskov)
    sourceDims = 1
    targetDims = 1
    source = []         % observations, one per row
    target = []
end

methods
    function obj = BF_MutualInfoCalculator(method,k)
        if ~ismember(method,{'kraskov1','kraskov2','gaussian'})
            error('BF_MutualInfoCalculator:method', ...
                'Unknown mutual information estimation method ''%s''.',method);
        end
        obj.method = method;
        if nargin > 1 && ~isempty(k)
            obj.k = k;
        end
    end

    function initialise(obj,sourceDims,targetDims)
        % Set the dimensions and clear the observations.
        obj.sourceDims = sourceDims;
        obj.targetDims = targetDims;
        obj.source = [];
        obj.target = [];
    end

    function setObservations(obj,source,target)
        % One observation per row; a vector is one column.
        if isvector(source) && obj.sourceDims == 1, source = source(:); end
        if isvector(target) && obj.targetDims == 1, target = target(:); end
        if size(source,2) ~= obj.sourceDims || size(target,2) ~= obj.targetDims
            error('BF_MutualInfoCalculator:dims', ...
                'Observations have %u and %u columns; initialised for %u and %u.', ...
                size(source,2),size(target,2),obj.sourceDims,obj.targetDims);
        end
        if size(source,1) ~= size(target,1)
            error('BF_MutualInfoCalculator:rows', ...
                'Source and target have %u and %u observations.',size(source,1),size(target,1));
        end
        obj.source = double(source);
        obj.target = double(target);
    end

    function mi = computeAverageLocalOfObservations(obj)
        switch obj.method
            case {'kraskov1','kraskov2'}
                if obj.targetDims ~= 1
                    error('BF_MutualInfoCalculator:dims','The Kraskov estimator takes a one-column target.');
                end
                mi = KSG_mi(obj.source,obj.target,obj.k,1 + strcmp(obj.method,'kraskov2'));
            case 'gaussian'
                X = [obj.source,obj.target];
                if any(std(X,0,1) == 0)          % a constant variable (JIDT stopped)
                    mi = NaN; return
                end
                R = corrcoef(X);
                s = 1:obj.sourceDims;
                t = obj.sourceDims + (1:obj.targetDims);
                mi = 0.5*log(det(R(s,s))*det(R(t,t))/max(det(R),eps));
        end
    end
end
end
