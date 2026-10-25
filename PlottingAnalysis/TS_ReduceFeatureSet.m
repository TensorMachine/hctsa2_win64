function [reducedIDs,clusterInfo] = TS_ReduceFeatureSet(whatData,distThreshold,useSpearman,doPlot)
% TS_ReduceFeatureSet   Unsupervised clustering to a reduced feature set
%
% Clusters the features by complete linkage on absolute-correlation distances,
% 1 - |r|, cut at distThreshold: every pair of features within a cluster then has
% |r| >= 1 - distThreshold. Each cluster is represented by its most central
% feature (the smallest summed distance to the cluster's other members).
%
%---INPUTS:
% whatData, the hctsa data (default 'norm'; anything TS_LoadData accepts)
% distThreshold, where to cut the dendrogram (default 0.2, i.e. |r| >= 0.8)
% useSpearman, Spearman (true, default) or Pearson (false) correlations
% doPlot, draw BF_ClusterDown's dendrogram and reordered correlation matrix
%       (default true; it holds several copies of the full matrix)
% ([] for an input gives its default.)
%
%---OUTPUTS:
% reducedIDs, a set of Operation IDs that forms the new reduced feature set: the
%       central feature of each cluster, largest cluster first
% clusterInfo, a table with one row per feature of whatData (its order): ID, Name,
%       Cluster (1 = the largest), ClusterSize, Central (true for the feature in
%       reducedIDs)
%
% A correlation is computed over the series where both features have a value. A
% pair that shares no such series, or where a feature is constant on them, counts
% as uncorrelated (distance 1). A single feature is one cluster.
%
% HCTSA2 2.5.1: help text and defaults for []; NaN-safe correlations (pdist made
% every distance involving a NaN NaN, which linkage cannot use); a distance matrix,
% not a vector, for BF_ClusterDown (upstream issue P37: it indexed the vector as a
% matrix and stopped); the doPlot input; the clusterInfo output; a single feature.

% ------------------------------------------------------------------------------
% Copyright (C) 2013-2026, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite the following two papers:
%
% (1) B.D. Fulcher and N.S. Jones, "hctsa: A Computational Framework for Automated
% Time-Series Phenotyping Using Massive Feature Extraction, Cell Systems 5: 527 (2017).
% DOI: 10.1016/j.cels.2017.10.001
%
% (2) B.D. Fulcher, M.A. Little, N.S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013).
% DOI: 10.1098/rsif.2013.0048
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Check inputs:
%-------------------------------------------------------------------------------
if nargin < 1 || isempty(whatData)
    whatData = 'norm';
end
if nargin < 2 || isempty(distThreshold)
    distThreshold = 0.2;
end
if nargin < 3 || isempty(useSpearman)
    useSpearman = true;
end
if nargin < 4 || isempty(doPlot)
    doPlot = true;
end

%-------------------------------------------------------------------------------
% Load in data:
[TS_DataMat,~,Operations] = TS_LoadData(whatData);

% Compute correlation distances between all pairs of features on the data
% (HCTSA2: over the series both features have, for features with missing values):
if useSpearman
    corrType = 'Spearman';
else
    corrType = 'Pearson';
end
if size(TS_DataMat,2) < 2
    % (HCTSA2: BF_ClusterDown reads a 1x1 distance matrix as a distance vector)
    cluster_Groupi = {1};
else
    R = corr(TS_DataMat,'Type',corrType);
    hasNaN = any(isnan(TS_DataMat),1);
    if any(hasNaN)
        R(hasNaN,:) = corr(TS_DataMat(:,hasNaN),TS_DataMat,'Type',corrType,'Rows','pairwise');
        R(:,hasNaN) = R(hasNaN,:)';
    end
    R(isnan(R)) = 0; % undefined (no shared series, or constant on them): uncorrelated
    R = min(max(R,-1),1); % (corr can return 1 + eps for identical features)
    distMat = 1 - (R + R')/2; % correlation distances, exactly symmetric
    clear('R');
    distMat(1:size(distMat,1)+1:end) = 0;

    % Cluster down with abscorr linkage clustering:
    [~,cluster_Groupi] = BF_ClusterDown(distMat,...
                    'clusterThreshold',distThreshold,...
                    'whatDistance','abscorr',...
                    'linkageMeth','complete',...
                    'doPlot',doPlot);
    clear('distMat');
end

% Determine cluster centers:
clusterCenterInd = cellfun(@(x)x(1),cluster_Groupi);

% Convert to feature IDs:
reducedIDs = Operations.ID(clusterCenterInd);

% Each feature's cluster:
numFeatures = height(Operations);
cluster = zeros(numFeatures,1);
clusterSize = zeros(numFeatures,1);
for j = 1:numel(cluster_Groupi)
    cluster(cluster_Groupi{j}) = j;
    clusterSize(cluster_Groupi{j}) = numel(cluster_Groupi{j});
end
central = false(numFeatures,1);
central(clusterCenterInd) = true;
clusterInfo = table(Operations.ID,Operations.Name,cluster,clusterSize,central, ...
                'VariableNames',{'ID','Name','Cluster','ClusterSize','Central'});

% Give user feedback:
fprintf(1,'Reduced to set of %u features using complete linkage clustering on abscorr distances (%.2f)\n',...
                length(reducedIDs),distThreshold);

end
