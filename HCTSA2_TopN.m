function T = HCTSA2_TopN(results,N,varargin)
% HCTSA2_TopN   Narrow SeriesTable down to the N best-separating features.
%
%   T = HCTSA2_TopN(results,20)
%   T = HCTSA2_TopN(results,20,'Metadata',{'Name','Group','Length'})
%   T = HCTSA2_TopN('hctsa2_run_.../HCTSA2_Results.mat',50)
%
% RESULTS is either the struct returned by HCTSA2_Analyse or the path to the
% HCTSA2_Results.mat it wrote.
%
% Returns the metadata columns followed by the N features with the best
% separability rank, in rank order — a table you can read on screen, unlike the
% full SeriesTable (one column per feature: up to 6,392 in 3.0.0).
%
% Features with no rank (the ranking did not score them, e.g. because they
% failed to compute) are excluded rather than silently filling the list.

if nargin < 2 || isempty(N), N = 20; end

p = inputParser;
addParameter(p,'Metadata',{'Name','Group'},@iscell);
parse(p,varargin{:});
meta = p.Results.Metadata;

if ischar(results) || isstring(results)
    D = load(char(results));
    if ~isfield(D,'FeatureTable') || ~isfield(D,'SeriesTable')
        error('HCTSA2_TopN:badFile','%s has no FeatureTable/SeriesTable.',results);
    end
    F = D.FeatureTable; S = D.SeriesTable;
else
    F = results.FeatureTable; S = results.SeriesTable;
end

if ~any(strcmp(F.Properties.VariableNames,'Rank'))
    error('HCTSA2_TopN:noRank', ...
        ['FeatureTable has no Rank column -- the ranking did not run. ' ...
         'Check the analysis log.']);
end

ranked = F(isfinite(F.Rank),:);
if isempty(ranked)
    error('HCTSA2_TopN:noRankedFeatures', ...
        'No features carry a separability rank; nothing to select.');
end
ranked = sortrows(ranked,'Rank','ascend');
N = min(N,height(ranked));

% Each feature's SeriesTable column is recorded in FeatureTable.Column (2.2.1+).
% Results from 2.2.0 and earlier lack it: there the names are re-derived, which
% maps the 66 MF_GP_* features whose names are longer than 63 characters to the
% wrong columns (re-run HCTSA2_Analyse with 'ExistingHCTSA' to fix those files).
have = string(S.Properties.VariableNames);
if any(strcmp(F.Properties.VariableNames,'Column'))
    wanted = string(ranked.Column(1:N));
else
    warning('HCTSA2_TopN:oldResults',['These results predate FeatureTable.Column; ' ...
        'features with names over 63 characters may map to the wrong column.']);
    wanted = matlab.lang.makeValidName(string(ranked.Name(1:N)));
end
[tf,loc] = ismember(wanted,have);
if ~all(tf)
    missing = sum(~tf);
    warning('HCTSA2_TopN:missingColumns', ...
        ['%u of the top %u features are not columns of SeriesTable ' ...
         '(it was built with MaxFeatures).'],missing,N);
    wanted = wanted(tf); loc = loc(tf);
end

meta = meta(ismember(meta,S.Properties.VariableNames));
T = [S(:,meta), S(:,loc)];
end
