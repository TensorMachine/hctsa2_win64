function results = HCTSA2_Analyse(dataFile,varargin)
% HCTSA2_Analyse   End-to-end hctsa 2.0 pipeline: data in, two tables and plots out.
%
% This script is a little jank because it was written by Claude but it
% gives you the general flow of how to run HCTSA.
%
%   results = HCTSA2_Analyse('mydata.mat')
%   results = HCTSA2_Analyse('mydata.mat','Groups',g,'OutputDir','run1')
%
% Takes a .mat containing a time-series data matrix and class labels, runs the
% whole hctsa pipeline, and returns two consolidated tables:
%
%   results.FeatureTable  one row per FEATURE      -- how well it separates the
%                         classes, its rank, and its per-group means
%   results.SeriesTable   one row per TIME SERIES  -- its group and its feature
%                         values, ready for downstream modelling
%
% hctsa already provides the analysis (TS_TopFeatures, TS_Classify and friends).
% This wrapper orchestrates them, keeps going when individual steps fail, and
% builds the two tables -- which hctsa does not produce, scattering results
% across HCTSA.mat, figure windows and function return values instead.
%
%---INPUT
% dataFile   a .mat file. The data matrix is found automatically: the largest
%            2-D numeric array, or the first cell array of vectors. Rows are
%            time series unless the matrix is obviously transposed (more
%            columns than rows is assumed to be series-by-time).
%
%---OPTIONS (name/value)
% 'DataVar'        name of the data variable, if auto-detection picks wrong
% 'Groups'         class labels: numeric vector, cell of strings, or the name of
%                  a variable in the file. Required for separability analysis.
% 'Labels'         per-series names (default: Series_1, Series_2, ...)
% 'OutputDir'      where to write everything (default: hctsa2_run_<timestamp>)
% 'FeatureSet'     'hctsa' (all ~7700) or 'catch22' (fast, 22). Default 'hctsa'
% 'NumTopFeatures' how many features to plot and rank in detail (default 40)
% 'MaxFeatures'    cap on feature columns in SeriesTable (default Inf = all)
% 'TestStatistic'  '' (default) lets TS_TopFeatures choose: 'ustat' for two
%                  classes, 'classification' for more. Setting 'ustat'
%                  explicitly ERRORS on multi-class data.
% 'DoPlots'        true (default) or false
% 'Parallel'       use the Parallel Computing Toolbox for extraction (default false)
% 'RestrictFeatures' names (string array or cellstr) of the only features to keep,
%                  e.g. the safe subset from HCTSA2_SamplingSafeFeatures. Applied
%                  after extraction, before ranking, so the separability
%                  statistics are computed on the restricted set only.
% 'ExistingHCTSA'  path to an already-computed HCTSA.mat, to skip extraction
%
%---OUTPUT
% results.FeatureTable, results.SeriesTable, results.classificationAccuracy,
% results.files, results.skipped (from HCTSA2_ReportSkipped)

%-------------------------------------------------------------------------------
%% Parse
%-------------------------------------------------------------------------------
p = inputParser;
addParameter(p,'DataVar','',@ischar);
addParameter(p,'Groups',[]);
addParameter(p,'Labels',{});
addParameter(p,'OutputDir','',@ischar);
addParameter(p,'FeatureSet','hctsa',@ischar);
addParameter(p,'NumTopFeatures',40,@isnumeric);
addParameter(p,'MaxFeatures',Inf,@isnumeric);
addParameter(p,'TestStatistic','',@ischar);   % '' = let TS_TopFeatures choose
addParameter(p,'DoPlots',true,@islogical);
addParameter(p,'Parallel',false,@islogical);
addParameter(p,'RestrictFeatures',[]);   % names from HCTSA2_SamplingSafeFeatures
addParameter(p,'ExistingHCTSA','',@ischar);
% Long-format table columns: name or index; default to positional 1..4
addParameter(p,'IdVar',[]);
addParameter(p,'TimeVar',[]);
addParameter(p,'ValueVar',[]);
addParameter(p,'ClassVar',[]);
parse(p,varargin{:});
opt = p.Results;

if isempty(opt.OutputDir)
    opt.OutputDir = ['hctsa2_run_' datestr(now,'yyyymmdd_HHMMSS')];
end
if exist(opt.OutputDir,'dir') ~= 7, mkdir(opt.OutputDir); end

results = struct();
results.options = opt;
results.files = struct();
t0 = tic;

% Remember which figures were already open, so only figures this run creates
% get saved into the output folder.
preExistingFigs = findobj('Type','figure');

fprintf('\n================ hctsa 2.0 analysis ================\n');
fprintf('output: %s\n\n',opt.OutputDir);

%-------------------------------------------------------------------------------
%% 0. Installation check -- so missing binaries are known before, not after
%-------------------------------------------------------------------------------
try
    chk = HCTSA2_CheckInstall(false);
    fprintf('Installation: %u of %u features available (%.0f%%)\n', ...
        chk.featuresAvailable,chk.featuresTotal, ...
        100*chk.featuresAvailable/chk.featuresTotal);
    if chk.featuresAvailable < chk.featuresTotal
        fprintf('  (run HCTSA2_CheckInstall for detail; the run will continue regardless)\n');
    end
    results.install = chk;
catch
    fprintf('Installation check unavailable.\n');
end

%-------------------------------------------------------------------------------
%% 1. Ingest
%-------------------------------------------------------------------------------
if isempty(opt.ExistingHCTSA)
    fprintf('\n-1- Reading %s\n',dataFile);
    [tsData,labels,groups] = HCTSA2_Ingest(dataFile,opt);
    nSeries = numel(tsData);
    fprintf('    %u time series, lengths %u..%u\n',nSeries, ...
        min(cellfun(@numel,tsData)),max(cellfun(@numel,tsData)));

    [groupLabels,groupNames] = localGroups(groups,nSeries);
    results.groupNames = groupNames;
    fprintf('    %u classes: %s\n',numel(groupNames),strjoin(groupNames,', '));

    % hctsa takes groups through the keyword field
    keywords = cell(nSeries,1);
    for i = 1:nSeries
        keywords{i} = groupNames{groupLabels(i)};
    end

    inpFile = fullfile(opt.OutputDir,'INP_ts.mat');
    timeSeriesData = tsData(:)'; 
    labels = labels(:)'; 
    keywords = keywords(:)'; 
    save(inpFile,'timeSeriesData','labels','keywords');
    results.files.inputFile = inpFile;

    %---------------------------------------------------------------------------
    %% 2. Extract
    %---------------------------------------------------------------------------
    hctsaFile = fullfile(opt.OutputDir,'HCTSA.mat');
    % Scale guard: at tens of thousands of series the full library is a
    % multi-day job, so say so before starting rather than after.
    nFeatEst = 7700; if strcmpi(opt.FeatureSet,'catch22'), nFeatEst = 22; end
    secPerSeries = nFeatEst/7700 * 20;          % ~20 s/series for the full set
    estHours = nSeries*secPerSeries/3600;
    fprintf('\n-2- Extracting features (%s set: ~%u features x %u series)\n', ...
        opt.FeatureSet,nFeatEst,nSeries);
    if estHours > 2
        fprintf(2,'    Rough estimate: %.0f hours single-threaded.\n',estHours);
        fprintf(2,'    Consider ''Parallel'',true, or ''FeatureSet'',''catch22'' for a first pass.\n');
        if ~opt.Parallel
            fprintf(2,'    Proceeding serially in 10 s -- Ctrl-C to abort.\n');
            pause(10);
        end
    end
    TS_Init(inpFile,opt.FeatureSet,[false,false,false],hctsaFile);
    TS_Compute(opt.Parallel,[],[],'missing',hctsaFile,false);
    results.files.hctsaFile = hctsaFile;

    fprintf('\n-2b- What did not compute:\n');
    try
        results.skipped = HCTSA2_ReportSkipped(hctsaFile,true);
    catch err
        fprintf('    (report unavailable: %s)\n',err.message);
    end
else
    hctsaFile = opt.ExistingHCTSA;
    results.files.hctsaFile = hctsaFile;
    fprintf('\n-1/2- Using existing %s\n',hctsaFile);
    groupNames = {};
end

%-------------------------------------------------------------------------------
%% 3. Label and normalise
%-------------------------------------------------------------------------------
fprintf('\n-3- Labelling groups and normalising\n');
try
    TS_LabelGroups(hctsaFile,{},true,false);
catch err
    fprintf('    TS_LabelGroups: %s\n',err.message);
end
% Restrict to a caller-supplied feature subset (irregular sampling, etc.)
if ~isempty(opt.RestrictFeatures)
    try
        Ops = TS_GetFromData(hctsaFile,'Operations');
        nBefore = height(Ops);
        % TS_FilterData works on operation IDs, not names, so map first
        wanted = string(opt.RestrictFeatures);
        [tf,~] = ismember(string(Ops.Name),wanted);
        keepIDs = Ops.ID(tf);
        if isempty(keepIDs)
            fprintf(2,'    RestrictFeatures matched no operations; keeping all.\n');
        else
            TS_FilterData(hctsaFile,[],keepIDs,hctsaFile);
            fprintf('    restricted to %u of %u features\n',numel(keepIDs),nBefore);
            if numel(keepIDs) < numel(wanted)
                fprintf('    (%u requested names had no matching operation)\n', ...
                    numel(wanted)-numel(keepIDs));
            end
        end
    catch err
        fprintf(2,'    RestrictFeatures failed (%s); continuing with all features.\n', ...
            err.message);
    end
end

normFile = '';
try
    normFile = TS_Normalize('mixedSigmoid',[0.7,1.0],hctsaFile,true);
    results.files.normFile = normFile;
catch err
    fprintf('    TS_Normalize failed (%s); continuing with raw values.\n',err.message);
    normFile = hctsaFile;
end

%-------------------------------------------------------------------------------
%% 4. Importance and class separability
%-------------------------------------------------------------------------------
if isempty(opt.TestStatistic)
    fprintf('\n-4- Ranking features by class separability (statistic chosen automatically)\n');
else
    fprintf('\n-4- Ranking features by class separability (%s)\n',opt.TestStatistic);
end
ifeat = []; testStat = [];
whatPlots = {};
if opt.DoPlots, whatPlots = {'histogram','distributions'}; end
try
    % Pass an empty struct: TS_TopFeatures only needs cfnParams for
    % 'classification', and builds its own default when the struct has no
    % fields. Constructing one here would just add a failure mode.
    [ifeat,testStat] = TS_TopFeatures(normFile,opt.TestStatistic,struct(), ...
        'whatPlots',whatPlots,'numTopFeatures',opt.NumTopFeatures,'numFeaturesDistr',16);
catch err
    fprintf('    TS_TopFeatures failed: %s\n',err.message);
    fprintf('    Feature ranking will be omitted from FeatureTable.\n');
end

%-------------------------------------------------------------------------------
%% 5. Classification accuracy
%-------------------------------------------------------------------------------
fprintf('\n-5- Classification accuracy\n');
try
    % TS_Classify builds default cfnParams itself when given []
    results.classificationAccuracy = TS_Classify(normFile,[],0,'doPlot',opt.DoPlots);
    fprintf('    mean accuracy: %.1f%%\n',results.classificationAccuracy);
catch err
    fprintf('    TS_Classify failed: %s\n',err.message);
    results.classificationAccuracy = NaN;
end

%-------------------------------------------------------------------------------
%% 6. Build the two tables
%-------------------------------------------------------------------------------
fprintf('\n-6- Building output tables\n');
[results.FeatureTable,results.SeriesTable] = localTables(normFile,ifeat,testStat,opt);
fprintf('    FeatureTable: %u rows x %u cols\n', ...
    height(results.FeatureTable),width(results.FeatureTable));
fprintf('    SeriesTable : %u rows x %u cols\n', ...
    height(results.SeriesTable),width(results.SeriesTable));

outMat = fullfile(opt.OutputDir,'HCTSA2_Results.mat');
FeatureTable = results.FeatureTable; %#ok<NASGU>
SeriesTable  = results.SeriesTable;  %#ok<NASGU>
save(outMat,'FeatureTable','SeriesTable','-v7.3');
results.files.resultsFile = outMat;

try
    writetable(results.FeatureTable,fullfile(opt.OutputDir,'FeatureTable.csv'));
    if width(results.SeriesTable) < 2000
        writetable(results.SeriesTable,fullfile(opt.OutputDir,'SeriesTable.csv'));
    else
        fprintf('    (SeriesTable too wide for CSV; in the .mat only)\n');
    end
catch
end

if opt.DoPlots
    figs = setdiff(findobj('Type','figure'),preExistingFigs);
    nSaved = 0;
    for i = 1:numel(figs)
        try
            nm = get(figs(i),'Name');
            if isempty(nm)
                nm = sprintf('figure_%u',i);
            else
                nm = matlab.lang.makeValidName(nm);
            end
            saveas(figs(i),fullfile(opt.OutputDir,[nm '.png']));
            nSaved = nSaved + 1;
        catch
        end
    end
    if nSaved > 0
        fprintf('    saved %u figure(s)\n',nSaved);
    end
end

fprintf('\n=================================================\n');
fprintf('Done in %s. Everything is in %s\n',BF_TheTime(toc(t0)),opt.OutputDir);
fprintf('=================================================\n\n');
end

%===============================================================================
function [groupLabels,groupNames] = localGroups(groups,n)
if iscell(groups) || isa(groups,'string')
    groups = cellstr(groups);
    groupNames = unique(groups,'stable');
    groupLabels = zeros(n,1);
    for i = 1:numel(groupNames)
        groupLabels(strcmp(groups,groupNames{i})) = i;
    end
else
    groups = double(groups(:));
    u = unique(groups,'stable');
    groupNames = arrayfun(@(g) sprintf('class%g',g),u,'UniformOutput',false);
    groupLabels = zeros(n,1);
    for i = 1:numel(u)
        groupLabels(groups == u(i)) = i;
    end
end
groupNames = groupNames(:)';
end

%===============================================================================
function [FeatureTable,SeriesTable] = localTables(normFile,ifeat,testStat,opt)
D = load(normFile);
TS_DataMat = D.TS_DataMat;
Operations = D.Operations;
TimeSeries = D.TimeSeries;
[nTS,nOps] = size(TS_DataMat);

% group index per series
grp = repmat("unknown",nTS,1);
try
    if any(strcmp(TimeSeries.Properties.VariableNames,'Group'))
        g = TimeSeries.Group;
        grp = string(g);   % string() converts categorical, cellstr and char alike
    else
        grp = string(TimeSeries.Keywords);
    end
catch
end

%--- FeatureTable: one row per feature
FeatureTable = table();
FeatureTable.Name = string(Operations.Name);
try FeatureTable.Keywords = string(Operations.Keywords); catch, end
try FeatureTable.ID = Operations.ID; catch, end

% separability statistic and rank, if TS_TopFeatures succeeded
stat = nan(nOps,1); rank = nan(nOps,1);
if ~isempty(ifeat) && ~isempty(testStat)
    stat(ifeat) = testStat(ifeat);
    rank(ifeat) = (1:numel(ifeat))';
end
FeatureTable.Separability = stat;
FeatureTable.Rank = rank;

% fraction of series on which the feature actually computed
FeatureTable.PctComputed = 100*sum(isfinite(TS_DataMat),1)'/nTS;

% per-group mean
ug = unique(grp,'stable');
for i = 1:numel(ug)
    sel = (grp == ug(i));
    nm = matlab.lang.makeValidName("mean_"+ug(i));
    FeatureTable.(nm) = mean(TS_DataMat(sel,:),1,'omitnan')';
end
FeatureTable = sortrows(FeatureTable,'Rank','ascend','MissingPlacement','last');

%--- SeriesTable: one row per time series
SeriesTable = table();
try SeriesTable.Name = string(TimeSeries.Name); catch, SeriesTable.Name = "series"+(1:nTS)'; end
SeriesTable.Group = grp;
try SeriesTable.Length = TimeSeries.Length; catch, end

keep = 1:nOps;
if ~isfinite(opt.MaxFeatures) && nTS*nOps > 5e7
    fprintf(2,'    SeriesTable will be %u x %u (~%.1f GB). Set ''MaxFeatures'' to trim.\n', ...
        nTS,nOps,8*nTS*nOps/1e9);
end
if isfinite(opt.MaxFeatures) && opt.MaxFeatures < nOps
    if ~isempty(ifeat)
        keep = ifeat(1:min(opt.MaxFeatures,numel(ifeat)));
    else
        keep = 1:opt.MaxFeatures;
    end
    fprintf('    (SeriesTable restricted to %u features)\n',numel(keep));
end
featNames = matlab.lang.makeValidName(string(Operations.Name(keep)));
% Make the feature names unique AGAINST the metadata columns already present,
% not just among themselves. Without the exclusion list a feature named "Name",
% "Group" or "Length" would produce a duplicate variable name and the
% horizontal concatenation below would error. No such feature exists today --
% checked across all 7,693 -- but the feature set is upstream's to change, and
% the exclusion list costs nothing.
featNames = matlab.lang.makeUniqueStrings(featNames,SeriesTable.Properties.VariableNames);
SeriesTable = [SeriesTable, array2table(TS_DataMat(:,keep),'VariableNames',cellstr(featNames))];
end
