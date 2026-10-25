function results = HCTSA2_Analyse(data,varargin)
% HCTSA2_Analyse   End-to-end hctsa2 pipeline: data in, three tables and plots out.
%
%   results = HCTSA2_Analyse('mydata.mat')
%   results = HCTSA2_Analyse('mydata.mat','Groups',g,'OutputDir','run1')
%   results = HCTSA2_Analyse(X,'Groups',g)              % matrix, cell array or table
%   results = HCTSA2_Analyse([],'ExistingHCTSA','run1\HCTSA.mat','Groups',g2)
%
% Runs the whole hctsa pipeline unattended and returns three consolidated tables:
%
%   results.FeatureTable  one row per FEATURE      -- whether it was kept, how
%                         well it separates the classes, its rank, per-group
%                         statistics
%   results.SeriesTable   one row per TIME SERIES  -- its group and its feature
%                         values, ready for downstream modelling
%   results.FeatureAUC    one row per kept FEATURE -- how well it separates every
%                         pair of classes and each class from the rest (AUC), with
%                         q-values; written to FeatureAUC.csv
%
% hctsa already provides the analysis (TS_TopFeatures, TS_Classify and friends).
% This wrapper orchestrates them, keeps going when individual steps fail, says
% what it dropped, and builds the tables -- which hctsa does not produce,
% scattering results across HCTSA.mat, figure windows and function return values.
%
% It never asks a question and never modifies its input: everything is written
% to OutputDir (an ExistingHCTSA file is copied there first). Re-running into an
% OutputDir that holds the same data and feature selection resumes the
% extraction where it stopped; progress is saved every 10 minutes.
%
%---INPUT
% data       a .mat file, a .csv file, or the data itself (a numeric matrix, a
%            cell array of vectors, or a table). In a .mat the data are found
%            automatically: a table, the first cell array of vectors, or the
%            largest 2-D numeric array. Rows are time series unless the matrix is
%            obviously transposed (more rows than columns is read as time-by-series).
%            Long-format data (id, time, value[, class]; one row per observation)
%            is detected. NaN/Inf at either end of a series is padding and is
%            trimmed; a series with NaN/Inf inside it is dropped and listed.
%            Use [] with 'ExistingHCTSA'.
%
%---OPTIONS (name/value)
% Data
% 'Groups'         class labels: numeric, logical, categorical or text, one per
%                  series, or (for a .mat) the name of a variable in the file.
%                  Missing labels (NaN, <undefined>, '') leave a series unlabelled.
%                  Required for separability analysis. With 'ExistingHCTSA', relabels.
% 'Labels'         per-series names (default: from the file, else Series_1, ...)
% 'DataVar'        name of the data variable, if auto-detection picks wrong
% 'Format'         'auto' (default), 'long' or 'wide': how to read a table
% 'IdVar','TimeVar','ValueVar','ClassVar'  long-format columns (name or index)
% 'ExistingHCTSA'  an already-computed HCTSA.mat, to skip extraction (copied)
%
% Features and extraction
% 'FeatureSet'     'hctsa' (all 6,392 in 3.0.0), 'catch22' or 'catch24'. Default 'hctsa'
% 'Features'       names of the only features to compute (wildcards allowed,
%                  e.g. 'SB_*'); applied BEFORE extraction. 'RestrictFeatures' is
%                  the same option under its 2.2.0 name
% 'ExcludeFeatures' names of features not to compute, e.g. 'NL_DVV*' (the most
%                  expensive operation in the library)
% 'Parallel'       use the Parallel Computing Toolbox for extraction (default false)
% 'ParallelMode'   'series' (default): each worker computes whole time series --
%                      near-linear speed-up for many series;
%                  'operations': the operations of one series are spread over the
%                      workers -- for a few very long series.
% 'Resume'         true (default): an OutputDir holding the same data resumes;
%                  false: refuse an OutputDir that already holds a run
%
% Normalisation
% 'Normalization'  'mixedSigmoid' (default), 'scaledRobustSigmoid', 'robustSigmoid',
%                  'scaledSigmoid', 'sigmoid', 'zscore', 'maxmin' or 'none'
% 'FilterThresholds' [series feature]: the fraction of good values a time series
%                  and a feature need to be kept. Default [0.7 0.7]: features that
%                  fail on up to 30% of series are kept, so they can show whether
%                  they separate the classes ([0.7 1] is hctsa's default)
% 'ClassVarFilter' also drop features that are constant within any class (false)
% 'NaNToZero'      after normalisation, set the remaining NaNs to 0 in HCTSA_N.mat
%                  (false). 0 is the BOTTOM of the sigmoid scales (0.5 is the
%                  median), so a failed feature reads as "lowest value".
%
% Ranking
% 'TestStatistic'  '' (default): 'ustatP' (-log10 p of a two-sided Mann-Whitney
%                  test) for two classes; 'meanAUC' for more: the mean of the AUCs
%                  of every pair of classes (Hand & Till 2001), ties broken by the
%                  lowest pair's AUC, then by ID. Also 'meanAUC' for two classes,
%                  'classification' (2.3.0's default for more) or any statistic
%                  TS_TopFeatures accepts. (3.0.0: 'ustat' is upstream's Mann-Whitney
%                  U minus its chance value n1*n2/2, ranked by magnitude; 2.2.1-2.5.x
%                  used the signed z statistic. 'ustatP' is unchanged.)
% 'NumTopFeatures' how many features to plot and list (default 40)
% 'RedundancyThreshold' among the top features, mark each one whose |Spearman
%                  correlation| with a better-ranked one is at least this (0.8);
%                  0 turns it off. Also the cut for 'ReduceFeatureSet'
% 'RedundancyTopN' how many of the top-ranked features that check covers (200)
% 'ReduceFeatureSet' also reduce the WHOLE feature set (false): every feature of
%                  HCTSA_N.mat is clustered by complete linkage on 1 - |Spearman
%                  rho| (hctsa's TS_ReduceFeatureSet), so every pair within a cluster
%                  has |rho| >= RedundancyThreshold (which must then be between 0
%                  and 1). Each cluster keeps one feature: its best-ranked member,
%                  or its most central one when it has no ranked member or the run
%                  has no ranking (fewer than 2 classes). Writes HCTSA_N_reduced.mat:
%                  HCTSA_N.mat with only the kept features, plus a reducedFeatureSet
%                  table (every feature's cluster and the feature kept for it). Adds
%                  Cluster, ClusterSize and Representative (the feature kept for the
%                  cluster; "" for that feature itself) to FeatureTable and
%                  FeatureAUC. Two kept features can still have |rho| above the
%                  threshold: complete linkage can split a correlated pair between
%                  clusters. Correlations use the series both features have (with
%                  'NaNToZero', the zeros count as values). Independent of the top-N
%                  check above
%
% Classification (a rough indication of how well the features separate the
% classes, not a tuned model)
% 'Classifier'     default 'tree' (quick, multi-class); any classifier TS_Classify
%                  accepts: 'svm-linear-lowdim' (hctsa's default), 'svm-linear',
%                  'svm-rbf', 'knn', 'linear', 'diaglinear', 'logistic' (2 classes)...
%                  A tree splits a node only if it holds at least 10 series, and
%                  with 3+ classes hctsa trains one tree per pair of classes, so
%                  with a handful of series per class the estimate sits at chance
%                  (HCTSA2_TestRun's 12 series: 30%, against 90% from
%                  'svm-linear-lowdim'). Unlike the tree, most classifiers drop
%                  every series that has a NaN, which FilterThresholds < 1 allows.
% 'NumFolds'       cross-validation folds (default: 2, 5 or 10 by the size of the
%                  smallest class, hctsa's HowManyFolds); 0 = in-sample
% 'NumRepeats'     repeats of the cross-validation (default 1)
% 'NumNulls'       label-shuffled nulls for a permutation p-value (default 0)
%
% Output
% 'OutputDir'      where to write everything (default: hctsa2_run_<timestamp>)
% 'MaxFeatures'    cap on feature columns in SeriesTable (default Inf = all)
% 'DoPlots'        true (default) or false
% 'ShowFigures'    false draws figures invisibly (they are still saved)
% 'CloseFigures'   close the figures this run made after saving them (true)
% 'FigureFormat'   'png' (default), 'fig' or 'both'
%
%---OUTPUT
% results.FeatureTable, results.SeriesTable, results.FeatureAUC,
% results.DroppedSeries, results.classification (accuracy, chance level, folds,
% permutation p-value), results.classificationAccuracy, results.files,
% results.skipped (from HCTSA2_ReportSkipped), results.runInfo (options, versions,
% step status). With 'ReduceFeatureSet', results.reducedFeatureSet (the table saved
% in HCTSA_N_reduced.mat) and results.files.reducedFile.
% A log of the console output is written to OutputDir/HCTSA2_Analyse_log.txt.
%
% FeatureAUC (2+ classes) is computed on the raw values of the labelled series
% that normalisation kept, and sorted by Rank. An AUC is the chance that a series
% of one class has a higher value than a series of the other (ties count half),
% folded so that 0.5 is no separation and 1 is perfect. Folded, it sits above 0.5
% by chance: with no real difference it averages about 0.68 at 4 series per class,
% 0.61 at 10 and 0.57 at 20 (two classes) -- the q-values say which are real.
%   MeanAUC, MinAUC      the mean and the lowest of the pairwise AUCs; NaN if the
%                        feature computed on no series of some class
%   KW_p, KW_q           Kruskal-Wallis test across all classes; its q-value
%                        (Benjamini-Hochberg over the features)
%   AUC_A_vs_B, Higher_A_vs_B, q_A_vs_B   per pair of classes: the AUC, the class
%                        with the higher values, the Mann-Whitney q-value (BH over
%                        the features, pair by pair)
%   AUC_A_vs_rest, Higher_A_vs_rest       (3+ classes) one class against the rest
%   PctComputed, N_A     % of those series the feature computed on; per class, on
%                        how many (the class size, if on all)
%   RedundantWith, RedundancyRho, Distinct   from FeatureTable, when the check ran:
%                        Distinct is "yes", "no" (it duplicates RedundantWith) or
%                        "not checked"
%   Cluster, ClusterSize, Representative     from FeatureTable, with 'ReduceFeatureSet'
% The p-values are approximations (Mann-Whitney: normal, with tie and continuity
% corrections, as ranksum 'approx'; Kruskal-Wallis: chi-square, tie-corrected),
% rough for classes of only a few series.

%-------------------------------------------------------------------------------
%% Parse and validate everything before any work
%-------------------------------------------------------------------------------
isText = @(s) ischar(s) || (isa(s,'string') && isscalar(s));
isFlag = @(x) (islogical(x) || isnumeric(x)) && isscalar(x);
isCount = @(x) isnumeric(x) && isscalar(x) && x >= 0;

p = inputParser;
addParameter(p,'DataVar','',isText);
addParameter(p,'Groups',[]);
addParameter(p,'Labels',{});
addParameter(p,'Format','auto',isText);
addParameter(p,'IdVar',[]);
addParameter(p,'TimeVar',[]);
addParameter(p,'ValueVar',[]);
addParameter(p,'ClassVar',[]);
addParameter(p,'ExistingHCTSA','',isText);
addParameter(p,'FeatureSet','hctsa',isText);
addParameter(p,'Features',[]);
addParameter(p,'RestrictFeatures',[]);
addParameter(p,'ExcludeFeatures',[]);
addParameter(p,'Parallel',false,isFlag);
addParameter(p,'ParallelMode','series',isText);
addParameter(p,'Resume',true,isFlag);
addParameter(p,'Normalization','mixedSigmoid',isText);
addParameter(p,'FilterThresholds',[0.7 0.7],@(x) isnumeric(x) && numel(x) == 2 && all(x >= 0 & x <= 1));
addParameter(p,'ClassVarFilter',false,isFlag);
addParameter(p,'NaNToZero',false,isFlag);
addParameter(p,'TestStatistic','',isText);
addParameter(p,'NumTopFeatures',40,@(x) isCount(x) && x >= 1);
addParameter(p,'RedundancyThreshold',0.8,@(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
addParameter(p,'RedundancyTopN',200,@(x) isCount(x) && x >= 1);
addParameter(p,'ReduceFeatureSet',false,isFlag);
addParameter(p,'Classifier','tree',isText);
addParameter(p,'NumFolds',[],@(x) isempty(x) || isCount(x));
addParameter(p,'NumRepeats',1,@(x) isCount(x) && x >= 1);
addParameter(p,'NumNulls',0,isCount);
addParameter(p,'OutputDir','',isText);
addParameter(p,'MaxFeatures',Inf,@(x) isCount(x) && x >= 1);
addParameter(p,'DoPlots',true,isFlag);
addParameter(p,'ShowFigures',true,isFlag);
addParameter(p,'CloseFigures',true,isFlag);
addParameter(p,'FigureFormat','png',isText);
parse(p,varargin{:});
opt = p.Results;

for f = {'DataVar','Format','ExistingHCTSA','FeatureSet','ParallelMode','Normalization', ...
         'TestStatistic','Classifier','OutputDir','FigureFormat'}
    opt.(f{1}) = char(opt.(f{1}));
end
for f = {'Parallel','Resume','ClassVarFilter','NaNToZero','ReduceFeatureSet','DoPlots', ...
         'ShowFigures','CloseFigures'}
    opt.(f{1}) = logical(opt.(f{1}));
end
if opt.ReduceFeatureSet && (opt.RedundancyThreshold == 0 || opt.RedundancyThreshold == 1)
    error('HCTSA2_Analyse:reduce', ['ReduceFeatureSet needs a RedundancyThreshold between 0 and 1 ' ...
        '(0 would put every feature in one cluster, 1 would cut the dendrogram at height 0).']);
end
opt.FilterThresholds = double(opt.FilterThresholds(:)');
opt.Format = localPick(opt.Format,{'auto','long','wide'},'Format');
opt.FeatureSet = localPick(opt.FeatureSet,{'hctsa','catch22','catch24'},'FeatureSet');
opt.ParallelMode = localPick(opt.ParallelMode,{'series','operations'},'ParallelMode');
opt.Normalization = localPick(opt.Normalization,{'mixedSigmoid','scaledRobustSigmoid', ...
    'robustSigmoid','scaledSigmoid','sigmoid','zscore','maxmin','none'},'Normalization');
opt.FigureFormat = localPick(opt.FigureFormat,{'png','fig','both'},'FigureFormat');
opt.Classifier = localPick(opt.Classifier,{'tree','knn','linear','linclass','diaglinear', ...
    'svm','svm-linear','svm-linear-lowdim','svm-linear-highdim','svm-rbf','logistic', ...
    'fast-linear'},'Classifier');
twoClassStats = {'ustat','ranksum','ustatP','ranksumP','ustatPsign','ranksumPsign', ...
    'ustatExactP','ranksumExactP','ttest','tstat','ttestP','tstatP','svmBeta'};
if ~isempty(opt.TestStatistic)
    opt.TestStatistic = localPick(opt.TestStatistic,[twoClassStats, ...
        {'meanAUC','classification','medianClassify','nearestMedian'}],'TestStatistic');
end
if ~isempty(opt.RestrictFeatures)
    if ~isempty(opt.Features)
        error('HCTSA2_Analyse:features','Give Features or RestrictFeatures, not both.');
    end
    opt.Features = opt.RestrictFeatures;
end
opt.Features = localNameList(opt.Features,'Features');
opt.ExcludeFeatures = localNameList(opt.ExcludeFeatures,'ExcludeFeatures');
useExisting = ~isempty(opt.ExistingHCTSA);
if useExisting && ~isfile(opt.ExistingHCTSA)
    error('HCTSA2_Analyse:noFile','ExistingHCTSA file %s not found.',opt.ExistingHCTSA);
end
if ~useExisting && isempty(data)
    error('HCTSA2_Analyse:noData','Give data, or an ExistingHCTSA file.');
end

if isempty(opt.OutputDir)
    opt.OutputDir = ['hctsa2_run_' datestr(now,'yyyymmdd_HHMMSS')];
end
if exist(opt.OutputDir,'dir') ~= 7, mkdir(opt.OutputDir); end

%-------------------------------------------------------------------------------
%% Session state this function changes, and restores on exit
%-------------------------------------------------------------------------------
logFile = fullfile(opt.OutputDir,'HCTSA2_Analyse_log.txt');
prevDiary = get(0,'Diary'); prevDiaryFile = get(0,'DiaryFile');
diary(logFile);
restoreDiary = onCleanup(@() localRestoreDiary(prevDiary,prevDiaryFile)); %#ok<NASGU>
rngOnEntry = rng;
restoreRng = onCleanup(@() rng(rngOnEntry)); %#ok<NASGU>
if ~opt.ShowFigures
    prevVisible = get(groot,'DefaultFigureVisible');
    set(groot,'DefaultFigureVisible','off');
    restoreVisible = onCleanup(@() set(groot,'DefaultFigureVisible',prevVisible)); %#ok<NASGU>
end
% Remember which figures were already open, so only figures this run creates
% are saved (and closed).
preExistingFigs = findobj('Type','figure');

results = struct();
results.options = opt;
results.files = struct('log',logFile);
results.classificationAccuracy = NaN;
results.classification = struct();
steps = struct('step',{},'ok',{},'note',{});
t0 = tic;

fprintf('\n================ hctsa2 analysis ================\n');
fprintf('%s\noutput: %s\n\n',datestr(now),opt.OutputDir);

%-------------------------------------------------------------------------------
%% 0. Installation check -- so missing binaries are known before, not after
%-------------------------------------------------------------------------------
chk = [];
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

hctsaFile = fullfile(opt.OutputDir,'HCTSA.mat');
droppedIngest = struct('Name',{{}},'Reason',{{}});

if ~useExisting
    %---------------------------------------------------------------------------
    %% 1. Ingest
    %---------------------------------------------------------------------------
    if ischar(data) || isa(data,'string')
        fprintf('\n-1- Reading %s\n',char(data));
    else
        fprintf('\n-1- Reading the %s given\n',class(data));
    end
    ing = struct('DataVar',opt.DataVar,'Labels',{opt.Labels},'Groups',{opt.Groups}, ...
        'IdVar',{opt.IdVar},'TimeVar',{opt.TimeVar},'ValueVar',{opt.ValueVar}, ...
        'ClassVar',{opt.ClassVar},'Format',opt.Format);
    [tsData,labels,groups,droppedIngest] = HCTSA2_Ingest(data,ing);
    nSeries = numel(tsData);
    fprintf('    %u time series, lengths %u..%u\n',nSeries, ...
        min(cellfun(@numel,tsData)),max(cellfun(@numel,tsData)));
    [groupIdx,groupNames] = localGroups(groups,nSeries);
    localReportGroups(groupIdx,groupNames);
    localCheckAgainstClasses(opt,groupIdx,groupNames,twoClassStats);
    results.groupNames = groupNames;

    % Keywords record the group (hctsa splits keywords on commas, so none here)
    keywords = repmat({''},nSeries,1);
    keywords(groupIdx > 0) = strrep(groupNames(groupIdx(groupIdx > 0)),',',';');

    %---------------------------------------------------------------------------
    %% 2. Extract (or resume)
    %---------------------------------------------------------------------------
    inpFile = fullfile(opt.OutputDir,'INP_ts.mat');
    specFile = fullfile(opt.OutputDir,'HCTSA2_RunSpec.mat');
    inp = struct('timeSeriesData',{tsData(:)'},'labels',{labels(:)'},'keywords',{keywords(:)'});
    spec = struct('FeatureSet',opt.FeatureSet,'Features',{opt.Features}, ...
        'ExcludeFeatures',{opt.ExcludeFeatures});
    if isfile(hctsaFile)
        if ~opt.Resume
            error('HCTSA2_Analyse:outputDirInUse', ...
                '%s already holds a run (HCTSA.mat). Use a new OutputDir, or Resume.',opt.OutputDir);
        end
        same = isfile(inpFile) && isfile(specFile) && isequaln(load(inpFile),inp);
        if same
            old = load(specFile);
            same = isfield(old,'spec') && isequaln(old.spec,spec);
        end
        if ~same
            error('HCTSA2_Analyse:outputDirInUse', ['%s already holds a run on different ' ...
                'data or a different feature selection. Use a new OutputDir.'],opt.OutputDir);
        end
        fprintf('\n-2- Resuming the extraction in %s\n',hctsaFile);
        localSelectFeatures(hctsaFile,opt);   % in case the last run stopped before it
    else
        save(inpFile,'-struct','inp');
        save(specFile,'spec');
        fprintf('\n-2- Extracting features (%s set)\n',opt.FeatureSet);
        TS_Init(inpFile,opt.FeatureSet,[false,false,false],hctsaFile);
        localSelectFeatures(hctsaFile,opt);
    end
    results.files.inputFile = inpFile;
    results.files.hctsaFile = hctsaFile;
    localCompute(hctsaFile,opt);

    fprintf('\n-2b- What did not compute:\n');
    try
        results.skipped = HCTSA2_ReportSkipped(hctsaFile,true,chk);
    catch err
        fprintf('    (report unavailable: %s)\n',err.message);
    end

    % Group labels, written straight into the run's own HCTSA.mat (no prompt)
    localWriteGroups(hctsaFile,groupIdx,groupNames,[]);
    steps(end+1) = struct('step','extract','ok',true,'note','');
else
    %---------------------------------------------------------------------------
    %% 1/2. An existing HCTSA.mat: copied, never modified
    %---------------------------------------------------------------------------
    src = opt.ExistingHCTSA;
    fprintf('\n-1/2- Using existing %s\n',src);
    if ~strcmp(localFullPath(src),localFullPath(hctsaFile))
        if isfile(hctsaFile)
            error('HCTSA2_Analyse:outputDirInUse','%s already holds an HCTSA.mat.',opt.OutputDir);
        end
        copyfile(src,hctsaFile);
        fprintf('    copied to %s\n',hctsaFile);
    end
    results.files.hctsaFile = hctsaFile;
    localSelectFeatures(hctsaFile,opt);
    TimeSeries = TS_GetFromData(hctsaFile,'TimeSeries');
    nSeries = height(TimeSeries);
    if ~isempty(opt.Labels)
        newNames = cellstr(opt.Labels(:));
        if numel(newNames) ~= nSeries || numel(unique(newNames)) ~= nSeries
            error('HCTSA2_Analyse:labels','Labels needs %u unique names.',nSeries);
        end
        TimeSeries.Name = newNames;
    end
    if ~isempty(opt.Groups)
        if isText(opt.Groups)
            error('HCTSA2_Analyse:groups',['With ExistingHCTSA, Groups must be the labels ' ...
                'themselves (one per series), not a variable name.']);
        end
        [groupIdx,groupNames] = localGroups(opt.Groups,nSeries);
    elseif ismember('Group',TimeSeries.Properties.VariableNames) && ...
            iscategorical(TimeSeries.Group) && ~all(isundefined(TimeSeries.Group))
        [groupIdx,groupNames] = localGroups(TimeSeries.Group,nSeries);
    else
        % one keyword per series is taken as its group; anything else needs Groups
        kw = SUB_cell2cellcell(TimeSeries.Keywords);
        if all(cellfun(@numel,kw) == 1)
            [groupIdx,groupNames] = localGroups(cellfun(@(k) k{1},kw,'UniformOutput',false),nSeries);
        else
            fprintf(2,['    Series carry several keywords each, so the groups are ambiguous: ' ...
                'pass ''Groups''. Continuing unlabelled.\n']);
            groupIdx = zeros(nSeries,1); groupNames = {};
        end
    end
    localReportGroups(groupIdx,groupNames);
    localCheckAgainstClasses(opt,groupIdx,groupNames,twoClassStats);
    results.groupNames = groupNames;
    localWriteGroups(hctsaFile,groupIdx,groupNames,TimeSeries);
    try
        results.skipped = HCTSA2_ReportSkipped(hctsaFile,true,chk);
    catch err
        fprintf('    (report unavailable: %s)\n',err.message);
    end
    steps(end+1) = struct('step','load existing','ok',true,'note','');
end

%-------------------------------------------------------------------------------
%% 3. Normalise
%-------------------------------------------------------------------------------
fprintf('\n-3- Normalising (%s, filter [%g %g]%s)\n',opt.Normalization, ...
    opt.FilterThresholds,localIf(opt.ClassVarFilter,', class-variance filter',''));
normFile = '';
try
    normFile = TS_Normalize(opt.Normalization,opt.FilterThresholds,hctsaFile,opt.ClassVarFilter);
    if opt.NaNToZero
        N = load(normFile,'TS_DataMat','TS_Quality');
        isNaN = isnan(N.TS_DataMat);
        TS_DataMat = N.TS_DataMat; TS_DataMat(isNaN) = 0;
        % TS_LoadData turns cells with TS_Quality > 0 back into NaN, so these are
        % marked good here (their real quality stays in HCTSA.mat)
        TS_Quality = N.TS_Quality; TS_Quality(isNaN) = 0; %#ok<NASGU>
        save(normFile,'TS_DataMat','TS_Quality','-append');
        fprintf('    NaNToZero: %u NaN values (%.2f%%) set to 0 in %s\n', ...
            nnz(isNaN),100*mean(isNaN(:)),normFile);
    end
    results.files.normFile = normFile;
    steps(end+1) = struct('step','normalise','ok',true,'note','');
catch err
    fprintf(2,'    TS_Normalize failed (%s).\n    CONTINUING WITH UN-NORMALISED VALUES.\n',err.message);
    normFile = hctsaFile;
    steps(end+1) = struct('step','normalise','ok',false,'note',err.message);
end

% Classes that survive normalisation (a class can lose every member)
TSn = TS_GetFromData(normFile,'TimeSeries');
if ismember('Group',TSn.Properties.VariableNames)
    classesNow = categories(removecats(TSn.Group));
else
    classesNow = {};
end
lostClasses = setdiff(groupNames,classesNow);
if ~isempty(lostClasses)
    fprintf(2,'    Every series of class(es) %s was dropped by normalisation.\n', ...
        strjoin(lostClasses,', '));
end
numClasses = numel(classesNow);

%-------------------------------------------------------------------------------
%% 4. Rank features by class separability
%-------------------------------------------------------------------------------
ifeat = []; testStat = []; statName = opt.TestStatistic; sep = [];
whatPlots = {};
if opt.DoPlots, whatPlots = {'histogram','distributions'}; end
if numClasses < 2
    fprintf('\n-4- Ranking skipped: needs at least 2 labelled classes (have %u)\n',numClasses);
    steps(end+1) = struct('step','rank','ok',false,'note','fewer than 2 classes');
else
    if isempty(statName)
        statName = localIf(numClasses == 2,'ustatP','meanAUC');
    end
    fprintf('\n-4- Ranking features by class separability (%s)\n',statName);
    try
        sep = localClassStats(hctsaFile,normFile);   % for FeatureAUC, and 'meanAUC'
    catch err
        fprintf(2,'    AUC statistics failed: %s\n',err.message);
        steps(end+1) = struct('step','auc','ok',false,'note',err.message);
    end
    try
        if strcmp(statName,'meanAUC')
            if isempty(sep), error('HCTSA2_Analyse:noAUC','no AUC statistics'); end
            [ifeat,testStat] = localRankByAUC(sep,normFile,opt);
        else
            % An empty struct: TS_TopFeatures builds its own (in-sample) classifier
            % settings for 'classification'.
            [ifeat,testStat] = TS_TopFeatures(normFile,statName,struct(), ...
                'whatPlots',whatPlots,'numTopFeatures',opt.NumTopFeatures, ...
                'numFeaturesDistr',min(16,opt.NumTopFeatures));
        end
        steps(end+1) = struct('step','rank','ok',true,'note','');
    catch err
        fprintf(2,'    Ranking failed: %s\n',err.message);
        steps(end+1) = struct('step','rank','ok',false,'note',err.message);
    end
end
results.separabilityStatistic = statName;

%-------------------------------------------------------------------------------
%% 4b. Reduce the whole feature set (optional, 'ReduceFeatureSet')
%-------------------------------------------------------------------------------
reduced = [];
if opt.ReduceFeatureSet
    fprintf('\n-4b- Reducing the feature set: clusters of |rho| >= %g\n',opt.RedundancyThreshold);
    if strcmp(normFile,hctsaFile)
        fprintf(2,'    Skipped: there is no HCTSA_N.mat (normalisation failed).\n');
        steps(end+1) = struct('step','reduce','ok',false,'note','no HCTSA_N.mat');
    else
        try
            [reduced,results.files.reducedFile] = localReduce(normFile,ifeat,statName,opt);
            results.reducedFeatureSet = reduced;
            steps(end+1) = struct('step','reduce','ok',true,'note','');
        catch err
            reduced = [];
            fprintf(2,'    Reducing the feature set failed: %s\n',err.message);
            steps(end+1) = struct('step','reduce','ok',false,'note',err.message);
        end
    end
end

%-------------------------------------------------------------------------------
%% 5. Classification accuracy
%-------------------------------------------------------------------------------
if numClasses < 2
    fprintf('\n-5- Classification skipped: needs at least 2 labelled classes\n');
    steps(end+1) = struct('step','classify','ok',false,'note','fewer than 2 classes');
else
    fprintf('\n-5- Classification accuracy (%s classifier)\n',opt.Classifier);
    try
        [results.classification,results.classificationAccuracy] = ...
            localClassify(normFile,opt,TSn,classesNow);
        steps(end+1) = struct('step','classify','ok',true,'note','');
    catch err
        fprintf(2,'    TS_Classify failed: %s\n',err.message);
        steps(end+1) = struct('step','classify','ok',false,'note',err.message);
    end
end

if opt.DoPlots && numClasses >= 2
    try
        TS_PlotLowDim(normFile,'pca',true,[],[]);
    catch err
        fprintf('    (PCA plot not drawn: %s)\n',err.message);
    end
end

%-------------------------------------------------------------------------------
%% 6. Build the tables
%-------------------------------------------------------------------------------
fprintf('\n-6- Building output tables\n');
[results.FeatureTable,results.SeriesTable,results.DroppedSeries,results.FeatureAUC] = ...
    localTables(hctsaFile,normFile,ifeat,testStat,opt,droppedIngest,sep,reduced);
fprintf('    FeatureTable: %u rows x %u cols (%u ranked)\n',height(results.FeatureTable), ...
    width(results.FeatureTable),sum(isfinite(results.FeatureTable.Rank)));
fprintf('    SeriesTable : %u rows x %u cols\n', ...
    height(results.SeriesTable),width(results.SeriesTable));
fprintf('    FeatureAUC  : %u rows x %u cols\n', ...
    height(results.FeatureAUC),width(results.FeatureAUC));
if height(results.DroppedSeries) > 0
    fprintf('    DroppedSeries: %u series not in SeriesTable (see results.DroppedSeries)\n', ...
        height(results.DroppedSeries));
end

%-------------------------------------------------------------------------------
%% 7. Save, export, figures
%-------------------------------------------------------------------------------
results.runInfo = localRunInfo(opt,steps,t0);
outMat = fullfile(opt.OutputDir,'HCTSA2_Results.mat');
FeatureTable = results.FeatureTable; %#ok<NASGU>
SeriesTable  = results.SeriesTable;  %#ok<NASGU>
DroppedSeries = results.DroppedSeries; %#ok<NASGU>
FeatureAUC = results.FeatureAUC; %#ok<NASGU>
classification = results.classification; %#ok<NASGU>
runInfo = results.runInfo; %#ok<NASGU>
save(outMat,'FeatureTable','SeriesTable','DroppedSeries','FeatureAUC','classification', ...
    'runInfo','-v7.3');
results.files.resultsFile = outMat;

localWriteCsv(results.FeatureTable,fullfile(opt.OutputDir,'FeatureTable.csv'));
if height(results.FeatureAUC) > 0
    localWriteCsv(results.FeatureAUC,fullfile(opt.OutputDir,'FeatureAUC.csv'));
end
if width(results.SeriesTable) < 2000
    localWriteCsv(results.SeriesTable,fullfile(opt.OutputDir,'SeriesTable.csv'));
else
    fprintf(['    (SeriesTable has %u columns: too wide for a CSV, so it is in the .mat only;\n' ...
             '     HCTSA2_TopN gives a narrower table you can write yourself)\n'],width(results.SeriesTable));
end
if height(results.DroppedSeries) > 0
    localWriteCsv(results.DroppedSeries,fullfile(opt.OutputDir,'DroppedSeries.csv'));
end

if opt.DoPlots
    localSaveFigures(setdiff(findobj('Type','figure'),preExistingFigs),opt);
end

%-------------------------------------------------------------------------------
%% Summary
%-------------------------------------------------------------------------------
results.runInfo = localRunInfo(opt,steps,t0);
failed = steps(~[steps.ok]);
fprintf('\n=================================================\n');
if isempty(failed)
    fprintf('Done in %s. Everything is in %s\n',BF_TheTime(toc(t0)),opt.OutputDir);
else
    fprintf(2,'Finished in %s with %u step(s) incomplete:\n',BF_TheTime(toc(t0)),numel(failed));
    for i = 1:numel(failed)
        fprintf(2,'  %-10s %s\n',failed(i).step,failed(i).note);
    end
    fprintf('Results so far are in %s\n',opt.OutputDir);
end
fprintf('=================================================\n\n');
end

%===============================================================================
function v = localPick(v,allowed,optName)
% Case-insensitive exact match against the allowed values (no abbreviations).
i = find(strcmpi(v,allowed),1);
if isempty(i)
    error('HCTSA2_Analyse:badOption','''%s'' is not a valid %s. Valid: %s.', ...
        v,optName,strjoin(allowed,', '));
end
v = allowed{i};
end

%===============================================================================
function names = localNameList(names,optName)
% Feature names as a row cellstr; [] or {} = none.
if isempty(names), names = {}; return; end
if ischar(names) || isa(names,'string'), names = cellstr(names); end
if ~iscellstr(names) %#ok<ISCLSTR>
    error('HCTSA2_Analyse:badOption','%s must be feature names (text).',optName);
end
names = names(:)';
end

%===============================================================================
function out = localIf(cond,a,b)
if cond, out = a; else, out = b; end
end

%===============================================================================
function localRestoreDiary(prevDiary,prevFile)
diary('off');
if strcmp(prevDiary,'on')
    diary(prevFile);
end
end

%===============================================================================
function f = localFullPath(f)
% Absolute, canonical path of a file that may not exist yet.
[d,n,e] = fileparts(f);
if isempty(d), d = pwd; end
if exist(d,'dir') == 7
    old = cd(d); d = pwd; cd(old);
end
f = fullfile(d,[n e]);
end

%===============================================================================
function [groupIdx,groupNames] = localGroups(groups,n)
% Group index per series (0 = unlabelled) and the group names, in order of
% first appearance. Categorical labels keep their names; NaN, <undefined> and
% '' mean unlabelled.
if isa(groups,'categorical')
    labelsText = cellstr(groups(:));                 % <undefined> -> ''
    labelsText(isundefined(groups(:))) = {''};
elseif iscell(groups) || isa(groups,'string')
    labelsText = cellstr(groups(:));
else
    g = double(groups(:));
    labelsText = arrayfun(@(x) sprintf('class%g',x),g,'UniformOutput',false);
    labelsText(isnan(g)) = {''};
end
if numel(labelsText) ~= n
    error('HCTSA2_Analyse:groupCount','%u group labels for %u time series.',numel(labelsText),n);
end
labelsText = strtrim(labelsText);
isLabelled = ~cellfun(@isempty,labelsText);
groupNames = unique(labelsText(isLabelled),'stable');
[~,idx] = ismember(labelsText(isLabelled),groupNames);
groupIdx = zeros(n,1);
groupIdx(isLabelled) = idx;
groupNames = groupNames(:)';
end

%===============================================================================
function localReportGroups(groupIdx,groupNames)
counts = accumarray(groupIdx(groupIdx > 0),1,[numel(groupNames),1]);
fprintf('    %u classes: %s\n',numel(groupNames), ...
    strjoin(cellfun(@(g,c) sprintf('%s (%u)',g,c),groupNames(:),num2cell(counts), ...
    'UniformOutput',false),', '));
if any(groupIdx == 0)
    fprintf('    %u series unlabelled (kept, but left out of ranking and classification)\n', ...
        sum(groupIdx == 0));
end
if any(counts < 2)
    fprintf(2,'    Class(es) with a single series: %s -- too few to rank or classify on.\n', ...
        strjoin(groupNames(counts < 2),', '));
end
end

%===============================================================================
function localCheckAgainstClasses(opt,groupIdx,groupNames,twoClassStats)
% Settings that cannot work with this number of classes fail now, before extraction.
nc = numel(groupNames);
if ismember(opt.TestStatistic,twoClassStats) && nc ~= 2
    error('HCTSA2_Analyse:testStatistic','''%s'' needs exactly 2 classes; there are %u.', ...
        opt.TestStatistic,nc);
end
if ismember(opt.Classifier,{'logistic','svm-linear-highdim'}) && nc > 2
    error('HCTSA2_Analyse:classifier','''%s'' needs 2 classes; there are %u.',opt.Classifier,nc);
end
if nc >= 2 && strcmp(opt.Classifier,'tree')
    counts = accumarray(groupIdx(groupIdx > 0),1,[nc,1]);
    pairSums = counts + counts';
    pairSums(logical(eye(nc))) = Inf;
    if min(pairSums(:)) < 12
        fprintf(['    Note: some pair of classes has fewer than 12 series between them, so the\n' ...
                 '    default tree classifier will likely sit at chance; consider ' ...
                 '''Classifier'',''svm-linear-lowdim''.\n']);
    end
end
end

%===============================================================================
function localWriteGroups(hctsaFile,groupIdx,groupNames,TimeSeries)
% Write the Group column of TimeSeries into the run's HCTSA.mat (TS_LabelGroups
% does the same from keywords, but asks questions; this cannot).
if isempty(TimeSeries)
    TimeSeries = TS_GetFromData(hctsaFile,'TimeSeries');
end
if isempty(groupNames)
    TimeSeries.Group = categorical(nan(height(TimeSeries),1));
else
    TimeSeries.Group = categorical(groupIdx(:),1:numel(groupNames),groupNames);
end
save(hctsaFile,'TimeSeries','-append');
end

%===============================================================================
function localSelectFeatures(hctsaFile,opt)
% Keep only the requested features (before extraction for a new run).
if isempty(opt.Features) && isempty(opt.ExcludeFeatures), return; end
Ops = TS_GetFromData(hctsaFile,'Operations');
names = Ops.Name;
keep = true(height(Ops),1);
if ~isempty(opt.Features)
    [keep,unmatched] = localMatch(names,opt.Features);
    if ~isempty(unmatched)
        fprintf(2,'    Features: %u name(s) matched nothing: %s\n',numel(unmatched), ...
            strjoin(unmatched(1:min(10,end)),', '));
    end
end
if ~isempty(opt.ExcludeFeatures)
    [drop,unmatched] = localMatch(names,opt.ExcludeFeatures);
    keep = keep & ~drop;
    if ~isempty(unmatched)
        fprintf(2,'    ExcludeFeatures: %u name(s) matched nothing: %s\n',numel(unmatched), ...
            strjoin(unmatched(1:min(10,end)),', '));
    end
end
if ~any(keep)
    error('HCTSA2_Analyse:noFeatures','The feature selection leaves no features.');
end
if all(keep)
    fprintf('    feature selection keeps all %u features\n',height(Ops));
    return
end
TS_Subset(hctsaFile,[],Ops.ID(keep),true,hctsaFile);
fprintf('    feature selection: %u of %u features\n',sum(keep),height(Ops));
end

%===============================================================================
function [hit,unmatched] = localMatch(names,patterns)
% Which names match any pattern (exact, or with * and ? wildcards).
hit = false(numel(names),1);
isUnmatched = false(1,numel(patterns));
for i = 1:numel(patterns)
    if any(patterns{i} == '*' | patterns{i} == '?')
        rx = ['^' regexptranslate('wildcard',patterns{i}) '$'];
        m = ~cellfun(@isempty,regexp(names,rx,'once'));
    else
        m = strcmp(names,patterns{i});
    end
    hit = hit | m(:);
    isUnmatched(i) = ~any(m);
end
unmatched = patterns(isUnmatched);
end

%===============================================================================
function localCompute(hctsaFile,opt)
% TS_Compute with a progress bar and autosave; times the first pending series
% to give a real estimate of the rest.
if opt.Parallel
    parallelHow = opt.ParallelMode;
    fprintf('    in parallel over %s\n',parallelHow);
else
    parallelHow = 'none';
end
autosaveSeconds = 600;
TimeSeries = TS_GetFromData(hctsaFile,'TimeSeries');
Q = TS_GetFromData(hctsaFile,'TS_Quality');
pending = find(any(isnan(Q),2));
if isempty(pending)
    fprintf('    every feature is already computed\n');
    return
end
fprintf('    %u of %u series to compute\n',numel(pending),height(TimeSeries));
if numel(pending) >= 2
    t = tic;
    TS_Compute(parallelHow,TimeSeries.ID(pending(1)),[],'missing',hctsaFile,'fast',false,autosaveSeconds);
    t1 = toc(t);
    workers = 1;
    if strcmp(parallelHow,'series') && exist('gcp','file') == 2
        pool = gcp('nocreate');
        if ~isempty(pool), workers = pool.NumWorkers; end
    end
    eta = t1*(numel(pending)-1)/workers;
    fprintf('    first series: %s; estimate for the remaining %u: %s\n', ...
        BF_TheTime(t1),numel(pending)-1,BF_TheTime(eta));
    if eta > 2*3600
        fprintf(2,['    That is long. Progress is saved every 10 minutes and a rerun resumes.\n' ...
            '    To shorten it: ''Parallel'',true; ''ExcludeFeatures'',''NL_DVV*''; ' ...
            'or ''FeatureSet'',''catch22'' for a first look.\n']);
    end
end
TS_Compute(parallelHow,[],[],'missing',hctsaFile,'fast',false,autosaveSeconds);
end

%===============================================================================
function [C,acc] = localClassify(normFile,opt,TSn,classesNow)
cfnParams = GiveMeDefaultClassificationParams(normFile);
if isempty(fieldnames(cfnParams))
    error('HCTSA2_Analyse:unlabelled','no group labels');
end
cfnParams.whatClassifier = opt.Classifier;
cfnParams = UpdateClassifierText(cfnParams);
if ~isempty(opt.NumFolds), cfnParams.numFolds = opt.NumFolds; end
cfnParams.numRepeats = opt.NumRepeats;

% Most classifiers drop every series with a NaN: say so when there are NaNs
X = TS_GetFromData(normFile,'TS_DataMat');
Q = TS_GetFromData(normFile,'TS_Quality');
X(Q > 0) = NaN;
labelled = ~isundefined(TSn.Group);
hasNaN = any(isnan(X(labelled,:)),2);
if any(hasNaN) && ~strcmp(opt.Classifier,'tree')
    fprintf(2,['    %u of %u labelled series have a NaN feature value; the %s classifier ' ...
        'leaves those series out.\n    Use the tree, or ''NaNToZero'',true.\n'], ...
        sum(hasNaN),sum(labelled),opt.Classifier);
end

fprintf('    %u-fold cross-validation, %u repeat(s), %s; chance level %.1f%%\n', ...
    cfnParams.numFolds,cfnParams.numRepeats,cfnParams.whatLoss,100/numel(classesNow));
if opt.NumNulls > 0
    [acc,nullStats] = TS_Classify(normFile,cfnParams,opt.NumNulls,'doPlot',opt.DoPlots);
else
    acc = TS_Classify(normFile,cfnParams,0,'doPlot',opt.DoPlots); % (no nullStats output then)
    nullStats = [];
end
C = struct('classifier',opt.Classifier,'accuracy',acc,'measure',cfnParams.whatLoss, ...
    'numFolds',cfnParams.numFolds,'numRepeats',cfnParams.numRepeats, ...
    'chance',100/numel(classesNow),'numNulls',opt.NumNulls,'nullMean',NaN,'pValue',NaN);
fprintf('    mean %s: %.1f%% (chance %.1f%%)\n',cfnParams.whatLoss,acc,C.chance);
if opt.NumNulls > 0 && ~isempty(nullStats)
    C.nullMean = mean(nullStats);
    C.pValue = mean(nullStats >= acc);
    fprintf('    permutation test: null mean %.1f%%, p = %.3g (%u nulls)\n', ...
        C.nullMean,C.pValue,opt.NumNulls);
end
end

%===============================================================================
function sep = localClassStats(hctsaFile,normFile)
% localSeparation on the raw values of the labelled series and the features that
% normalisation kept (the data of FeatureTable's two-class AUC column).
R = load(hctsaFile,'TS_DataMat','TS_Quality','Operations','TimeSeries');
N = load(normFile,'Operations','TimeSeries');
[~,col] = ismember(N.Operations.ID,R.Operations.ID);
rows = ismember(R.TimeSeries.ID,N.TimeSeries.ID) & ~isundefined(R.TimeSeries.Group);
grp = removecats(R.TimeSeries.Group(rows));
X = R.TS_DataMat(rows,col);
X(~isfinite(X) | R.TS_Quality(rows,col) > 0) = NaN;
sep = localSeparation(X,double(grp),numel(categories(grp)));
sep.classes = categories(grp)';
sep.ID = N.Operations.ID;
sep.Name = N.Operations.Name;
sep.numSeries = sum(rows);
end

%===============================================================================
function S = localSeparation(X,g,K)
% How well each feature (a column of X; NaN = not computed) separates the K classes
% of g (1..K): every pair's AUC and Mann-Whitney p, each class against the rest,
% and a Kruskal-Wallis p. One sort per feature: the per-class counts in each group
% of tied values give every pair's U statistic (ties count half) and the rank sums.
% p-values: normal (with the continuity correction, as ranksum 'approx') and
% chi-square approximations, both tie-corrected.
nf = size(X,2);
S.pairs = nchoosek(1:K,2);
pa = S.pairs(:,1); pb = S.pairs(:,2);
iab = sub2ind([K K],pa,pb);
nP = numel(pa);
S.n = zeros(K,nf);                                   % series computed, per class
S.auc = nan(nP,nf); S.higher = nan(nP,nf); S.p = nan(nP,nf);
S.aucRest = nan(K,nf); S.higherRest = nan(K,nf);
S.pKW = nan(1,nf);
for j = 1:nf
    ok = ~isnan(X(:,j));
    if ~any(ok), continue; end
    [xs,ix] = sort(X(ok,j));
    gj = g(ok);
    tie = cumsum([true; diff(xs) ~= 0]);             % tie group of each sorted value
    c = accumarray([tie gj(ix)],1,[tie(end) K]);     % counts per tie group and class
    below = [zeros(1,K); cumsum(c(1:end-1,:),1)];    % per class, counts below each group
    U = c'*(below + c/2);                            % U(a,b) = #(a > b) + #(a == b)/2
    n = sum(c,1)'; N = sum(n);
    S.n(:,j) = n;
    % pairs
    nn = n(pa).*n(pb); m = n(pa) + n(pb);
    t = c(:,pa) + c(:,pb);
    v = nn/12.*((m+1) - sum(t.^3 - t,1)'./(m.*(m-1)));   % variance of U, tie-corrected
    w = U(iab) - nn/2;
    p = erfc(abs(w - sign(w)/2)./sqrt(2*v));
    p(v == 0) = 1;                                   % every value tied
    A = U(iab)./nn;
    A(nn == 0) = NaN; p(nn == 0) = NaN;
    S.auc(:,j) = max(A,1-A); S.higher(:,j) = sign(A - 0.5); S.p(:,j) = p;
    % each class against the rest
    A = (sum(U,2) - diag(U))./(n.*(N - n));
    S.aucRest(:,j) = max(A,1-A); S.higherRest(:,j) = sign(A - 0.5);
    % all classes: Kruskal-Wallis
    if all(n > 0)
        tg = sum(c,2);
        R = c'*(cumsum([0; tg(1:end-1)]) + (tg+1)/2); % rank sums (ties: average rank)
        H = 12/(N*(N+1))*sum(R.^2./n) - 3*(N+1);
        C = 1 - sum(tg.^3 - tg)/(N^3 - N);
        if C > 0
            S.pKW(j) = gammainc(max(H,0)/C/2,(K-1)/2,'upper');
        else
            S.pKW(j) = 1;
        end
    end
end
S.meanAUC = mean(S.auc,1);                           % NaN if a pair is not covered
S.minAUC = min(S.auc,[],1);
S.minAUC(any(isnan(S.auc),1)) = NaN;
end

%===============================================================================
function q = localBH(p)
% Benjamini-Hochberg q-values along each row of p (NaN: not tested, left out).
q = nan(size(p));
for i = 1:size(p,1)
    ok = find(~isnan(p(i,:)));
    m = numel(ok);
    [ps,ix] = sort(p(i,ok));
    qs = min(1,ps*m./(1:m));
    q(i,ok(ix)) = fliplr(cummin(fliplr(qs)));
end
end

%===============================================================================
function [ifeat,stat] = localRankByAUC(sep,normFile,opt)
% The features (indices into the normalised file's Operations) by MeanAUC, best
% first, ties by MinAUC; prints the top ones and draws TS_TopFeatures' two figures.
stat = sep.meanAUC(:);
% (rounded: equal means of different pairs' AUCs can differ in the last bit)
[~,ifeat] = sortrows([-round(stat,12),-round(sep.minAUC(:),12),sep.ID(:)]);
ifeat = ifeat(isfinite(stat(ifeat)));
if isempty(ifeat)
    error('HCTSA2_Analyse:noAUC','no feature computed on series of every class');
end
fprintf('    %u of %u features ranked by mean pairwise AUC (0.5 = no separation, 1 = perfect)\n', ...
    numel(ifeat),numel(stat));
for i = 1:min(opt.NumTopFeatures,numel(ifeat))
    j = ifeat(i);
    fprintf('[%u] %s -- MeanAUC %.3f, MinAUC %.3f\n',sep.ID(j),sep.Name{j},stat(j),sep.minAUC(j));
end
if ~opt.DoPlots || isempty(ifeat), return; end
try
    f = figure('color','w'); hold('on')
    f.Position(3:4) = [559,278];
    colors = BF_GetColorMap('spectral',5,1);
    h = histogram(stat,'BinLimits',[0.5 1],'Normalization','probability', ...
        'FaceColor',colors{4},'EdgeColor','k');
    top = max(h.Values);
    lc = plot([0.5 0.5],[0 top],'--','color',colors{1},'LineWidth',2);
    lm = plot(mean(stat,'omitnan')*[1 1],[0 top],'--','color',colors{5},'LineWidth',2);
    xlabel(sprintf('Mean pairwise AUC across %u features',numel(ifeat)))
    ylabel('Probability')
    legend([h,lc,lm],{'features','no separation','mean'});
    % class distributions of the top features (normalised values)
    D = struct();
    [D.TS_DataMat,D.TimeSeries,D.Operations] = TS_LoadData(normFile);
    [D.TS_DataMat,D.TimeSeries] = FilterLabeledTimeSeries(D.TS_DataMat,D.TimeSeries);
    D.TimeSeries.Group = removecats(D.TimeSeries.Group);
    nd = min([16,opt.NumTopFeatures,numel(ifeat)]);
    f = figure('color','w');
    f.Position(3:4) = [1353,857];
    for i = 1:nd
        subplot(ceil(nd/4),4,i);
        TS_SingleFeature(D,sep.ID(ifeat(i)),true,false,sprintf('MeanAUC %.3f',stat(ifeat(i))),false);
    end
catch err
    fprintf('    (ranking figures not drawn: %s)\n',err.message);
end
end

%===============================================================================
function [red,reducedFile] = localReduce(normFile,ifeat,statName,opt)
% 'ReduceFeatureSet' (2.5.1): clusters every feature of HCTSA_N.mat with
% TS_ReduceFeatureSet (complete linkage on 1 - |Spearman rho|, cut at
% 1 - RedundancyThreshold), keeps one feature per cluster (localRepresentatives)
% and writes HCTSA_N_reduced.mat with TS_Subset: HCTSA_N.mat with only the kept
% features, every other variable carried over, plus the table red.
[~,clu] = TS_ReduceFeatureSet(normFile,1 - opt.RedundancyThreshold,true,false);
S = load(normFile);
rankN = nan(height(S.Operations),1);
rankN(ifeat) = (1:numel(ifeat))';
[~,pos] = ismember(clu.ID,S.Operations.ID);
rnk = rankN(pos);
rep = localRepresentatives(clu.Cluster,clu.Central,rnk);

red = table();
red.ID = clu.ID;
red.Name = string(clu.Name);
red.Cluster = clu.Cluster;           % 1 = the largest
red.ClusterSize = clu.ClusterSize;
red.Central = clu.Central;           % TS_ReduceFeatureSet's choice
red.Rank = rnk;
red.Kept = rep == (1:height(clu))';
red.RepresentativeID = clu.ID(rep);
red.Representative = red.Name(rep);
if isempty(ifeat)
    statName = '';
    rule = 'its most central member (the run has no ranking)';
else
    rule = sprintf('its best-ranked member (by %s), or its most central where none is ranked',statName);
end
[~,nm,ext] = fileparts(normFile);
red.Properties.Description = sprintf(['Complete-linkage clusters of the features of %s on ' ...
    '1 - |Spearman rho|, cut at %g: every pair within a cluster has |rho| >= %g. ' ...
    'Kept for each cluster: %s.'],[nm ext],1 - opt.RedundancyThreshold,opt.RedundancyThreshold,rule);
red.Properties.UserData = struct('RedundancyThreshold',opt.RedundancyThreshold, ...
    'distThreshold',1 - opt.RedundancyThreshold,'rankedBy',statName);
S.reducedFeatureSet = red;

reducedFile = regexprep(normFile,'\.mat$','_reduced.mat');
TS_Subset(S,[],red.ID(red.Kept),true,reducedFile);
fprintf('    %u features -> %u clusters (largest %u, %u single features); kept %s\n', ...
    height(red),sum(red.Kept),max(red.ClusterSize),sum(red.ClusterSize == 1),rule);
if ~isempty(ifeat)
    fprintf('    (%u of the kept features are not their cluster''s most central member)\n', ...
        sum(red.Kept & ~red.Central));
end
fprintf('    %s: %u features\n',reducedFile,sum(red.Kept));
end

%===============================================================================
function rep = localRepresentatives(cluster,central,rnk)
% The feature kept for each cluster: its best-ranked member (lowest rnk), or its
% central member when no member is ranked (rnk all NaN). rep(k) is the index of
% the feature kept for feature k's cluster. (Plain arrays: the Octave check runs it.)
rep = zeros(size(cluster));
for c = unique(cluster(:))'
    m = find(cluster == c);
    [best,i] = min(rnk(m));          % (min skips NaN)
    if isnan(best)
        i = find(central(m),1);
    end
    rep(m) = m(i);
end
end

%===============================================================================
function [FeatureTable,SeriesTable,DroppedSeries,FeatureAUC] = localTables(hctsaFile, ...
                                        normFile,ifeat,testStat,opt,droppedIngest,sep,reduced)
% Raw values (every series and feature of the run) and normalised values (what
% survived TS_Normalize).
R = load(hctsaFile,'TS_DataMat','TS_Quality','Operations','TimeSeries');
X = R.TS_DataMat; X(~isfinite(X)) = NaN; X(R.TS_Quality > 0) = NaN;
Ops = R.Operations; TS = R.TimeSeries;
[nTS,nOps] = size(X);
N = load(normFile,'TS_DataMat','TS_Quality','Operations','TimeSeries');
XN = N.TS_DataMat; XN(N.TS_Quality > 0) = NaN;
[inNorm,posN] = ismember(Ops.ID,N.Operations.ID);
keptRow = ismember(TS.ID,N.TimeSeries.ID);
grp = TS.Group;
labelled = keptRow & ~isundefined(grp);
cats = categories(removecats(grp(labelled)));

%--- FeatureTable: one row per feature of the run
FeatureTable = table();
FeatureTable.Name = string(Ops.Name);
FeatureTable.Column = repmat("",nOps,1);
try FeatureTable.Keywords = string(Ops.Keywords); catch, end
FeatureTable.ID = Ops.ID;
FeatureTable.Status = repmat("",nOps,1);
FeatureTable.PctComputed = 100*mean(~isnan(X),1)';

rank = nan(nOps,1); stat = nan(nOps,1);
if ~isempty(ifeat)
    rankN = nan(height(N.Operations),1);
    rankN(ifeat) = (1:numel(ifeat))';
    rank(inNorm) = rankN(posN(inNorm));
    stat(inNorm) = testStat(posN(inNorm));
end
FeatureTable.Rank = rank;
FeatureTable.Separability = stat;

% Why each feature is in or out (the order of TS_Normalize's own filters)
Xk = X(keptRow,:);
goodFrac = mean(~isnan(Xk),1)';
sdAll = std(Xk,0,1,'omitnan')';
status = repmat("dropped: constant after normalisation",nOps,1);
status(sdAll < 10*eps) = "dropped: constant";
if opt.ClassVarFilter && ~isempty(cats)
    for c = 1:numel(cats)
        sel = keptRow & grp == cats{c};
        if sum(sel) >= 2
            status(std(X(sel,:),0,1,'omitnan')' < 10*eps & sdAll >= 10*eps) = ...
                "dropped: constant within " + string(cats{c});
        end
    end
end
status(goodFrac < opt.FilterThresholds(2)) = "dropped: computed on " + ...
    string(round(100*goodFrac(goodFrac < opt.FilterThresholds(2)))) + "% of series";
status(inNorm) = "kept";
status(isfinite(rank)) = "ranked";
FeatureTable.Status = status;

% Two classes: AUC (0.5 = no separation, 1 = perfect) and which class is higher
if numel(cats) == 2
    auc = nan(nOps,1); dirn = repmat("",nOps,1);
    for j = 1:nOps
        x = X(labelled,j); g = grp(labelled);
        ok = ~isnan(x); x = x(ok); g = g(ok);
        in1 = (g == cats{1}); n1 = sum(in1); n0 = numel(x) - n1;
        if n1 == 0 || n0 == 0, continue; end
        r = tiedrank(x);
        a = (sum(r(in1)) - n1*(n1+1)/2)/(n1*n0);   % P(class 1 > class 2)
        auc(j) = max(a,1-a);
        if a > 0.5, dirn(j) = string(cats{1}); elseif a < 0.5, dirn(j) = string(cats{2}); end
    end
    FeatureTable.AUC = auc;
    FeatureTable.HigherIn = dirn;
end

% Per-group statistics: normalised mean, raw median and IQR
for c = 1:numel(cats)
    sel = labelled & grp == cats{c};
    m = nan(nOps,1);
    selN = ismember(N.TimeSeries.ID,TS.ID(sel));
    m(inNorm) = mean(XN(selN,posN(inNorm)),1,'omitnan')';
    nm = matlab.lang.makeUniqueStrings(matlab.lang.makeValidName( ...
        {['mean_' cats{c}],['median_' cats{c}],['iqr_' cats{c}]}), ...
        FeatureTable.Properties.VariableNames,63);
    FeatureTable.(nm{1}) = m;
    FeatureTable.(nm{2}) = median(X(sel,:),1,'omitnan')';
    FeatureTable.(nm{3}) = iqr(X(sel,:),1)';
end

% Redundancy among the top features: the better-ranked feature each one
% duplicates (|Spearman rho| >= RedundancyThreshold), or "" if none
if opt.RedundancyThreshold > 0 && ~isempty(ifeat)
    K = min(opt.RedundancyTopN,numel(ifeat));
    top = ifeat(1:K);
    rho = abs(corr(XN(:,top),'type','Spearman','rows','pairwise'));
    redundantWith = repmat("",nOps,1); redundantRho = nan(nOps,1);
    distinct = repmat("not checked",nOps,1);
    reps = [];
    for k = 1:K
        r = find(rho(k,reps) >= opt.RedundancyThreshold,1);
        j = find(inNorm & posN == top(k),1);
        if isempty(r)
            reps(end+1) = k; %#ok<AGROW>
            distinct(j) = "yes";
        else
            redundantWith(j) = string(N.Operations.Name{top(reps(r))});
            redundantRho(j) = rho(k,reps(r));
            distinct(j) = "no";
        end
    end
    FeatureTable.RedundantWith = redundantWith;
    FeatureTable.RedundancyRho = redundantRho;
    FeatureTable.Distinct = distinct;
    fprintf('    top %u features: %u distinct, %u redundant with a better-ranked one (|rho| >= %g)\n', ...
        K,numel(reps),K-numel(reps),opt.RedundancyThreshold);
end

% The whole feature set's clusters ('ReduceFeatureSet'): each feature's cluster,
% its size, and the feature kept for it ("" for that feature itself)
if ~isempty(reduced)
    [inRed,posR] = ismember(Ops.ID,reduced.ID);
    clNum = nan(nOps,1); clSize = nan(nOps,1); clRep = repmat("",nOps,1);
    clNum(inRed) = reduced.Cluster(posR(inRed));
    clSize(inRed) = reduced.ClusterSize(posR(inRed));
    clRep(inRed) = reduced.Representative(posR(inRed));
    clRep(ismember(Ops.ID,reduced.ID(reduced.Kept))) = "";
    FeatureTable.Cluster = clNum;
    FeatureTable.ClusterSize = clSize;
    FeatureTable.Representative = clRep;
end

%--- FeatureAUC: how well each kept feature separates the classes (see the help)
FeatureAUC = table();
if ~isempty(sep)
    [~,jF] = ismember(sep.ID,FeatureTable.ID);
    cl = sep.classes;
    newNames = @(c,T) matlab.lang.makeUniqueStrings(matlab.lang.makeValidName(c), ...
        T.Properties.VariableNames,63);
    FeatureAUC.Rank = FeatureTable.Rank(jF);
    FeatureAUC.Name = FeatureTable.Name(jF);
    FeatureAUC.ID = sep.ID;
    FeatureAUC.MeanAUC = sep.meanAUC';
    FeatureAUC.MinAUC = sep.minAUC';
    FeatureAUC.KW_p = sep.pKW';
    FeatureAUC.KW_q = localBH(sep.pKW)';
    q = localBH(sep.p);
    for k = 1:size(sep.pairs,1)
        a = cl{sep.pairs(k,1)}; b = cl{sep.pairs(k,2)};
        nm = newNames({['AUC_' a '_vs_' b],['Higher_' a '_vs_' b],['q_' a '_vs_' b]},FeatureAUC);
        FeatureAUC.(nm{1}) = sep.auc(k,:)';
        FeatureAUC.(nm{2}) = localWhich(sep.higher(k,:)',a,b);
        FeatureAUC.(nm{3}) = q(k,:)';
    end
    if numel(cl) > 2                    % (with two classes these repeat the pair)
        for k = 1:numel(cl)
            nm = newNames({['AUC_' cl{k} '_vs_rest'],['Higher_' cl{k} '_vs_rest']},FeatureAUC);
            FeatureAUC.(nm{1}) = sep.aucRest(k,:)';
            FeatureAUC.(nm{2}) = localWhich(sep.higherRest(k,:)',cl{k},'rest');
        end
    end
    FeatureAUC.PctComputed = 100*sum(sep.n,1)'/sep.numSeries;
    for k = 1:numel(cl)
        nm = newNames({['N_' cl{k}]},FeatureAUC);
        FeatureAUC.(nm{1}) = sep.n(k,:)';
    end
    for v = {'RedundantWith','RedundancyRho','Distinct','Cluster','ClusterSize','Representative'}
        if ismember(v{1},FeatureTable.Properties.VariableNames)
            FeatureAUC.(v{1}) = FeatureTable.(v{1})(jF);
        end
    end
    FeatureAUC = sortrows(FeatureAUC,'Rank','ascend','MissingPlacement','last');
end

%--- SeriesTable: one row per time series that survived normalisation
SeriesTable = table();
SeriesTable.Name = string(N.TimeSeries.Name);
SeriesTable.Group = string(N.TimeSeries.Group);
SeriesTable.Group(ismissing(SeriesTable.Group)) = "";
try SeriesTable.Length = N.TimeSeries.Length; catch, end
nOpsN = height(N.Operations);
keep = 1:nOpsN;
if ~isfinite(opt.MaxFeatures) && height(SeriesTable)*nOpsN > 5e7
    fprintf(2,'    SeriesTable will be %u x %u (~%.1f GB). Set ''MaxFeatures'' to trim.\n', ...
        height(SeriesTable),nOpsN,8*height(SeriesTable)*nOpsN/1e9);
end
if isfinite(opt.MaxFeatures) && opt.MaxFeatures < nOpsN
    if ~isempty(ifeat)
        keep = ifeat(1:min(opt.MaxFeatures,numel(ifeat)));
    else
        keep = 1:opt.MaxFeatures;
    end
    fprintf('    (SeriesTable restricted to %u features)\n',numel(keep));
end
cols = localColumnNames(N.Operations.Name(keep),SeriesTable.Properties.VariableNames);
SeriesTable = [SeriesTable, array2table(XN(:,keep),'VariableNames',cols)];
[~,jRaw] = ismember(N.Operations.ID(keep),Ops.ID);
FeatureTable.Column(jRaw) = string(cols);

FeatureTable = sortrows(FeatureTable,'Rank','ascend','MissingPlacement','last');

%--- DroppedSeries: at ingest, and by normalisation
dropN = find(~keptRow);
DroppedSeries = table();
DroppedSeries.Name = [string(droppedIngest.Name(:)); string(TS.Name(dropN))];
DroppedSeries.Stage = [repmat("ingest",numel(droppedIngest.Name),1); ...
                       repmat("normalisation",numel(dropN),1)];
DroppedSeries.Reason = [string(droppedIngest.Reason(:)); ...
    "only " + string(round(100*mean(~isnan(X(dropN,:)),2))) + ...
    "% of features computed (needs " + string(100*opt.FilterThresholds(1)) + "%)"];
end

%===============================================================================
function cols = localColumnNames(names,reserved)
% Valid, unique table-variable names of at most 63 characters. A long name keeps
% its first 24 and last 38 characters: the end is what tells the many
% MF_GP_* variants apart (plain truncation made 66 of them collide).
cols = regexprep(cellstr(names(:))','[^A-Za-z0-9_]','_');
for i = 1:numel(cols)
    if isempty(regexp(cols{i},'^[A-Za-z]','once')), cols{i} = ['x' cols{i}]; end
    if numel(cols{i}) > 63
        cols{i} = [cols{i}(1:24) '_' cols{i}(end-37:end)];
    end
end
cols = matlab.lang.makeUniqueStrings(cols,reserved,63);
end

%===============================================================================
function s = localWhich(h,a,b)
% a where h > 0, b where h < 0, "" otherwise (a tie, or not computed).
s = repmat("",numel(h),1);
s(h > 0) = string(a);
s(h < 0) = string(b);
end

%===============================================================================
function localWriteCsv(T,f)
try
    writetable(T,f);
catch err
    fprintf(2,'    Could not write %s (%s) -- is it open in another program?\n',f,err.message);
end
end

%===============================================================================
function localSaveFigures(figs,opt)
nSaved = 0; used = {};
for i = 1:numel(figs)
    try
        nm = get(figs(i),'Name');
        if isempty(nm)
            nm = sprintf('figure_%u',i);
        else
            nm = matlab.lang.makeValidName(nm);
        end
        nm = matlab.lang.makeUniqueStrings(nm,used);
        used{end+1} = nm; %#ok<AGROW>
        base = fullfile(opt.OutputDir,nm);
        if ismember(opt.FigureFormat,{'png','both'}), saveas(figs(i),[base '.png']); end
        if ismember(opt.FigureFormat,{'fig','both'}), savefig(figs(i),[base '.fig']); end
        nSaved = nSaved + 1;
    catch err
        fprintf('    (figure not saved: %s)\n',err.message);
    end
end
if nSaved > 0
    fprintf('    saved %u figure(s)\n',nSaved);
end
if opt.CloseFigures
    close(figs(isgraphics(figs)));
end
end

%===============================================================================
function info = localRunInfo(opt,steps,t0)
info = struct();
info.options = opt;
info.date = datestr(now);
info.matlab = version;
info.computer = computer;
try
    g = TS_AddGitInfo();
    info.git = g;
catch
    info.git = struct();
end
info.steps = steps;
info.seconds = toc(t0);
end
