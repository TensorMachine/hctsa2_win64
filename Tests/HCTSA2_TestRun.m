function HCTSA2_TestRun(outDir)
% HCTSA2_TestRun   Plumbing check for HCTSA2_Analyse on a real MATLAB installation.
%
%   >> HCTSA2_TestRun              % output in Tests/testrun_output
%   >> HCTSA2_TestRun(outDir)      % output in a folder of your choice
%
% Run after startup (which puts Tests/ on the path), from any folder. It reads
% hctsa2_testdata.mat beside this file and writes its testrun_* folders into its
% output folder, deleting the ones a previous run left there. (HCTSA2 3.0.0: it
% lives in Tests/ and is a function; the default output folder, Tests/testrun_output,
% is not tracked by git, so neither the hctsa2 top level nor the current folder
% fills up.)
%
% Uses hctsa2_testdata.mat: 12 synthetic series (~300 points each) in three
% classes with deliberately distinct dynamics --
%
%   class 1  white noise            lag-1 autocorrelation ~0,     first zero of ac = 1
%   class 2  AR(1), phi = 0.9       lag-1 autocorrelation ~0.9,   first zero  21-31
%   class 3  noisy sine, period 24  lag-1 autocorrelation ~0.8,   first zero   6-7
%
% Classes 2 and 3 are both strongly autocorrelated on purpose: no single feature
% separates all three, so the ranking stage has to actually work rather than
% latching onto one obvious statistic.
%
% The file stores a LONG-FORMAT NUMERIC MATRIX (id, time, value, class) with its
% rows shuffled, which is what a real export usually looks like.
%
%   stage 1  the matrix, catch22, plots on; any keyboard prompt is a FAIL.
%            FeatureAUC is recomputed here, independently (pairwise and
%            one-vs-rest AUCs by counting, ranksum, kruskalwallis, BH q-values),
%            and the ranking must follow MeanAUC
%   stage 2  the same data as a MATLAB table (different ingestion code), ranked
%            by 'classification': the values must equal stage 1's
%   stage 3  parallel extraction over series and over operations (Parallel
%            Computing Toolbox only): the values must equal stage 1's
%   stage 4  two classes with separating features in BOTH directions, passed as
%            a matrix in memory: every strong separator must rank at the top;
%            FeatureAUC's single pair must equal FeatureTable's AUC
%   stage 5  unbalanced classes (30 and 8 series) in a NaN-padded matrix: the
%            padding is trimmed and the accuracy is a number, not NaN
%   stage 6  resume: half of a finished run is blanked out, and the same call
%            must resume it and reproduce the values
%   stage 7  'Features' and 'ExcludeFeatures': only the selected features are
%            computed
%   stage 8  'ExistingHCTSA' with new group labels: the input file is untouched
%   stage 9  'ReduceFeatureSet' on stage 4's file, with one feature made constant
%            (dropped by normalisation) and two missing a fifth of their values:
%            HCTSA_N_reduced.mat holds exactly one feature per cluster with
%            HCTSA_N.mat's values; the clusters equal an independent complete
%            linkage on |Spearman rho| computed pair by pair on the shared series;
%            the kept feature is the best-ranked member (with classes) or the most
%            central (one class); the option changes nothing else;
%            TS_ReduceFeatureSet's figure gives the same clusters, and a single
%            feature is one cluster
%   stage 10 HCTSA2_ComputeFeatures: stage 4's catch22 values for two of its
%            series come back bit for bit, by ID, by name and by field name, in the
%            requested order and with NaN padding; NaN inside a series gives NaN;
%            unknown features are errors; the full library's 6,392 field names;
%            (3.0.0) 'NormalizeLike' with stage 4's run: its two series come back
%            with their HCTSA_N.mat values, bit for bit, and a run without the
%            stored normalisation, or a feature it filtered out, is an error;
%            features with known values (mean, median, standard deviation) and two
%            outputs of one master operation equal a direct call
%   stage 11 (3.0.0) TS_TopFeatures with 5 label-shuffled nulls on stage 4's file: the
%            p-values are in (0,1] and the q-values (BF_mafdr, in place of the
%            Bioinformatics Toolbox's mafdr) equal the Benjamini-Hochberg definition,
%            computed here pair by pair; the global random stream is restored
%
% The accuracy printed in stages 1-3 is at chance (about 30%) because 12 series
% are too few for HCTSA2_Analyse's default decision trees to split (see
% 'Classifier' in its help). This script checks the plumbing, not the accuracy.
%
% Timing: about a minute. Stages 1-9 compute the catch22 set; stage 10 reads the
% full library's tables and computes five cheap features of it; stage 11 ranks
% stage 4's catch22 features six times. The full library
% is exercised by the Windows validation package.

fprintf('\n============ HCTSA2 pipeline check ============\n');

testDir = fileparts(mfilename('fullpath'));
dataFile = fullfile(testDir,'hctsa2_testdata.mat');
if exist(dataFile,'file') ~= 2
    error('HCTSA2_TestRun:noData','%s not found (it belongs beside this file).',dataFile);
end
if nargin < 1 || isempty(outDir)
    outDir = fullfile(testDir,'testrun_output');
end
if ~isfolder(outDir), mkdir(outDir); end
startDir = cd(outDir);
restoreDir = onCleanup(@() cd(startDir)); %#ok<NASGU>
fprintf('Output folder: %s\n',pwd);
for d = {'testrun_catch22','testrun_table','testrun_parallel_series', ...
         'testrun_parallel_operations','testrun_twoclass','testrun_unbalanced', ...
         'testrun_resume','testrun_select','testrun_existing','testrun_reduce_src', ...
         'testrun_reduce','testrun_reduce_off','testrun_reduce_oneclass'}
    if exist(d{1},'dir') == 7, rmdir(d{1},'s'); end
end
nFail = 0;
rs = RandStream('mt19937ar','Seed',1);   % this script's data; the global stream is untouched

%% ---------------------------------------------------------------------------
%% Stage 1 -- long-format numeric matrix, catch22 (fast); no prompts allowed
%% ---------------------------------------------------------------------------
fprintf('\n--- stage 1: long-format matrix, catch22 feature set ---\n');
% Any input() reached anywhere below now throws (HCTSA2_Analyse must not ask).
promptDir = tempname; mkdir(promptDir);
fid = fopen(fullfile(promptDir,'input.m'),'w');
fprintf(fid,'function varargout = input(varargin)\nerror(''HCTSA2_TestRun:prompt'',''HCTSA2_Analyse asked for keyboard input'');\n');
fclose(fid);
w = warning('off','all');     % the shadowing is the point: no "same name as a built-in"
addpath(promptDir);
warning(w);
t0 = tic;
try
    R1 = HCTSA2_Analyse(dataFile,'FeatureSet','catch22','OutputDir','testrun_catch22','DoPlots',true);
catch err
    R1 = [];
    fprintf(2,'    FAIL: %s\n',err.message);
    nFail = nFail + 1;
end
rmpath(promptDir); rmdir(promptDir,'s');
fprintf('\n  stage 1 finished in %.0f s\n',toc(t0));
nFail = nFail + checkResult(R1,'stage 1',12,3);
nFail = nFail + checkAUC(R1,'stage 1','meanAUC');

%% ---------------------------------------------------------------------------
%% Stage 2 -- the same data as a MATLAB table
%% ---------------------------------------------------------------------------
% This is the format most real data arrives in, and it reaches different code in
% HCTSA2_Ingest than the matrix above, so it is worth running both.
fprintf('\n--- stage 2: same data as a table ---\n');
S = load(dataFile);
T = table(S.LongData(:,1),S.LongData(:,2),S.LongData(:,3),S.LongData(:,4), ...
          'VariableNames',{'id','time','value','class'});
save('hctsa2_testdata_table.mat','T');
R2 = HCTSA2_Analyse('hctsa2_testdata_table.mat','FeatureSet','catch22', ...
        'OutputDir','testrun_table','DoPlots',false,'TestStatistic','classification');
nFail = nFail + checkResult(R2,'stage 2',12,3);
% The two runs used identical numbers by two different routes, so their feature
% values should agree exactly. If they do not, one of the ingestion paths is
% mangling the data -- which is precisely the sort of fault that produces
% plausible-looking output and no error.
nFail = nFail + compareRuns(R1,R2,'the two ingestion routes');

%% ---------------------------------------------------------------------------
%% Stage 3 -- parallel extraction, both modes (Parallel Computing Toolbox only)
%% ---------------------------------------------------------------------------
% HCTSA2_Analyse passes 'ParallelMode' down to TS_Compute. Either way the values
% must be stage 1's: catch22 is deterministic compiled code.
if license('test','Distrib_Computing_Toolbox') && exist('gcp','file') == 2
    for mode = {'series','operations'}
        fprintf('\n--- stage 3: parallel extraction over %s ---\n',mode{1});
        R3 = HCTSA2_Analyse(dataFile,'FeatureSet','catch22', ...
                'OutputDir',['testrun_parallel_' mode{1}],'DoPlots',false, ...
                'Parallel',true,'ParallelMode',mode{1});
        nFail = nFail + checkResult(R3,['stage 3, ' mode{1}],12,3);
        nFail = nFail + compareRuns(R1,R3,['serial and parallel (' mode{1} ') extraction']);
    end
else
    fprintf('\n--- stage 3 skipped: no Parallel Computing Toolbox ---\n');
end

%% ---------------------------------------------------------------------------
%% Stage 4 -- two classes, separators in both directions (in-memory matrix)
%% ---------------------------------------------------------------------------
% AR(1) with phi = +0.6 against phi = -0.6: autocorrelation features are higher
% in class A, increment-based ones (e.g. successive differences) in class B.
% Before 2.2.1 the default two-class statistic ranked one direction last.
fprintf('\n--- stage 4: two classes, separators in both directions ---\n');
X4 = [arSeries(rs,0.6,20,250); arSeries(rs,-0.6,20,250)];
g4 = [repmat({'A'},20,1); repmat({'B'},20,1)];
R4 = HCTSA2_Analyse(X4,'Groups',g4,'FeatureSet','catch22','OutputDir','testrun_twoclass', ...
        'DoPlots',false);
nFail = nFail + checkResult(R4,'stage 4',40,2);
F = R4.FeatureTable;
strong = F.AUC >= 0.95;
dirs = unique(F.HigherIn(strong));
if numel(dirs) < 2
    fprintf(2,'    FAIL: the data did not give strong separators in both directions (%s)\n', ...
        strjoin(cellstr(dirs),', '));
    nFail = nFail + 1;
elseif max(F.Rank(strong)) > sum(F.AUC >= 0.8)
    fprintf(2,'    FAIL: a strong separator (AUC >= 0.95) ranks %u, below weaker ones\n', ...
        max(F.Rank(strong)));
    nFail = nFail + 1;
else
    fprintf('    OK  : %u strong separators (higher in %s), all ranked in the top %u\n', ...
        sum(strong),strjoin(cellstr(dirs),' and '),max(F.Rank(strong)));
end
nFail = nFail + checkAUC(R4,'stage 4','ustatP');
A4 = R4.FeatureAUC;
[~,j4] = ismember(A4.ID,F.ID);
if isequaln(A4.AUC_A_vs_B,F.AUC(j4)) && isequal(A4.Higher_A_vs_B,F.HigherIn(j4)) && ...
        ~any(endsWith(A4.Properties.VariableNames,'_vs_rest'))
    fprintf('    OK  : FeatureAUC''s A vs B equals FeatureTable''s AUC and HigherIn\n');
else
    fprintf(2,'    FAIL: FeatureAUC''s A vs B differs from FeatureTable''s AUC\n'); nFail = nFail + 1;
end

%% ---------------------------------------------------------------------------
%% Stage 5 -- unbalanced classes in a NaN-padded matrix
%% ---------------------------------------------------------------------------
fprintf('\n--- stage 5: 30 and 8 series, ragged lengths padded with NaN ---\n');
X5 = [arSeries(rs,0.6,30,250); arSeries(rs,-0.6,8,250)];
for i = 1:5:38, X5(i,200:end) = NaN; end      % shorter series, NaN-padded
g5 = [ones(30,1); 2*ones(8,1)];
R5 = HCTSA2_Analyse(X5,'Groups',g5,'FeatureSet','catch22','OutputDir','testrun_unbalanced', ...
        'DoPlots',false);
nFail = nFail + checkResult(R5,'stage 5',38,2);
if isfinite(R5.classificationAccuracy)
    fprintf('    OK  : %s %.1f%% with %u folds\n',R5.classification.measure, ...
        R5.classificationAccuracy,R5.classification.numFolds);
else
    fprintf(2,'    FAIL: classification accuracy is NaN\n'); nFail = nFail + 1;
end
if height(R5.DroppedSeries) == 0 && all(R5.SeriesTable.Length(1:5:38) == 199)
    fprintf('    OK  : padding trimmed, no series dropped\n');
else
    fprintf(2,'    FAIL: padded series were dropped or not trimmed\n'); nFail = nFail + 1;
end

%% ---------------------------------------------------------------------------
%% Stage 6 -- resume an interrupted extraction
%% ---------------------------------------------------------------------------
fprintf('\n--- stage 6: resume ---\n');
args6 = {'Groups',g4,'FeatureSet','catch22','OutputDir','testrun_resume','DoPlots',false};
R6a = HCTSA2_Analyse(X4,args6{:});
f6 = fullfile('testrun_resume','HCTSA.mat');
M = load(f6,'TS_DataMat','TS_Quality','TS_CalcTime');
M.TS_DataMat(21:end,:) = NaN; M.TS_Quality(21:end,:) = NaN; M.TS_CalcTime(21:end,:) = NaN;
save(f6,'-struct','M','-append');                % as if stopped after 20 series
log6 = evalc('R6b = HCTSA2_Analyse(X4,args6{:});');
if contains(log6,'Resuming') && contains(log6,'20 of 40 series to compute')
    fprintf('    OK  : resumed, computing only the 20 missing series\n');
else
    fprintf(2,'    FAIL: the second call did not resume the 20 missing series\n'); nFail = nFail + 1;
end
nFail = nFail + compareRuns(R6a,R6b,'the complete and the resumed run');

%% ---------------------------------------------------------------------------
%% Stage 7 -- feature selection before extraction
%% ---------------------------------------------------------------------------
fprintf('\n--- stage 7: Features and ExcludeFeatures ---\n');
R7 = HCTSA2_Analyse(dataFile,'FeatureSet','catch22','OutputDir','testrun_select', ...
        'DoPlots',false,'Features',{'mode_*','outlier_*','acf_timescale'}, ...
        'ExcludeFeatures','mode_10');
ops7 = TS_GetFromData(fullfile('testrun_select','HCTSA.mat'),'Operations');
want7 = {'mode_5','outlier_timing_neg','outlier_timing_pos','acf_timescale'};
if isequal(sort(ops7.Name(:)),sort(want7(:)))
    fprintf('    OK  : computed exactly the %u selected features\n',numel(want7));
else
    fprintf(2,'    FAIL: computed %s\n',strjoin(ops7.Name,', ')); nFail = nFail + 1;
end

%% ---------------------------------------------------------------------------
%% Stage 8 -- ExistingHCTSA with new labels leaves the input file untouched
%% ---------------------------------------------------------------------------
fprintf('\n--- stage 8: ExistingHCTSA, relabelled ---\n');
src8 = fullfile('testrun_twoclass','HCTSA.mat');
before = dir(src8); contentBefore = load(src8);
g8 = [repmat({'P'},10,1); repmat({'Q'},30,1)];
R8 = HCTSA2_Analyse([],'ExistingHCTSA',src8,'Groups',g8,'OutputDir','testrun_existing', ...
        'DoPlots',false);
after = dir(src8);
if before.bytes == after.bytes && before.datenum == after.datenum && ...
        isequaln(contentBefore,load(src8))
    fprintf('    OK  : %s unchanged\n',src8);
else
    fprintf(2,'    FAIL: %s was modified\n',src8); nFail = nFail + 1;
end
if isequal(sort(R8.groupNames),{'P','Q'}) && isfile(fullfile('testrun_existing','HCTSA_N.mat'))
    fprintf('    OK  : relabelled to P/Q; all output in testrun_existing\n');
else
    fprintf(2,'    FAIL: new labels not used, or output written elsewhere\n'); nFail = nFail + 1;
end

%% ---------------------------------------------------------------------------
%% Stage 9 -- 'ReduceFeatureSet': the whole feature set clustered (2.5.1)
%% ---------------------------------------------------------------------------
% Stage 4's file again, with the last feature made constant (normalisation drops
% it) and the first two missing on series 1-8 and 5-12 (kept: 80% computed), so
% that HCTSA_N.mat has missing values, as real runs do. The threshold is lowered
% to 0.5 so that catch22's 22 features form clusters of several members.
fprintf('\n--- stage 9: ReduceFeatureSet ---\n');
mkdir('testrun_reduce_src');
src9 = fullfile('testrun_reduce_src','HCTSA.mat');
M9 = load(fullfile('testrun_twoclass','HCTSA.mat'));
M9.TS_DataMat(:,end) = 0.5;
M9.TS_DataMat(1:8,1) = NaN;  M9.TS_Quality(1:8,1) = 1;
M9.TS_DataMat(5:12,2) = NaN; M9.TS_Quality(5:12,2) = 1;
save(src9,'-struct','M9','-v7.3');
args9 = {'ExistingHCTSA',src9,'DoPlots',false,'RedundancyThreshold',0.5};
R9 = HCTSA2_Analyse([],args9{:},'Groups',g4,'OutputDir','testrun_reduce','ReduceFeatureSet',true);
nFail = nFail + checkReduced(R9,'stage 9, two classes',0.5,true);
R9off = HCTSA2_Analyse([],args9{:},'Groups',g4,'OutputDir','testrun_reduce_off');
nFail = nFail + compareReduceOff(R9,R9off);
R9one = HCTSA2_Analyse([],args9{:},'Groups',repmat({'P'},40,1), ...
        'OutputDir','testrun_reduce_oneclass','ReduceFeatureSet',true);
nFail = nFail + checkReduced(R9one,'stage 9, one class',0.5,false);
for thr0 = [0 1]
    try
        HCTSA2_Analyse([],'ExistingHCTSA',src9,'Groups',g4,'OutputDir','testrun_reduce_off', ...
            'DoPlots',false,'ReduceFeatureSet',true,'RedundancyThreshold',thr0);
        fprintf(2,'    FAIL: ReduceFeatureSet with RedundancyThreshold %g was accepted\n',thr0);
        nFail = nFail + 1;
    catch err
        if strcmp(err.identifier,'HCTSA2_Analyse:reduce')
            fprintf('    OK  : ReduceFeatureSet with RedundancyThreshold %g is refused\n',thr0);
        else
            fprintf(2,'    FAIL: RedundancyThreshold %g: %s\n',thr0,err.message); nFail = nFail + 1;
        end
    end
end

%% ---------------------------------------------------------------------------
%% Stage 10 -- HCTSA2_ComputeFeatures: chosen features of one series (2.5.1)
%% ---------------------------------------------------------------------------
fprintf('\n--- stage 10: HCTSA2_ComputeFeatures ---\n');
nFail = nFail + checkComputeFeatures(R4,rs);

%% ---------------------------------------------------------------------------
%% Stage 11 -- TS_TopFeatures' p- and q-values with nulls (3.0.0)
%% ---------------------------------------------------------------------------
fprintf('\n--- stage 11: TS_TopFeatures with nulls ---\n');
nFail = nFail + checkTopFeaturesNulls(fullfile('testrun_twoclass','HCTSA_N.mat'));

%% ---------------------------------------------------------------------------
fprintf('\n-------------------------------------------------\n');
if nFail == 0
    fprintf('All stages passed.\n');
else
    fprintf(2,'%u check(s) FAILED (see above).\n',nFail);
end
fprintf('=================================================\n\n');
end

%% ---------------------------------------------------------------------------
function nFail = checkTopFeaturesNulls(fN)
% (3.0.0) Upstream's TS_TopFeatures, with 'numNulls' > 0, gives each feature a p-value
% from pooled label-shuffled nulls and adjusts them for the false discovery rate
% (Benjamini-Hochberg). Upstream calls the Bioinformatics Toolbox's mafdr for that; this
% fork calls BF_mafdr. The q-values must equal the definition, computed here pair by
% pair: q_i = min over p_j >= p_i of min(1, m*p_j/r_j), with r_j the number of p-values
% <= p_j. The global random stream, which the shuffles use, is restored afterwards.
nFail = 0;
s0 = rng;
restoreRng = onCleanup(@() rng(s0)); %#ok<NASGU>
try
    rng(7,'twister');
    [~,testStat,nulls,~,pVals,qVals] = TS_TopFeatures(fN,'ustatP',struct(), ...
        'whatPlots',{},'numNulls',5);
    ok = ~isnan(testStat(:));                  % TS_TopFeatures' own isValid
    p = pVals(ok); p = p(:); q = qVals(ok); q = q(:);
    m = numel(p);
    qb = zeros(m,1);
    for i = 1:m
        qb(i) = Inf;
        for j = 1:m
            if p(j) >= p(i)
                qb(i) = min(qb(i),min(1,m*p(j)/sum(p <= p(j))));
            end
        end
    end
    worst = max([0; abs(q - qb)]);
    good = m > 0 && size(nulls,2) == 5 && all(p > 0 & p <= 1) && ...
        all(isnan(pVals(~ok))) && all(isnan(qVals(~ok))) && worst < 1e-12 && ...
        all(q >= p*(1 - 1e-12));               % q >= p, up to the round-off of m*p/m
    if good
        fprintf(['    OK  : %u features with 5 nulls: p-values in (0,1], q-values equal the ' ...
            'Benjamini-Hochberg definition (worst %.2g), %u significant at q < 0.05\n'], ...
            m,worst,sum(q < 0.05));
    else
        fprintf(2,['    FAIL: p- or q-values wrong (%u features, %u nulls, worst q difference ' ...
            '%.3g)\n'],m,size(nulls,2),worst);
        nFail = 1;
    end
catch err
    fprintf(2,'    FAIL: TS_TopFeatures with nulls: %s\n',err.message);
    nFail = 1;
end
end

%% ---------------------------------------------------------------------------
function nFail = checkResult(R,label,nSeries,nClasses)
% Confirm the pipeline produced the tables, with usable content.
fprintf('\n  [%s] checking the returned structure\n',label);
nFail = 0;
if ~isstruct(R), fprintf(2,'    FAIL: no struct returned\n'); nFail = 1; return; end
for f = {'FeatureTable','SeriesTable','DroppedSeries'}
    if ~isfield(R,f{1}) || ~istable(R.(f{1}))
        fprintf(2,'    FAIL: %s missing or not a table\n',f{1}); nFail = nFail + 1;
    else
        fprintf('    OK  : %-13s %u rows x %u columns\n',f{1},height(R.(f{1})),width(R.(f{1})));
    end
end
if nFail > 0, return; end

n = height(R.SeriesTable);
nc = numel(unique(R.SeriesTable.Group));
if n == nSeries && nc == nClasses
    fprintf('    OK  : %u series in %u classes\n',n,nc);
else
    fprintf(2,'    FAIL: %u series in %u classes (expected %u in %u)\n',n,nc,nSeries,nClasses);
    nFail = nFail + 1;
end

% Every feature of the run is listed with a status, and the columns map back
isNum = varfun(@isnumeric,R.SeriesTable,'OutputFormat','uniform');
M = table2array(R.SeriesTable(:,isNum));
nanFrac = mean(~isfinite(M(:)));
if nanFrac > 0.2
    fprintf(2,'    WARN: %.0f%% of feature values are NaN or Inf\n',100*nanFrac);
end
F = R.FeatureTable;
mapped = F.Column ~= "";
if all(F.Status ~= "") && all(ismember(F.Column(mapped),R.SeriesTable.Properties.VariableNames))
    fprintf('    OK  : %u features listed (%u ranked); every Column is in SeriesTable\n', ...
        height(F),sum(isfinite(F.Rank)));
else
    fprintf(2,'    FAIL: FeatureTable Status/Column incomplete\n'); nFail = nFail + 1;
end
if isfinite(F.Rank(1))
    fprintf('    OK  : top-ranked feature is %s\n',F.Name(1));
end

% FeatureAUC: a row per kept feature and an AUC column per pair; also as a CSV
try
    A = R.FeatureAUC;
    v = A.Properties.VariableNames;
    nPairs = sum(startsWith(v,'AUC_') & ~endsWith(v,'_vs_rest'));
    csv = fullfile(R.options.OutputDir,'FeatureAUC.csv');
    nCsv = NaN;
    if isfile(csv), nCsv = height(readtable(csv)); end
    [~,jA] = ismember(A.ID,F.ID);
    nKept = sum(ismember(F.Status,["kept","ranked"]));
    sameRank = all(jA > 0) && isequaln(A.Rank,F.Rank(jA));
    if height(A) == nKept && nPairs == nClasses*(nClasses-1)/2 && nCsv == height(A) && sameRank
        fprintf('    OK  : FeatureAUC %u rows x %u columns (%u pairs), the same ranks; FeatureAUC.csv\n', ...
            height(A),width(A),nPairs);
    else
        fprintf(2,['    FAIL: FeatureAUC has %u rows (%u kept features), %u pair columns, ' ...
            '%g CSV rows; ranks as in FeatureTable: %u\n'],height(A),nKept,nPairs,nCsv,sameRank);
        nFail = nFail + 1;
    end
catch err
    fprintf(2,'    FAIL: FeatureAUC: %s\n',err.message); nFail = nFail + 1;
end
end

%% ---------------------------------------------------------------------------
function nFail = checkAUC(R,label,statName)
% FeatureAUC recomputed from HCTSA.mat by other means: AUCs by counting pairs,
% p-values from the Statistics Toolbox (ranksum 'approx', kruskalwallis), q-values
% by the definition of Benjamini-Hochberg; and the ranking and redundancy columns.
fprintf('\n  [%s] checking FeatureAUC against an independent computation\n',label);
try
    nFail = aucChecks(R,statName);
catch err
    fprintf(2,'    FAIL: %s\n',err.message);
    nFail = 1;
end
end

function nFail = aucChecks(R,statName)
nFail = 0;
A = R.FeatureAUC;
H = load(R.files.hctsaFile,'TS_DataMat','TS_Quality','Operations','TimeSeries');
N = load(R.files.normFile,'TimeSeries');
rows = ismember(H.TimeSeries.ID,N.TimeSeries.ID) & ~isundefined(H.TimeSeries.Group);
g = removecats(H.TimeSeries.Group(rows));
cl = categories(g); K = numel(cl);
[~,col] = ismember(A.ID,H.Operations.ID);
X = H.TS_DataMat(rows,col);
X(~isfinite(X) | H.TS_Quality(rows,col) > 0) = NaN;
pairs = nchoosek(1:K,2);
nF = height(A); nP = size(pairs,1);
auc = nan(nF,nP); hi = strings(nF,nP); p = nan(nF,nP);
aucR = nan(nF,K); hiR = strings(nF,K); pKW = nan(nF,1); n = zeros(nF,K);
for j = 1:nF
    x = X(:,j);
    for k = 1:K, n(j,k) = sum(g == cl{k} & ~isnan(x)); end
    for q = 1:nP
        xa = x(g == cl{pairs(q,1)} & ~isnan(x)); xb = x(g == cl{pairs(q,2)} & ~isnan(x));
        if isempty(xa) || isempty(xb), continue; end
        u = sum(sum((xa > xb') + (xa == xb')/2))/(numel(xa)*numel(xb));   % exact
        auc(j,q) = max(u,1-u);
        if u > 0.5, hi(j,q) = cl{pairs(q,1)}; elseif u < 0.5, hi(j,q) = cl{pairs(q,2)}; else, hi(j,q) = ""; end
        if all([xa; xb] == xa(1)), p(j,q) = 1; else, p(j,q) = ranksum(xa,xb,'method','approx'); end
    end
    for k = 1:K
        xa = x(g == cl{k} & ~isnan(x)); xb = x(g ~= cl{k} & ~isnan(x));
        if isempty(xa) || isempty(xb), continue; end
        u = sum(sum((xa > xb') + (xa == xb')/2))/(numel(xa)*numel(xb));
        aucR(j,k) = max(u,1-u);
        if u > 0.5, hiR(j,k) = cl{k}; elseif u < 0.5, hiR(j,k) = "rest"; else, hiR(j,k) = ""; end
    end
    ok = ~isnan(x);
    if all(n(j,:) > 0)
        if all(x(ok) == x(find(ok,1))), pKW(j) = 1; else, pKW(j) = kruskalwallis(x(ok),g(ok),'off'); end
    end
end
err = @(a,b) max([0; abs(a(:) - b(:))./(abs(b(:)) + 1e-3)],[],'omitnan');   % < 1e-9: 1e-9 relative, 1e-12 absolute
same = @(a,b) isequal(isnan(a),isnan(b));
d = [];
names = {};
for q = 1:nP
    s = [cl{pairs(q,1)} '_vs_' cl{pairs(q,2)}];
    d(end+1) = ternary(same(A.(['AUC_' s]),auc(:,q)),err(A.(['AUC_' s]),auc(:,q)),Inf); %#ok<AGROW>
    d(end+1) = ternary(isequal(A.(['Higher_' s]),hi(:,q)),0,Inf); %#ok<AGROW>
    qq = bh(p(:,q));
    d(end+1) = ternary(same(A.(['q_' s]),qq),err(A.(['q_' s]),qq),Inf); %#ok<AGROW>
    names = [names, {['AUC_' s],['Higher_' s],['q_' s]}]; %#ok<AGROW>
end
if K > 2
    for k = 1:K
        s = [cl{k} '_vs_rest'];
        d(end+1) = ternary(same(A.(['AUC_' s]),aucR(:,k)),err(A.(['AUC_' s]),aucR(:,k)),Inf); %#ok<AGROW>
        d(end+1) = ternary(isequal(A.(['Higher_' s]),hiR(:,k)),0,Inf); %#ok<AGROW>
        names = [names, {['AUC_' s],['Higher_' s]}]; %#ok<AGROW>
    end
end
meanA = mean(auc,2); minA = min(auc,[],2); minA(any(isnan(auc),2)) = NaN;
checks = {'MeanAUC',A.MeanAUC,meanA; 'MinAUC',A.MinAUC,minA; 'KW_p',A.KW_p,pKW; ...
          'KW_q',A.KW_q,bh(pKW); 'PctComputed',A.PctComputed,100*sum(n,2)/numel(g)};
for i = 1:size(checks,1)
    d(end+1) = ternary(same(checks{i,2},checks{i,3}),err(checks{i,2},checks{i,3}),Inf); %#ok<AGROW>
    names{end+1} = checks{i,1}; %#ok<AGROW>
end
for k = 1:K
    d(end+1) = ternary(isequal(A.(['N_' cl{k}]),n(:,k)),0,Inf); %#ok<AGROW>
    names{end+1} = ['N_' cl{k}]; %#ok<AGROW>
end
if max(d) < 1e-9
    fprintf('    OK  : %u features x %u columns agree (largest relative difference %.1e)\n', ...
        nF,numel(d),max(d));
else
    fprintf(2,'    FAIL: these columns disagree: %s\n',strjoin(names(d >= 1e-9),', '));
    nFail = nFail + 1;
end

% The ranking: the expected statistic; by MeanAUC (ties by MinAUC, then ID) for 'meanAUC'
r = A(isfinite(A.Rank),:);
if ~strcmp(R.separabilityStatistic,statName)
    fprintf(2,'    FAIL: ranked by %s, expected %s\n',R.separabilityStatistic,statName);
    nFail = nFail + 1;
elseif strcmp(statName,'meanAUC')
    key = [-round(r.MeanAUC,12),-round(r.MinAUC,12),r.ID];
    if issorted(r.Rank) && isequal(sortrows(key),key) && height(r) == sum(isfinite(A.MeanAUC))
        fprintf('    OK  : the %u ranked features are in MeanAUC order (ties by MinAUC, ID)\n',height(r));
    else
        fprintf(2,'    FAIL: the ranking does not follow MeanAUC\n'); nFail = nFail + 1;
    end
end
if R.options.DoPlots && contains(fileread(R.files.log),'ranking figures not drawn')
    fprintf(2,'    FAIL: the ranking figures were not drawn (see the log)\n'); nFail = nFail + 1;
end

% Redundancy: the top RedundancyTopN are checked, the rest are not
F = R.FeatureTable;
nChecked = min(R.options.RedundancyTopN,sum(isfinite(F.Rank)));
isTop = isfinite(F.Rank) & F.Rank <= nChecked;
okDistinct = all(ismember(F.Distinct(isTop),["yes","no"])) && all(F.Distinct(~isTop) == "not checked") ...
    && isequal(F.Distinct == "no",F.RedundantWith ~= "");
[~,jA] = ismember(A.ID,F.ID);
if okDistinct && isequal(A.Distinct,F.Distinct(jA)) && isequaln(A.RedundancyRho,F.RedundancyRho(jA))
    fprintf('    OK  : %u features checked for redundancy: %u distinct, %u redundant\n', ...
        nChecked,sum(F.Distinct == "yes"),sum(F.Distinct == "no"));
else
    fprintf(2,'    FAIL: Distinct/RedundantWith inconsistent, or FeatureAUC''s differ\n'); nFail = nFail + 1;
end
end

%% ---------------------------------------------------------------------------
function q = bh(p)
% Benjamini-Hochberg by its definition: q_i = min over p_j >= p_i of m*p_j/#{p <= p_j}.
q = nan(size(p));
ok = find(~isnan(p)); m = numel(ok);
for i = ok(:)'
    pj = p(ok(p(ok) >= p(i)));
    q(i) = min(1,min(m*pj./arrayfun(@(t) sum(p(ok) <= t),pj)));
end
end

%% ---------------------------------------------------------------------------
function out = ternary(c,a,b)
if c, out = a; else, out = b; end
end

%% ---------------------------------------------------------------------------
function nFail = compareRuns(R1,R2,what)
% Two runs on identical numbers; the features must match, NaN for NaN.
fprintf('\n  comparing %s\n',what);
nFail = 0;
try
    v1 = table2array(R1.SeriesTable(:,varfun(@isnumeric,R1.SeriesTable,'OutputFormat','uniform')));
    v2 = table2array(R2.SeriesTable(:,varfun(@isnumeric,R2.SeriesTable,'OutputFormat','uniform')));
    if ~isequal(size(v1),size(v2))
        fprintf(2,'    FAIL: different shapes, %s vs %s\n',mat2str(size(v1)),mat2str(size(v2)));
        nFail = 1; return
    end
    nanDiff = sum(isnan(v1(:)) ~= isnan(v2(:)));
    both = isfinite(v1) & isfinite(v2);
    d = max(abs(v1(both) - v2(both)));
    if nanDiff > 0
        fprintf(2,'    FAIL: %u values are NaN in one run only\n',nanDiff); nFail = 1;
    elseif ~any(both(:))
        fprintf(2,'    FAIL: no values to compare\n'); nFail = 1;
    elseif d == 0
        fprintf('    OK  : identical feature values (%u compared)\n',sum(both(:)));
    elseif d < 1e-10
        fprintf('    OK  : agree to %.1e (rounding only)\n',d);
    else
        fprintf(2,'    FAIL: they disagree by up to %.3e -- one is mangling the data\n',d);
        nFail = 1;
    end
catch err
    fprintf(2,'    FAIL: could not compare: %s\n',err.message); nFail = 1;
end
end

%% ---------------------------------------------------------------------------
function nFail = checkReduced(R,label,thr,ranked)
% 'ReduceFeatureSet': HCTSA_N_reduced.mat and the cluster columns, checked against
% HCTSA_N.mat and an independent complete-linkage clustering.
fprintf('\n  [%s] checking the reduced feature set\n',label);
nFail = 0;
try
    if ~isfield(R.files,'reducedFile') || ~isfile(R.files.reducedFile) || ~isfield(R,'reducedFeatureSet')
        nFail = report9(false,'','no HCTSA_N_reduced.mat, or no results.reducedFeatureSet'); return
    end
    red = R.reducedFeatureSet; F = R.FeatureTable;
    N = load(R.files.normFile); M = load(R.files.reducedFile);
    [~,~,OpsR] = TS_LoadData(R.files.reducedFile);
    [~,jN] = ismember(M.Operations.ID,N.Operations.ID);
    nFail = nFail + report9(isequal(sort(OpsR.ID),sort(red.ID(red.Kept))) && isequal(sort(red.ID),sort(N.Operations.ID)), ...
        sprintf('HCTSA_N_reduced.mat loads; it holds the %u kept features of %u',sum(red.Kept),height(red)), ...
        'the reduced file''s features are not the kept ones');
    nFail = nFail + report9(all(jN > 0) && isequaln(M.TS_DataMat,N.TS_DataMat(:,jN)) && ...
        isequaln(M.TS_Quality,N.TS_Quality(:,jN)) && isequaln(M.TimeSeries,N.TimeSeries) && ...
        all(ismember(fieldnames(N),fieldnames(M))) && isequaln(M.reducedFeatureSet,red), ...
        'its values, series and other variables are HCTSA_N.mat''s, plus reducedFeatureSet', ...
        'the reduced file differs from HCTSA_N.mat beyond dropping features');

    % the clusters, independently: |Spearman rho| pair by pair on the series both
    % features have (ranks, then Pearson by its formula; undefined = 0), complete linkage
    X = N.TS_DataMat; X(N.TS_Quality > 0) = NaN;
    nOp = size(X,2); rho = eye(nOp);
    for i = 1:nOp
        for j = i+1:nOp
            ok = ~isnan(X(:,i)) & ~isnan(X(:,j));
            a = tiedrank(X(ok,i)); b = tiedrank(X(ok,j));
            a = a - mean(a); b = b - mean(b);
            r = sum(a.*b)/sqrt(sum(a.^2)*sum(b.^2));
            if ~isfinite(r), r = 0; end
            rho(i,j) = abs(r); rho(j,i) = abs(r);
        end
    end
    nOpNaN = sum(any(isnan(X),1));
    nFail = nFail + report9(nOpNaN >= 2, ...
        sprintf('HCTSA_N.mat has missing values in %u features (pairwise correlations)',nOpNaN), ...
        'HCTSA_N.mat has no missing values: the pairwise path is not exercised');
    D = 1 - rho; D(1:nOp+1:end) = 0;
    T = cluster(linkage(squareform(D),'complete'),'cutoff',1 - thr,'criterion','distance');
    [~,jR] = ismember(N.Operations.ID,red.ID);
    C = red.Cluster(jR);                                   % in HCTSA_N's order
    samePartition = max(T) == max(C) && all(arrayfun(@(c) numel(unique(T(C == c))) == 1 && ...
        sum(T == T(find(C == c,1))) == sum(C == c), 1:max(C)));
    within = arrayfun(@(c) min(min(rho(C == c,C == c))), 1:max(C));
    nFail = nFail + report9(samePartition && all(within >= thr - 1e-12), ...
        sprintf('%u clusters, as complete linkage gives; every within-cluster |rho| >= %g (lowest %.3f)', ...
            max(C),thr,min(within)), ...
        'the clusters differ from an independent complete linkage, or a pair is below the threshold');
    nMulti = sum(accumarray(C,1) > 1);
    nFail = nFail + report9(nMulti > 0,sprintf('%u clusters have several members',nMulti), ...
        'every cluster is a single feature: the test data do not exercise the choice');

    % one kept feature per cluster, the right one, and consistent columns
    [~,jF] = ismember(red.ID,F.ID);
    one = accumarray(red.Cluster,double(red.Kept)) == 1;
    sizes = accumarray(red.Cluster,1);
    if ranked
        rnk = F.Rank(jF);
        right = true;
        for c = 1:max(red.Cluster)
            m = red.Cluster == c;
            if any(isfinite(rnk(m)))
                right = right && isequal(rnk(m & red.Kept),min(rnk(m)));
            else
                right = right && isequal(red.Central(m & red.Kept),true);
            end
        end
        rule = 'the best-ranked member (FeatureTable.Rank)';
    else
        right = isequal(red.Kept,red.Central) && all(isnan(red.Rank));
        rule = 'the most central member';
    end
    repOK = false;
    if all(one)
        keptName = red.Name(red.Kept); keptCl = red.Cluster(red.Kept);
        [~,k] = ismember(red.Cluster,keptCl);
        repOK = isequal(red.Representative,keptName(k)) && isequal(red.RepresentativeID(red.Kept),red.ID(red.Kept));
    end
    nFail = nFail + report9(all(one) && right && repOK && isequal(red.ClusterSize,sizes(red.Cluster)), ...
        sprintf('one feature kept per cluster: %s',rule),'the kept features or the cluster table are wrong');
    inN = ismember(F.ID,red.ID);
    nFail = nFail + report9(any(~inN) && all(isnan(F.Cluster(~inN))) && all(F.Representative(~inN) == ""), ...
        sprintf('%u feature(s) dropped by normalisation: no cluster, not in reducedFeatureSet',sum(~inN)), ...
        'no feature was dropped by normalisation, or a dropped one has a cluster');
    wantRep = red.Representative; wantRep(red.Kept) = "";
    nFail = nFail + report9(all(jF > 0) && isequal(F.Cluster(jF),red.Cluster) && isequal(F.ClusterSize(jF),red.ClusterSize) && ...
        isequal(F.Representative(jF),wantRep) && all(isnan(F.Cluster(~inN))), ...
        'FeatureTable''s Cluster, ClusterSize and Representative match ("" for the kept feature)', ...
        'FeatureTable''s cluster columns do not match reducedFeatureSet');
    if height(R.FeatureAUC) > 0
        A = R.FeatureAUC; [~,jA] = ismember(A.ID,F.ID);
        nFail = nFail + report9(isequaln(A.Cluster,F.Cluster(jA)) && isequal(A.Representative,F.Representative(jA)), ...
            'FeatureAUC carries the same cluster columns','FeatureAUC''s cluster columns differ');
    end

    % TS_ReduceFeatureSet itself: its centres, and the figure (drawn invisibly)
    [ids0,clu0] = TS_ReduceFeatureSet(R.files.normFile,1 - thr,true,false);
    prevVis = get(groot,'DefaultFigureVisible'); set(groot,'DefaultFigureVisible','off');
    figsBefore = findobj('Type','figure');
    try
        [ids1,clu1] = TS_ReduceFeatureSet(R.files.normFile,1 - thr,true,true);
    catch err
        set(groot,'DefaultFigureVisible',prevVis); rethrow(err);
    end
    set(groot,'DefaultFigureVisible',prevVis);
    newFigs = setdiff(findobj('Type','figure'),figsBefore);
    nAxes = numel(findobj(newFigs,'Type','axes'));
    close(newFigs);
    nFail = nFail + report9(isequal(sort(ids0),sort(red.ID(red.Central))) && isequal(ids0,ids1) && isequal(clu0,clu1) && ...
        numel(newFigs) == 1 && nAxes >= 2, ...
        'TS_ReduceFeatureSet: centres as recorded; with its figure (dendrogram and matrix) the same clusters', ...
        'TS_ReduceFeatureSet''s figure run differs, or did not draw one figure');
    S1 = struct('TS_DataMat',N.TS_DataMat(:,1),'TS_Quality',N.TS_Quality(:,1), ...
        'TimeSeries',N.TimeSeries,'Operations',N.Operations(1,:));
    [ids1f,clu1f] = TS_ReduceFeatureSet(S1,1 - thr,true,false);
    nFail = nFail + report9(isequal(ids1f,N.Operations.ID(1)) && height(clu1f) == 1 && ...
        clu1f.Cluster == 1 && clu1f.ClusterSize == 1 && clu1f.Central, ...
        'TS_ReduceFeatureSet: a single feature is one cluster','TS_ReduceFeatureSet fails on a single feature');
catch err
    nFail = nFail + report9(false,'',err.message);
end
end

function n = report9(ok,okMsg,failMsg)
% Print a stage 9 check; 1 if it failed.
if ok
    fprintf('    OK  : %s\n',okMsg);
else
    fprintf(2,'    FAIL: %s\n',failMsg);
end
n = double(~ok);
end

%% ---------------------------------------------------------------------------
function nFail = compareReduceOff(Ron,Roff)
% The option adds the cluster columns and the file, and changes nothing else.
fprintf('\n  [stage 9] the same run without ReduceFeatureSet\n');
nFail = 0;
newCols = {'Cluster','ClusterSize','Representative'};
try
    same = isequaln(removevars(Ron.FeatureTable,newCols),Roff.FeatureTable) && ...
        isequaln(removevars(Ron.FeatureAUC,newCols),Roff.FeatureAUC) && ...
        isequaln(Ron.SeriesTable,Roff.SeriesTable) && ...
        ~any(ismember(newCols,Roff.FeatureTable.Properties.VariableNames)) && ...
        ~isfield(Roff,'reducedFeatureSet') && ~isfield(Roff.files,'reducedFile') && ...
        ~isfile(fullfile('testrun_reduce_off','HCTSA_N_reduced.mat'));
    if same
        fprintf('    OK  : without it, the tables are the same less the three cluster columns, and no file\n');
    else
        fprintf(2,'    FAIL: ReduceFeatureSet changes more than its columns and file\n'); nFail = 1;
    end
catch err
    fprintf(2,'    FAIL: could not compare: %s\n',err.message); nFail = 1;
end
end

%% ---------------------------------------------------------------------------
function nFail = checkComputeFeatures(R4,rs)
% HCTSA2_ComputeFeatures against stage 4's HCTSA.mat (catch22), direct calls and
% known values (the hctsa library).
nFail = 0;
try
    H = load(fullfile('testrun_twoclass','HCTSA.mat'));
    X = H.TS_DataMat; X(H.TS_Quality > 0) = NaN;
    ops = H.Operations;
    req = ops.ID(end:-1:1);                                  % every feature, reversed
    [~,j] = ismember(req,ops.ID);
    F = R4.FeatureTable; [~,jF] = ismember(req,F.ID);
    okID = true; okName = true; okField = true; okPad = true; okCols = true;
    for s = [1 25]                                           % one series of each class
        x = H.TimeSeries.Data{s};
        [Fs,info] = HCTSA2_ComputeFeatures(x,req,'FeatureSet','catch22');
        okID = okID && isequaln(cell2mat(struct2cell(Fs)),X(s,j)') && isequal(info.ID,req) && ...
            isequal(fieldnames(Fs),cellstr(info.Field)) && all(info.Quality(isfinite(X(s,j))) == 0);
        okName = okName && isequaln(HCTSA2_ComputeFeatures(x,ops.Name(j),'FeatureSet','catch22'),Fs);
        okField = okField && isequaln(HCTSA2_ComputeFeatures(x,cellstr(info.Field),'FeatureSet','catch22'),Fs);
        w0 = warning;
        okPad = okPad && isequaln(HCTSA2_ComputeFeatures([NaN; NaN; x(:); Inf],req,'FeatureSet','catch22'),Fs) && ...
            isequal(warning,w0) && ...
            isequaln(HCTSA2_ComputeFeatures(x(:)',req,'FeatureSet','catch22'),Fs);
        hasCol = F.Column(jF) ~= "";
        okCols = okCols && isequal(info.Field(hasCol),F.Column(jF(hasCol)));
    end
    nFail = nFail + report9(okID, ...
        'by ID: stage 4''s catch22 values for two series, bit for bit, in the requested order', ...
        'by ID: the values differ from HCTSA.mat, or the order or the info table is wrong');
    nFail = nFail + report9(okName && okField,'by name and by field name: the same', ...
        'by name or by field name: not the same as by ID');
    nFail = nFail + report9(okPad,'NaN/Inf padding at the ends trimmed; a row vector works; warning settings left as they were', ...
        'padding or a row vector changes the values, or the warning settings changed');
    nFail = nFail + report9(okCols,'the fields are FeatureTable.Column', ...
        'the fields differ from FeatureTable.Column');

    % NaN inside: every feature NaN, a warning, Quality NaN
    x2 = H.TimeSeries.Data{1}; x2(50) = NaN;
    lastwarn('');
    [F2,info2] = HCTSA2_ComputeFeatures(x2,req,'FeatureSet','catch22');
    [~,wid] = lastwarn;
    nFail = nFail + report9(all(isnan(cell2mat(struct2cell(F2)))) && all(isnan(info2.Quality)) && ...
        strcmp(wid,'HCTSA2_ComputeFeatures:nonFinite'), ...
        'NaN inside the series: every feature NaN, with a warning','NaN inside the series: not handled');
    threw = '';
    try
        HCTSA2_ComputeFeatures(x2,{'no_such_feature'},'FeatureSet','catch22');
    catch err
        threw = err.identifier;
    end
    nFail = nFail + report9(strcmp(threw,'HCTSA2_ComputeFeatures:unknown'), ...
        'an unknown feature is an error','an unknown feature: no error');

    % The full library: field names (no computation: a series with NaN inside)
    w = warning('off','HCTSA2_ComputeFeatures:nonFinite');
    [Fh,infoH] = HCTSA2_ComputeFeatures([1; NaN; 2],1:6392);
    warning(w);
    fh = fieldnames(Fh);
    k = find(infoH.Name == "firstMax_mi-gaussian");
    nFail = nFail + report9(numel(fh) == 6392 && all(cellfun(@isvarname,fh)) && ...
        numel(unique(fh)) == 6392 && sum(infoH.Field ~= infoH.Name) == 72 && ...
        infoH.Field(k) == "firstMax_mi_gaussian", ...
        'hctsa library: 6,392 distinct valid field names, 72 of them differing from the name', ...
        'hctsa library: the field names are not as documented');

    % (3.0.0) 'NormalizeLike': stage 4's own series, normalised as stage 4 normalised
    % them, are its HCTSA_N.mat values (TS_Normalize applied the same stored transform)
    fN = fullfile('testrun_twoclass','HCTSA_N.mat');
    N = load(fN);
    okN = true; okRaw = true;
    for s = [1 25]
        rowN = find(N.TimeSeries.ID == H.TimeSeries.ID(s));
        x = H.TimeSeries.Data{s};
        [Fn,infoN] = HCTSA2_ComputeFeatures(x,N.Operations.ID,'FeatureSet','catch22', ...
            'NormalizeLike',fN);
        okN = okN && isequaln(cell2mat(struct2cell(Fn)),N.TS_DataMat(rowN,:)') && ...
            isequaln(HCTSA2_ComputeFeatures(x,N.Operations.ID,'FeatureSet','catch22', ...
            'NormalizeLike','testrun_twoclass'),Fn);
        [~,jr] = ismember(N.Operations.ID,ops.ID);
        okRaw = okRaw && isequaln(infoN.RawValue,X(s,jr)');
    end
    nFail = nFail + report9(okN && okRaw, ...
        'NormalizeLike: two of the run''s series get their HCTSA_N.mat values, bit for bit (file or folder); RawValue is the raw value', ...
        'NormalizeLike: the values differ from HCTSA_N.mat, or RawValue is wrong');
    % a run without the stored normalisation (before 3.0.0), and a filtered-out feature
    old = N.normalizationInfo; old = rmfield(old,'transform');
    threw = '';
    try
        HCTSA2_ComputeFeatures(x,N.Operations.ID(1),'FeatureSet','catch22','NormalizeLike',old);
    catch err
        threw = err.identifier;
    end
    okOld = strcmp(threw,'HCTSA2_ComputeFeatures:noTransform');
    gone = setdiff(ops.ID,N.Operations.ID);
    if isempty(gone)
        % the run's filter kept every feature: remove one from a copy of its transform
        cut = N.normalizationInfo;
        goneID = cut.transform.featureParams.ID(1);
        cut.transform.featureParams(1,:) = [];
        src = cut;
    else
        goneID = gone(1);
        src = fN;
    end
    threw = '';
    try
        HCTSA2_ComputeFeatures(x,goneID,'FeatureSet','catch22','NormalizeLike',src);
    catch err
        threw = err.identifier;
    end
    okGone = strcmp(threw,'HCTSA2_ComputeFeatures:notNormalized');
    nFail = nFail + report9(okOld && okGone, ...
        sprintf(['NormalizeLike: a run without the stored normalisation is an error; ' ...
        'a feature the run filtered out is an error (%u filtered out by the run%s)'], ...
        numel(gone),ternary(isempty(gone),'; one removed from a copy','')), ...
        'NormalizeLike: an old run or a filtered-out feature is not refused');

    % Known values, and two outputs of one master operation against a direct call
    x3 = 10 + 3*randn(rs,300,1);
    [F3,info3] = HCTSA2_ComputeFeatures(x3,{'mean','median','standard_deviation', ...
        'SB_MotifTwo_mean_dd','SB_MotifTwo_mean_u'});
    out = SB_MotifTwo(zscore(x3),'mean');
    nFail = nFail + report9(abs(F3.mean - mean(x3)) < 1e-12 && F3.median == median(x3) && ...
        abs(F3.standard_deviation - std(x3)) < 1e-12 && F3.SB_MotifTwo_mean_dd == out.dd && ...
        F3.SB_MotifTwo_mean_u == out.u && info3.MasterID(4) == info3.MasterID(5) && ...
        isequaln(HCTSA2_ComputeFeatures(x3,info3.ID),F3), ...
        'hctsa library: mean, median, standard deviation and SB_MotifTwo(x_z,''mean'') outputs as computed directly', ...
        'hctsa library: values differ from a direct computation');
catch err
    nFail = nFail + report9(false,'',err.message);
end
end

%% ---------------------------------------------------------------------------
function X = arSeries(rs,phi,n,L)
% n AR(1) series of length L, one per row, from the random stream rs.
X = zeros(n,L);
e = randn(rs,n,L);
for t = 2:L
    X(:,t) = phi*X(:,t-1) + e(:,t);
end
end
