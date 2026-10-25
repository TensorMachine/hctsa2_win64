function report = HCTSA2_CheckInstall(verbose)
% HCTSA2_CheckInstall   Report exactly what this installation can and cannot do.
%
%   HCTSA2_CheckInstall
%   report = HCTSA2_CheckInstall(false)
%
% Run after `startup`; the folder you run it from does not matter (HCTSA2 2.2.1:
% every path is taken from this file's own location). Checks MATLAB version, required
% toolboxes, every compiled component and the TISEAN executables, then reports how
% many of the library features (6,392 in 3.0.0) are available.
%
% Nothing here modifies the installation. It answers one question: if a feature
% comes back NaN, is that a real result, a missing binary, or a broken install?

if nargin < 1, verbose = true; end

report = struct();
root = fileparts(mfilename('fullpath'));
pass = @(b) tern(b,'   OK',' FAIL');

if verbose
    fprintf('\n========== hctsa2 installation check ==========\n');
    fprintf('MATLAB %s on %s\n\n',version,computer);
end

%-------------------------------------------------------------------------------
%% 1. MATLAB version and required toolboxes
%-------------------------------------------------------------------------------
v = sscanf(version('-release'),'%d');
report.matlabRelease = version('-release');
report.matlabOK = ~isempty(v) && v >= 2024;     % README/INSTALL: R2024a or later

report.statsToolbox  = license('test','Statistics_Toolbox')  == 1;
report.signalToolbox = license('test','Signal_Toolbox')      == 1;

if verbose
    fprintf('-- MATLAB and required toolboxes --\n');
    fprintf('  %s  MATLAB release %s (R2024a or later)\n',pass(report.matlabOK),report.matlabRelease);
    fprintf('  %s  Statistics and Machine Learning Toolbox\n',pass(report.statsToolbox));
    fprintf('  %s  Signal Processing Toolbox\n',pass(report.signalToolbox));
    fprintf('\n-- toolboxes hctsa2 no longer needs --\n');
    for t = {'Curve_Fitting_Toolbox','Wavelet_Toolbox','Econometrics_Toolbox', ...
             'Identification_Toolbox','Financial_Toolbox'}
        have = license('test',t{1}) == 1;
        fprintf('     %-26s %s\n',t{1},tern(have,'present (unused)','absent (fine)'));
    end
end

%-------------------------------------------------------------------------------
%% 2. Compiled components
%-------------------------------------------------------------------------------
% {label, probe functions (ALL must be built), features, hasMatlabFallback}
%
% Detecting a built MEX needs care, and both obvious tests are wrong.
%
%   exist(fn,'file') == 2  -> false OK.   Several components ship a same-named
%                             .m beside the C source. MS_complexitybs.m and
%                             nn_prepare.m are documentation-only stubs with no
%                             executable code, yet they satisfy this test.
%   exist(fn,'file') == 3  -> false FAIL. When a .m sits alongside the MEX,
%                             exist reports 2 for the .m and the built MEX is
%                             declared missing -- which is what happened to
%                             Michael_Small, OpenTSTOOL and gpml.
%
% mexIsBuilt() below asks the unambiguous question instead: is there a file
% named <fn>.<mexext> anywhere the code can see it?
%
% gpml remains a special case: solve_chol.m is a real MATLAB implementation, so
% its features work without the MEX, just more slowly.
%
% (HCTSA2 2.2.1: one entry per BUILD TARGET, as HCTSA2_BuildBinaries has them.
% Max_Little used to be probed through ML_fastdfa_core alone and credited with
% all 39 features, though they come from three targets -- fastdfa 1, rpde 10 and
% steps_bumps 28 -- and steps_bumps is C++, which needs a different compiler.
% 2.5.0: upstream f54e4d1 dropped SC_FastDFA and the fastdfa code.)
%
% Feature counts are for FeatureSets/INP_ops_hctsa.txt at upstream v3.0.0, 4581e712
% (hctsa2 3.0.0): the features of every operation that calls one of the
% component's MEX functions, directly, through the component's own .m files, or
% through another operation (EN_MSE -> EN_SampEn). Counted by
% docs/hctsa2_2.2.0/tools/component_features.py. No hctsa feature uses the
% catch22 MEX files; only the separate catch22/catch24 feature sets do, so they
% are checked but gate nothing here. (2.5.0's counts: rpde 10, Michael_Small 53,
% sampen 112, KSG 453, gpml 147; upstream v3.0.0's redundancy review removed
% features from each.)
C = {
 'Max_Little rpde',                 {'ML_close_ret'},               8, false
 'Max_Little steps_bumps (C++)',    {'ML_kvsteps_core'},           28, false
 'Michael_Small',                   {'MS_complexitybs','MS_nearest','MS_shannon'}, 47, false
 'Physionet sampen',                {'sampen_mex'},                91, false
 'KSG mutual information',          {'KSG_mi'},                   440, false
 'catch22 (catch22 feature sets)',  {'catch22_DN_HistogramMode_5'}, 0, false
 ... % OpenTSTOOL was removed with upstream v2.0.0's TSTOOL replacement. Its 619
 ... % features are now produced by native and TISEAN-based operations, so there
 ... % is no MEX component to check for.
 'gpml (MEX = speed only)',         {'solve_chol'},                30, true
};
report.mex = struct('name',{},'available',{},'features',{},'fallback',{},'missing',{});
% (2.4.0) Some KSG features also need TISEAN: a 'mi' time delay with a TISEAN program
% or ripser or with an 'fnn' embedding, 104 in all (2.5.0 and 3.0.0; 150 in 2.4.0, before
% NL_PoincareSection stopped calling TISEAN). They are in both counts, so the totals
% below subtract them once. The KSG features are those of the rule "'kraskov' or a 'mi'
% delay, plus NL_PersistentHomology", which the dependency runs (KSG_mi disabled)
% confirmed in 2.4.0 and 2.5.0.
report.ksgTiseanOverlap = 104;
mexOK = 0; mexBad = 0;
if verbose, fprintf('\n-- compiled components (MEX) --\n'); end
for i = 1:size(C,1)
    nm = C{i,1}; probes = C{i,2}; nf = C{i,3}; fb = C{i,4};
    built = cellfun(@(f) mexIsBuilt(f,root),probes);
    have = all(built);
    % The fallback entry (gpml) is counted with the pure-MATLAB features below,
    % so neither of these sums includes it.
    if ~fb
        if have, mexOK = mexOK + nf; else, mexBad = mexBad + nf; end
    end
    report.mex(end+1) = struct('name',nm,'available',have,'features',nf, ...
        'fallback',fb,'missing',strjoin(probes(~built),', ')); %#ok<AGROW>
    if verbose
        if have
            fprintf('  %s  %-32s %4d features\n',pass(true),nm,nf);
        elseif fb
            fprintf('  %s  %-32s %4d features (MATLAB fallback in use)\n','  ---',nm,nf);
        else
            fprintf('  %s  %-32s %4d features (missing: %s)\n',pass(false),nm,nf, ...
                strjoin(probes(~built),', '));
        end
    end
end

%-------------------------------------------------------------------------------
%% 3. TISEAN -- MEX files or executables
%-------------------------------------------------------------------------------
% Eleven TISEAN programs and ripser, all supplied by the executable build (2.5.0:
% poincare is no longer called; NL_PoincareSection computes in MATLAB).
tiseanProgs = {'d2','nstat_z','false_nearest','boxcount','lyap_r', ...
               'lyap_spec','nrlazy','c1','c2d','c2g','c2t','ripser'};
nProgs = numel(tiseanProgs);
% Where startup.m and BF_RunTisean look: Toolboxes/Tisean_bin on Windows (ripser
% included); elsewhere Toolboxes/Tisean_3.0.1/bin if it exists, else ~/bin, and
% ripser in Toolboxes/ripser/bin (upstream's compile_ripser) or the same folder.
if ispc
    tiseanDir = fullfile(root,'Toolboxes','Tisean_bin');
    exeExtn = '.exe';
    ripserDirs = {tiseanDir};
else
    tiseanDir = fullfile(root,'Toolboxes','Tisean_3.0.1','bin');
    if exist(tiseanDir,'dir') ~= 7
        tiseanDir = fullfile(getenv('HOME'),'bin');
    end
    exeExtn = '';
    ripserDirs = {fullfile(root,'Toolboxes','ripser','bin'), tiseanDir};
end
haveExe = false(1,nProgs);
for i = 1:nProgs
    if strcmp(tiseanProgs{i},'ripser')
        haveExe(i) = any(cellfun(@(d) isfile(fullfile(d,['ripser' exeExtn])),ripserDirs));
    else
        haveExe(i) = isfile(fullfile(tiseanDir,[tiseanProgs{i} exeExtn]));
    end
end
% (HCTSA2 3.0.0) A false_nearest built before 3.0.0 writes four columns, and NL_FNN
% now refuses it (it needs the fifth, the number of points behind each fraction), so
% every 'fnn' feature fails until TISEAN is rebuilt. It is run once on a short series
% to tell, and an old build counts as missing.
report.falseNearestOld = false;
kFnn = find(strcmp(tiseanProgs,'false_nearest'));
if haveExe(kFnn)
    report.falseNearestOld = localOldFalseNearest();
    haveExe(kFnn) = ~report.falseNearestOld;
end
report.tiseanDir = tiseanDir;
report.tiseanExeCount = sum(haveExe);
report.tisean = all(haveExe);

% 538 feature values come from the eleven operations that call these programs
% directly: NL_FNN 114, NL_d2 100, NL_Dimensions 87, NL_BoxCountEntropyRate 63,
% NL_LargestLyap 60, SY_nstat_z 42, NL_GPCorrSum 23, NL_c1 20,
% PP_SchreiberDenoise 14, NL_PersistentHomology 8 (ripser), NL_LyapSpec 7. A further
% 50 reach false_nearest only through an 'fnn' embedding dimension (NL_nlpe 15,
% NL_LocalDensity 14, NL_FractalDimensions 8, NL_RQA 8, NL_RecurrenceTimes 4,
% NL_TakensEstimator 1). (Counted from FeatureSets/INP_ops_hctsa.txt at upstream
% v3.0.0 by docs/hctsa2_2.2.0/tools/component_features.py; 2.5.0's: 604 and 51.)
%
% (HCTSA2 2.2.1: gated PER PROGRAM. The table says which programs each group of
% features needs, so a build without gfortran loses the 124 features of NL_c1
% and NL_d2 and not all 793, and a missing ripser costs its 8. The 51 'fnn'
% features used to be counted as pure MATLAB, which overstated what was
% available when false_nearest was missing. The program lists come from the
% BF_RunTisean calls in each operation; NL_Dimensions, NL_GPCorrSum and the
% second NL_LargestLyap entry also embed with 'fnn' in every registration.)
fnnLabel = ['NL_nlpe, NL_LocalDensity, NL_FractalDimensions, NL_RQA, ' ...
            'NL_RecurrenceTimes, NL_TakensEstimator (''fnn'' embedding only)'];
% {operation(s), feature values, programs needed, calls a program directly}
tiseanOps = {
 'NL_Dimensions',            87, {'boxcount','d2','false_nearest'},    true
 'NL_FNN',                  114, {'false_nearest'},                    true
 'NL_d2',                   100, {'d2','c2g','c2t'},                   true
 'NL_BoxCountEntropyRate',   63, {'boxcount'},                         true
 'NL_LargestLyap',           30, {'lyap_r'},                           true
 'NL_LargestLyap',           30, {'lyap_r','false_nearest'},           true
 'SY_nstat_z',               42, {'nstat_z'},                          true
 'NL_GPCorrSum',             23, {'d2','false_nearest'},               true
 'NL_c1',                    20, {'c1','c2d'},                         true
 'PP_SchreiberDenoise',      14, {'nrlazy'},                           true
 'NL_PersistentHomology',     8, {'ripser'},                           true
 'NL_LyapSpec',               7, {'lyap_spec'},                        true
 fnnLabel,                   50, {'false_nearest'},                    false
};
opFeatures = [tiseanOps{:,2}]';
opDirect = [tiseanOps{:,4}]';
opLost = cellfun(@(p) any(ismember(p,tiseanProgs(~haveExe))),tiseanOps(:,3));
report.tiseanFeatures = sum(opFeatures(opDirect));       % 538: the programs' own operations
report.tiseanFnnFeatures = sum(opFeatures(~opDirect));   % 50: through an 'fnn' embedding only
report.tiseanMissing = tiseanProgs(~haveExe);
report.tiseanUnavailable = sum(opFeatures(opLost));      % feature values lost to missing programs
tiseanTotal = report.tiseanFeatures + report.tiseanFnnFeatures;
if verbose
    fprintf('\n-- TISEAN --\n');
    fprintf('  %s  %-30s %2u of %u programs (TISEAN + ripser)\n', ...
            pass(report.tisean),tiseanDir,report.tiseanExeCount,nProgs);
    fprintf('       %4d feature values need them (%d directly, %d more through an ''fnn'' embedding)\n', ...
            tiseanTotal,report.tiseanFeatures,report.tiseanFnnFeatures);
    if report.falseNearestOld
        fprintf(2,['  FAIL  false_nearest is a build from before HCTSA2 3.0.0 (four output ' ...
            'columns): every ''fnn'' feature fails\n        until TISEAN is rebuilt ' ...
            '(HCTSA2_BuildTisean, then startup)\n']);
    end
    if ~report.tisean
        fprintf('        missing: %s\n',strjoin(report.tiseanMissing,', '));
        fprintf('        costs %d feature values: %s\n',report.tiseanUnavailable, ...
                strjoin(unique(tiseanOps(opLost,1),'stable')',', '));
        fprintf('        (build them with HCTSA2_BuildTisean, then run startup)\n');
    end
end

% (HCTSA2 2.4.0: no Java library to check. The JIDT jar is gone; its Kraskov
% estimator is the KSG MEX above.)

%-------------------------------------------------------------------------------
%% 4. Live smoke test on the pure-MATLAB core
%-------------------------------------------------------------------------------
report.smokeOK = false;
try
    rng(1);
    y = cumsum(randn(500,1));
    y = (y-mean(y))/std(y);
    a = CO_AutoCorr(y,1,'Fourier');
    b = DN_Mean(y);
    report.smokeOK = isfinite(a) && isfinite(b);
catch
end
if verbose
    fprintf('\n-- smoke test --\n');
    fprintf('  %s  pure-MATLAB core executes\n',pass(report.smokeOK));
end

%-------------------------------------------------------------------------------
%% Summary
%-------------------------------------------------------------------------------
% Read the library size from the feature set rather than hard-coding it. The
% previous constants (7702 total, 6403 pure MATLAB) were measured before the
% upstream v2.0.0 merge and became wrong the moment TSTOOL was removed -- the
% checker then reported "6475 of 7702" for a library that has 7693 features.
% A count that has to be remembered will eventually be wrong.
featuresTotal = countFeatureSet(root);
gated = 0;
for i = 1:numel(report.mex)
    if ~report.mex(i).available && ~report.mex(i).fallback
        gated = gated + report.mex(i).features;
    end
end
gated = gated + report.tiseanUnavailable;
ksgMissing = any(strncmp({report.mex.name},'KSG',3) & ~[report.mex.available]);
if ksgMissing && ~report.tisean          % the shared features, lost once not twice
    gated = gated - report.ksgTiseanOverlap;
end
avail = max(featuresTotal - gated,0);
report.featuresAvailable = avail;
report.featuresTotal = featuresTotal;
pureMatlab = featuresTotal - totalCompiledFeatures(report);

if verbose
    fprintf('\n=================================================\n');
    fprintf('  features available : %4d of %u  (%.1f%%)\n',avail,featuresTotal, ...
        100*avail/max(featuresTotal,1));
    fprintf('    pure MATLAB      : %4d  (always)\n',pureMatlab);
    fprintf('    from MEX         : %4d\n',mexOK);
    fprintf('    from TISEAN/ripser: %3d\n',tiseanTotal - report.tiseanUnavailable);
    fprintf('    (%u of them need both the KSG MEX and TISEAN, so are in both lines)\n', ...
        report.ksgTiseanOverlap);
    if mexBad > 0 || report.tiseanUnavailable > 0
        fprintf('  unavailable        : %4d  -- run HCTSA2_BuildBinaries (MEX) or HCTSA2_BuildTisean, or see INSTALL.md\n', ...
            mexBad + report.tiseanUnavailable);
    end
    fprintf('=================================================\n');
    if ~report.statsToolbox || ~report.signalToolbox
        fprintf(2,'\nA REQUIRED toolbox is missing. Most of the library will not run.\n');
    elseif avail == featuresTotal
        fprintf('\nComplete installation -- every library feature is available.\n');
    end
    fprintf('\n');
end
end

%-------------------------------------------------------------------------------
function s = tern(c,a,b)
if c, s = a; else, s = b; end
end

%-------------------------------------------------------------------------------
function old = localOldFalseNearest()
% True if false_nearest writes four columns (a build from before HCTSA2 3.0.0)
% rather than five. Runs it once on a short series through BF_RunTisean, as NL_FNN
% does; any failure to run it is left to the other checks (not counted as old).
old = false;
f = '';
try
    y = sin((1:300)'/5) + 0.1*cos((1:300)'/1.3);
    f = BF_WriteTempFile(y);
    [~,res] = BF_RunTisean(sprintf('false_nearest -d1 -m1 -M1,2 -t1 -V0 %s',f));
    line = strtrim(strtok(res,newline));
    if ~isempty(line) && ~isnan(str2double(strtok(line)))
        old = numel(strsplit(line)) < 5;
    end
catch
end
if ~isempty(f)
    try, BF_CleanTempFiles(f); catch, end
end
end

%-------------------------------------------------------------------------------
function tf = mexIsBuilt(fn,root)
% Is a compiled MEX for FN present?
%
% Both obvious tests are wrong when a same-named .m sits beside the MEX:
%   exist(fn,'file') == 2  -> false OK   (a .m stub passes)
%   exist(fn,'file') == 3  -> false FAIL (the .m is reported and the MEX missed)
% The second is what made Michael_Small, OpenTSTOOL and gpml report FAIL while
% their .mexw64 files were sitting in the folders.
%
% The reliable question is simply whether a file named <fn>.<mexext> exists.
% That is asked three ways, cheapest first, and none of them depends on MATLAB's
% function-precedence rules:
%
%   1. which(fn,'-all')                -- every match on the path
%   2. the directory of the .m sibling -- where a built MEX almost always lands
%   3. a recursive walk of Toolboxes/  -- written out rather than using dir('**'),
%                                         which recurses in MATLAB but not Octave
ext = ['.' mexext];
tf = false;

% 1. all path matches
try
    w = which(fn,'-all');
    if ischar(w), w = {w}; end
    for i = 1:numel(w)
        if endsWithI(w{i},ext), tf = true; return; end
    end
catch
end

% 2. beside the .m of the same name
try
    w = which(fn);
    if ~isempty(w)
        if endsWithI(w,ext), tf = true; return; end
        if isfile(fullfile(fileparts(w),[fn ext])), tf = true; return; end
    end
catch
end

% 3. anywhere under Toolboxes/
tf = findFileBelow(fullfile(root,'Toolboxes'),[fn ext],0);
end

%-------------------------------------------------------------------------------
function tf = endsWithI(str,suffix)
tf = numel(str) >= numel(suffix) && ...
     strcmpi(str(end-numel(suffix)+1:end),suffix);
end

function tf = findFileBelow(root,name,depth)
% Recursive search, written out because dir('**') is MATLAB-only. Depth is
% capped: the deepest MEX in this tree is four levels down.
tf = false;
if depth > 6 || exist(root,'dir') ~= 7, return; end
% isfile: `name` here ends in .mexw64 and exist would answer 3, not 2, so this
% branch never fired. Steps 1 and 2 above were carrying the search alone.
if isfile(fullfile(root,name)), tf = true; return; end
d = dir(root);
for i = 1:numel(d)
    if ~d(i).isdir || strcmp(d(i).name,'.') || strcmp(d(i).name,'..'), continue; end
    if findFileBelow(fullfile(root,d(i).name),name,depth+1), tf = true; return; end
end
end

%-------------------------------------------------------------------------------
function n = countFeatureSet(root)
% Number of features the library defines, read from the feature set on disk.
n = 0;
f = fullfile(root,'FeatureSets','INP_ops_hctsa.txt');
if exist(f,'file') ~= 2, n = NaN; return; end
fid = fopen(f);
L = fgetl(fid);
while ischar(L)
    t = strtrim(L);
    if ~isempty(t) && t(1) ~= '#', n = n + 1; end
    L = fgetl(fid);
end
fclose(fid);
end

function n = totalCompiledFeatures(report)
n = 0;
for i = 1:numel(report.mex)
    if ~report.mex(i).fallback, n = n + report.mex(i).features; end
end
n = n + report.tiseanFeatures + report.tiseanFnnFeatures;
% (2.4.0) Features that need both the KSG MEX and TISEAN -- a 'mi' time delay with a
% TISEAN program or an 'fnn' embedding -- are in both counts above.
n = n - report.ksgTiseanOverlap;
end
