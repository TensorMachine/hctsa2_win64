function report = HCTSA2_CheckProtected(verbose)
% HCTSA2_CheckProtected   Has a merge reverted anything it must not?
%
%   HCTSA2_CheckProtected
%   report = HCTSA2_CheckProtected(false)
%
% Run (after startup) after every merge from upstream. Every path below is taken
% from this file's own location, so the folder you run it from does not matter
% (HCTSA2 2.2.1: it used to look in the current folder, scan nothing from
% anywhere else and still report every item intact).
%
% Taking files from upstream hctsa means taking them from a project whose goals
% differ from ours: its 3.0.0-alpha gates 41 features on toolboxes this build
% removed, and its external-program runners assume Unix. Any of that arriving
% unnoticed would quietly undo the point of this build, and it would not
% announce itself -- the code would still run, just not on the target machine.
%
% This checks protection manifest v2 (docs/HCTSA2_MergeUpstream_a7c69e5.md,
% section 3). It is deliberately blunt: it fails on anything suspicious rather
% than trying to judge intent. docs/merge_a7c69e5/tools/guard.py runs the same
% structural checks on a git revision without MATLAB.

if nargin < 1, verbose = true; end

% The checks use paths relative to the root; work from there and come back.
startDir = pwd;
cd(fileparts(mfilename('fullpath')));
restoreDir = onCleanup(@() cd(startDir)); %#ok<NASGU>

report = struct('check',{},'status',{},'detail',{});
if verbose
    fprintf('\n=========== protected-changes guard (manifest v2) ===========\n');
end

codeDirs = {'Operations','PeripheryFunctions','Calculation','PlottingAnalysis','.'};
% (HCTSA2 2.2.1: the scans of sections 1 and 2 also cover every folder under
% Toolboxes/, which they used to skip -- 0 hits today, so it was a gap and not a
% problem.)
toolboxDirs = foldersBelow('Toolboxes',0);
scanDirs = [codeDirs toolboxDirs];

%-------------------------------------------------------------------------------
%% 0. The scans below see the whole tree
%-------------------------------------------------------------------------------
% A scan that finds no files passes every check that follows it. These floors
% are well under the real counts (200 operations, 137 periphery functions, about
% 250 files under Toolboxes/) so they trip only when the scans look in the wrong
% place, not when a merge adds or drops a few files.
nOps = numel(dir(fullfile('Operations','*.m')));
nBF  = numel(dir(fullfile('PeripheryFunctions','*.m')));
nTb  = sum(cellfun(@(d) numel(dir(fullfile(d,'*.m'))),toolboxDirs));
report(end+1) = rec(sprintf('scans see the tree (%u operations, %u periphery, %u toolbox files)', ...
                            nOps,nBF,nTb), nOps >= 100 && nBF >= 50 && nTb >= 100, ...
                    'the scans below are looking in the wrong folder');

%-------------------------------------------------------------------------------
%% 1. No toolbox we removed may be gated on again
%-------------------------------------------------------------------------------
removedToolboxes = {'curve_fitting_toolbox','wavelet_toolbox', ...
                    'econometrics_toolbox','identification_toolbox', ...
                    'financial_toolbox'};
hits = {};
for i = 1:numel(removedToolboxes)
    f = grepTree(scanDirs,['BF_CheckToolbox(''' removedToolboxes{i}]);
    for k = 1:numel(f), hits{end+1} = sprintf('%s -> %s',f{k},removedToolboxes{i}); end %#ok<AGROW>
end
report(end+1) = rec('no removed-toolbox gating',isempty(hits),strjoin(hits,'; '));
% install.m gates through a list, toolboxCodes = {...}, and BF_CheckToolbox(
% toolboxCodes{i},...), which the literal pattern above cannot see. Upstream's
% list names all five removed toolboxes; ours holds only signal_toolbox.
removedPat = '(?i)''(curve_fitting|wavelet|econometrics|identification|financial)_toolbox''';
report(end+1) = rec('install.m: no removed-toolbox gating', ...
                    exist('install.m','file') == 2 && ~lineHas('install.m',removedPat,true),'install.m');

%-------------------------------------------------------------------------------
%% 2. No calls to any function of the five removed toolboxes
%-------------------------------------------------------------------------------
% MathWorks' function reference lists for the five toolboxes (fetched
% 2026-09-22), minus names that are generic methods or also exist in base
% MATLAB, Statistics, Signal Processing or Control -- the same lists as
% docs/merge_a7c69e5/tools/toolbox_funcs.py. A name counts as a call when it is
% followed by "(" or taken as a handle "@name", outside comments and strings,
% and is neither a file in this repository nor defined in the same file
% (function, argument, assigned variable, loop variable).
% Tests/ is scanned too: upstream's own unit tests once computed an expected
% value with the Wavelet Toolbox (OperationsUnitTests, WL_DetailCoeffs).
bad = findToolboxCalls([scanDirs {'Tests'}],removedToolboxFunctions());
report(end+1) = rec('no calls to removed-toolbox functions',isempty(bad),strjoin(bad,'; '));

%-------------------------------------------------------------------------------
%% 3. The BF_* replacements are present
%-------------------------------------------------------------------------------
% Every PeripheryFunctions file this fork added (none exists upstream) and still uses.
% (2.5.0: BF_aic, BF_arxstruc, BF_iddata, BF_selstruc and BF_wvarchg were deleted, as no
% caller was left once upstream f54e4d1 computed those features with its own algorithms.
% 3.0.0: likewise BF_aicbic, BF_ar, BF_archtest, BF_armax, BF_garch, BF_garchEstimate,
% BF_garchInfer, BF_garchVar, BF_gjr, BF_modwtvar, BF_n4sid, BF_n4sidAIC, BF_pe,
% BF_predict and BF_ssLoss, once upstream v3.0.0 fitted GARCH and AR models and the MODWT
% variance itself. 3.0.0 also added BF_mafdr, the Bioinformatics Toolbox's mafdr for
% TS_TopFeatures' q-values.)
mustExist = {'BF_CleanTempFiles','BF_RunTisean','BF_SweepStaleTempFiles', ...
    'BF_appcoef', ...
    'BF_autocorr','BF_besttree','BF_boxcox','BF_bsplineBasis', ...
    'BF_centfrq','BF_cwt','BF_detcoef','BF_dwt','BF_fit','BF_fitEval', ...
    'BF_fitModel','BF_fitoptions','BF_fittype','BF_fnval', ...
    'BF_idwt','BF_intwave','BF_kpsstest','BF_lbqtest','BF_mafdr','BF_modwt', ...
    'BF_parcorr','BF_pptest','BF_pptestTable', ...
    'BF_scal2frq','BF_spap2','BF_vratiotest', ...
    'BF_wavedec','BF_wavefun','BF_wenergy','BF_wentropy','BF_wfbmesti', ...
    'BF_wfilters','BF_wmaxlev','BF_wnoisest','BF_wpdec','BF_wrcoef'};
missing = {};
for i = 1:numel(mustExist)
    if exist(fullfile('PeripheryFunctions',[mustExist{i} '.m']),'file') ~= 2
        missing{end+1} = mustExist{i}; %#ok<AGROW>
    end
end
report(end+1) = rec(sprintf('%u BF_* replacements present',numel(mustExist)), ...
                    isempty(missing),strjoin(missing,', '));

%-------------------------------------------------------------------------------
%% 4. Windows runtime and build
%-------------------------------------------------------------------------------
checks = {
 'startup.m is platform-aware',        'startup.m',                               'ispc'
 'TISEAN runner: stderr kept apart',   'PeripheryFunctions/BF_RunTisean.m',       'keepStderr'
 'TISEAN runner: Java timeout',        'PeripheryFunctions/BF_RunTisean.m',       'java.lang.ProcessBuilder'
 'TISEAN runner: not-installed error', 'PeripheryFunctions/BF_RunTisean.m',       'hctsa:tisean:notInstalled'
 'TISEAN runner: Fortran path limit',  'PeripheryFunctions/BF_RunTisean.m',       'FORTRAN_NAME_LIMIT'
 'TISEAN runner: LD_LIBRARY_PATH',     'PeripheryFunctions/BF_RunTisean.m',       'LD_LIBRARY_PATH'
 'short TISEAN temp names',            'PeripheryFunctions/BF_WriteTempFile.m',   'hc%06u'
 'atomic temp-name claim',             'PeripheryFunctions/BF_WriteTempFile.m',   'createNewFile'
 'stale temp sweep in TS_Compute',     'Calculation/TS_Compute.m',                'BF_SweepStaleTempFiles'
 'Tests/ on the path (HCTSA2 tools)',  'startup.m',                               'addfcn(''Tests'')'
 'TISEAN + ripser build, static',      'HCTSA2_BuildTisean.m',                    '-static'
 'ripser in the build',                'HCTSA2_BuildTisean.m',                    '''ripser'','
};
for i = 1:size(checks,1)
    ok = fileContains(checks{i,2},checks{i,3});
    report(end+1) = rec(checks{i,1},ok,checks{i,2}); %#ok<AGROW>
end
% Every operation that writes a TISEAN/ripser temp file must remove it on every
% exit; HCTSA2_AuditTiseanCleanup checks the source for the pattern.
auditOk = false;
try
    evalc('auditOk = HCTSA2_AuditTiseanCleanup(false);');
catch
end
report(end+1) = rec('temp files cleaned on every exit',auditOk,'HCTSA2_AuditTiseanCleanup');

%-------------------------------------------------------------------------------
%% 5. Third-party source fixes
%-------------------------------------------------------------------------------
T = 'Toolboxes/Tisean_3.0.1/';
srcFixes = {
 'MS_shannon comparator',          'Toolboxes/Michael_Small/MS_shannon.c',      'return (a > b) - (a < b);'
 'MS_complexitybs comparator',     'Toolboxes/Michael_Small/MS_complexitybs.c', 'return (a > b) - (a < b);'
 'solve_chol MinGW guard',         'Toolboxes/gpml/util/solve_chol.c',          '__MINGW64__'
 'gpml make: MinGW library path',  'Toolboxes/gpml/util/make.m',                'cc = ''mingw64'''
 'false_nearest delay fix',        [T 'source_c/false_nearest.c'],              ['vcomp[i]=i%comp;' char(10) '    vemb[i]=(i/comp)*delay;']
 'lyap_r length guard',            [T 'source_c/lyap_r.c'],                     'LYAP_R_NOT_ENOUGH_POINTS'
 'nstat_z guard',                  [T 'source_c/nstat_z.c'],                    'hctsa 2.0 guard'
 'boxcount guard',                 [T 'source_c/boxcount.c'],                   'hctsa 2.0 guard'
 'poincare guard',                 [T 'source_c/poincare.c'],                   'hctsa 2.0 guard'
 'lyap_spec cannot loop forever',  [T 'source_c/lyap_spec.c'],                  '(epsmin > 0.0) ? epsmin : 1./1000.'
 'd2 scramble RNG 64-bit (D3)',    [T 'source_c/d2.c'],                         'unsigned long long rnd,rndf'
 'd2 reads inside the series',     [T 'source_c/d2.c'],                         '(unsigned long)(sn+DELAY) < length'
 'TISEAN RNG 64-bit (D3)',         [T 'source_c/routines/rand.c'],              'unsigned long long'
 'c1 skips an m without points',   [T 'source_f/c1.f'],                         'if(nmax-(m-1)*id.lt.1) goto 10'
 'd1 search cannot loop forever',  [T 'source_f/d1.f'],                         'nsweep.gt.200'
 'd1 constant-series guard',       [T 'source_f/d1.f'],                         'if(.not.(eps.gt.0.)) then'
 'c1 512-char file names',         [T 'source_f/c1.f'],                         'character*512'
 'c2d 512-char file names',        [T 'source_f/c2d.f'],                        'character*512'
 'c2g 512-char file names',        [T 'source_f/c2g.f'],                        'character*512'
 'c2t 512-char file names',        [T 'source_f/c2t.f'],                        'character*512'
 'ripser LLP64 shift',             'Toolboxes/ripser/ripser.cpp',               'index_t(1) << (8 * sizeof(index_t)'
};
for i = 1:size(srcFixes,1)
    if exist(srcFixes{i,2},'file') ~= 2
        report(end+1) = rec(srcFixes{i,1},false,[srcFixes{i,2} ' missing']); %#ok<AGROW>
    else
        report(end+1) = rec(srcFixes{i,1},fileContains(srcFixes{i,2},srcFixes{i,3}),srcFixes{i,2}); %#ok<AGROW>
    end
end
% catch22 is vendored as plain files (the target machine is offline), not a
% submodule and not fetched by startup.
gm = '';
if exist('.gitmodules','file') == 2, gm = fileread('.gitmodules'); end
report(end+1) = rec('catch22 vendored, not a submodule', ...
    exist(fullfile('Toolboxes','catch22','C'),'dir') == 7 && isempty(strfind(gm,'catch22')), ... %#ok<STREMP>
    'Toolboxes/catch22');

%-------------------------------------------------------------------------------
%% 5a. Other code of ours
%-------------------------------------------------------------------------------
ours = {
 'ZG_hmm vectorised posteriors',   'Toolboxes/ZG_hmm/ZG_hmm.m',   'bb = beta(2:T,:).*B(2:T,:);'
 'EN_Randomize: no coeffnames()',  'Operations/EN_Randomize.m',   'if isfield(c, ''c'')'
 'Mutual information without JIDT', 'Operations/IN_Initialize_MI.m', 'miCalc = BF_MutualInfoCalculator(estMethod, k);'
 'q-values without Bioinformatics', 'PlottingAnalysis/TS_TopFeatures.m', 'qVals(isValid) = BF_mafdr(pVals(isValid),''BHFDR'',true);'
};
for i = 1:size(ours,1)
    report(end+1) = rec(ours{i,1},fileContains(ours{i,2},ours{i,3}),ours{i,2}); %#ok<AGROW>
end
report(end+1) = rec('hctsa2_testdata.mat present', ...
    isfile(fullfile('Tests','hctsa2_testdata.mat')),'Tests/hctsa2_testdata.mat');

%-------------------------------------------------------------------------------
%% 5b. OpenTSTOOL must stay gone
%-------------------------------------------------------------------------------
% OpenTSTOOL coming back would mean the TSTOOL removal had been reverted.
report(end+1) = rec('OpenTSTOOL stays removed', exist('Toolboxes/OpenTSTOOL','dir') ~= 7, ...
                    'Toolboxes/OpenTSTOOL should not exist');
noTSTL = isempty(grepTreeRegex({'Operations'},'^\s*function\s.*TSTL_'));
report(end+1) = rec('no TSTL_ operations', noTSTL, '');

%-------------------------------------------------------------------------------
%% 5c. Every function the library calls must exist
%-------------------------------------------------------------------------------
% Deleting an operation is only half the job: files that CALLED it have to be
% updated too. A v2.0.0-era phase deleted NL_MS_fnn but left NL_MS_nlpe and
% SD_SurrogateTest calling it, and the failure only surfaced when those
% operations were actually run -- the code parses perfectly well.
%
% This resolves every hctsa-style identifier called in Operations/,
% PeripheryFunctions/ and the HCTSA2_*.m tooling, and reports any that no
% longer exist anywhere.
dangling = findDanglingCalls();
report(end+1) = rec('no calls to deleted functions',isempty(dangling), ...
                    strjoin(dangling,'; '));

%-------------------------------------------------------------------------------
%% 6. Our tooling
%-------------------------------------------------------------------------------
% By file, in the folder each belongs in (HCTSA2 3.0.0: the test and debug tools are
% in Tests/; the tools a user works with stay at the top level).
tools = {'HCTSA2_Analyse','HCTSA2_BuildBinaries', ...
    'HCTSA2_BuildTisean','HCTSA2_CheckInstall','HCTSA2_CheckProtected','HCTSA2_ComputeFeatures', ...
    'HCTSA2_Ingest','HCTSA2_ProbeToolchain','HCTSA2_ReportSkipped', ...
    'HCTSA2_SamplingSafeFeatures','HCTSA2_TopN','HCTSA2_WhyTiseanFails'};
testTools = {'HCTSA2_AuditTiseanCleanup','HCTSA2_TestRun','HCTSA2_TestSolveChol', ...
    'HCTSA2_TestTisean','HCTSA2_TiseanProbe','HCTSA2_ValidateKnownSignals'};
absent = [tools(cellfun(@(t) ~isfile([t '.m']),tools)), ...
          strcat('Tests/',testTools(cellfun(@(t) ~isfile(fullfile('Tests',[t '.m'])),testTools)))];
report(end+1) = rec(sprintf('%u HCTSA2_* tools present (%u in Tests/)', ...
                    numel(tools) + numel(testTools),numel(testTools)),isempty(absent), ...
                    strjoin(absent,', '));

%-------------------------------------------------------------------------------
%% Report
%-------------------------------------------------------------------------------
nFail = sum(~[report.status]);
if verbose
    for i = 1:numel(report)
        if report(i).status
            fprintf('   OK   %s\n',report(i).check);
        else
            fprintf(2,'  FAIL  %-34s %s\n',report(i).check,report(i).detail);
        end
    end
    fprintf('-------------------------------------------------------------\n');
    if nFail == 0
        fprintf('  All %u protected items intact.\n',numel(report));
    else
        fprintf(2,'  %u of %u checks FAILED -- the merge has reverted something.\n', ...
            nFail,numel(report));
    end
    fprintf('=============================================================\n\n');
end
end

%-------------------------------------------------------------------------------
function names = removedToolboxFunctions()
T = struct();
T.curve_fitting = strsplit(strtrim([ ... % Curve Fitting Toolbox
    'aptknt argnames augknt aveknt bkbrk brk2knt bspligui bspline cfit chbpnt ' ...
    'coeffnames coeffvalues confint csape csapi csaps cscvn curveFitter ' ...
    'datastats dependnames differentiate excludedata fit fitoptions fittype ' ...
    'fn2fm fnbrk fnchg fncmb fnder fndir fnint fnjmp fnmin fnplt fnrfn fntlr ' ...
    'fnval fnxtr fnzeros formula getcurve indepnames integrate islinear knt2brk ' ...
    'knt2mlt newknt numargs numcoeffs optknt ppmak predint prepareCurveData ' ...
    'prepareSurfaceData probnames probvalues rpmak rscvn rsmak setoptions sfit ' ...
    'slvblk smooth spap2 spapi spaps spcol spcrv splinetool splpp spmak sprpp ' ...
    'spterms stcol stmak subplus tpaps ' ...
    ]));
T.wavelet = strsplit(strtrim([ ... % Wavelet Toolbox
    'allnodes appcoef appcoef2 bestlevt besttree biorfilt biorwavf blscalf ' ...
    'bswfun centerFrequencies centerPeriods centfrq cfs2wpt cgauwavf cmorwavf ' ...
    'coifwavf cqt csymwavf cwt cwtfilterbank cwtfreqbounds cwtft cwtft2 ' ...
    'cwtftinfo2 cwtmag2sig daubfactors dbaux dbwavf dddtree dddtreecfs ddencmp ' ...
    'depo2ind detcoef detcoef2 dgt dtfilters dtree dualtree dwpt dwt dwt2 ' ...
    'dwtfilterbank dwtleader dwtmode dwtpassbands dyaddown dyadup entrupd ewt ' ...
    'fbspwavf fejerkorovkin filters2lp framebounds gauswavf haart hanscalf icqt ' ...
    'icwt idddtree idualtree idwpt idwt idwt2 ihaart ilwt imodwpt imodwt ' ...
    'ind2depo intwave isBiorthogonal isbiorthwfb isnode isOrthogonal isorthwfb ' ...
    'istnode iswt iswt2 itqwt iwsst leaves liftingScheme liftingStep lwt mbscalf ' ...
    'mdwtdec mdwtrec measerr mexihat meyer meyeraux modwpt modwptdetails modwt ' ...
    'modwtcorr modwtmra modwtvar modwtxcorr morlet nodeasc nodedesc nodejoin ' ...
    'nodepar nodesplt noleaves ntree orthfilt otnodes pat2cwav plotdt ' ...
    'qbiorthfilt qfactor qmf qorthwavf rbiowavf removeHarmonics scal2frq ' ...
    'scaleSpectrum scalingfunctions shanwavf swt swt2 symaux symwavf tffilt ' ...
    'thselect timeSpectrum tnodes tqwt tqwtmra treedpth treeord upcoef upwlev ' ...
    'wavedec wavedec2 wavefun wavefun2 waveinfo waveletfamilies ' ...
    'waveletScattering wavemngr wavenames waverec waverec2 wavsupport wbmpen ' ...
    'wcoherence wdcbm wden wdencmp wdenoise wenergy wentropy wextend wfbm ' ...
    'wfbmesti wfilters wkeep wkeep1 wmaxlev wmulden wnoise wnoisest wp2wtree ' ...
    'wpcoef wpcutree wpdec wpdencmp wpfun wpjoin wprcoef wprec wpspectrum wpsplt ' ...
    'wptree wpviewcf wrcoef wrev wsst wsstridge wthcoef wthresh wtmm wvarchg ' ...
    ]));
T.econometrics = strsplit(strtrim([ ... % Econometrics Toolbox
    'adftest aicbic archtest arima arma2ar arma2ma armafevd armairf autocorr ' ...
    'bayeslm bkfilter bssm cffilter chowtest collintest corrplot crosscorr ' ...
    'cusumtest dssm dtmc econometricModeler egarch egcitest fgls garch garchfit ' ...
    'garchinfer garchpred garchset gctest gjr hac hfilter hpfilter i10test infer ' ...
    'jcitest jcontest kpsstest lagmatrix LagOp lbqtest lmctest lmtest lratiotest ' ...
    'mcmix msVAR parcorr pptest price2ret recessionplot recreg regARIMA ' ...
    'ret2price simsmooth ssm summarize tsVAR ugarch var2vec varm vec2var vecm ' ...
    'vratiotest waldtest ' ...
    ]));
T.identification = strsplit(strtrim([ ... % System Identification Toolbox
    'advice aic ar armax armaxOptions aroptions arx arxOptions arxRegul arxstruc ' ...
    'bj checkFeedback chgFreqUnit compareOptions cra delayest era etfe fcat fdel ' ...
    'findstates forecastOptions fpe fselect getcov getexp getpar getpvec ' ...
    'getTrend greyest iddata iddata2timetable idfilt idfrd idgrey idinput idlti ' ...
    'idnlarx idnlhw idpar idplot idpoly idproc idresamp idss idssdata idtf ' ...
    'impulseest iv4 ivar ivstruc ivx midprefs misdata n4sid n4sidOptions nkshift ' ...
    'nlarx nlhw noise2meas noisecnv oe pe pem peOptions pexcit polydata polyest ' ...
    'predictOptions present procest realdata recursiveAR recursiveARMA ' ...
    'recursiveARX resid residOptions selstruc setcov setpar setPolyFormat ' ...
    'setpvec simsd spa spafdr spectrumest ssest ssestOptions ssform ssregest ' ...
    'struc tfest ' ...
    ]));
T.financial = strsplit(strtrim([ ... % Financial Toolbox
    'bollinger boxcox chaikvolat elpm emaxdrawdown ewstats hhigh llow lpm macd ' ...
    'maxdrawdown movavg ret2tick rsindex sharpe stochosc tick2ret tsaccel tsmom ' ...
    'tsmovavg willpctr ' ...
    ]));
names = [T.curve_fitting, T.wavelet, T.econometrics, T.identification, T.financial];
names = unique(names(~cellfun(@isempty,names)));
end

%-------------------------------------------------------------------------------
function out = findToolboxCalls(dirs,names)
% "file -> name(" for every call to a removed-toolbox function (see section 2).
out = {};
forbidden = containers.Map(names,num2cell(true(1,numel(names))));
known = containers.Map('KeyType','char','ValueType','logical');
% (HCTSA2 3.0.0: FeatureSets/ too, which holds TS_GiveMeInputFiles)
for d = {'Operations','PeripheryFunctions','Calculation','PlottingAnalysis','Toolboxes','Database','FeatureSets'}
    collectNames(d{1},known,0);
end
L0 = dir('*.m');
for i = 1:numel(L0), [~,n] = fileparts(L0(i).name); known(n) = true; end
for d = 1:numel(dirs)
    L = dir(fullfile(dirs{d},'*.m'));
    for i = 1:numel(L)
        f = fullfile(dirs{d},L(i).name);
        try
            code = stripCode(fileread(f));
        catch
            continue
        end
        local = definedNames(code);
        toks = [regexp(code,'(?<![\w.])([A-Za-z]\w*)\s*\(','tokens'), ...
                regexp(code,'@([A-Za-z]\w*)','tokens')];
        for t = 1:numel(toks)
            nm = toks{t}{1};
            if isKey(forbidden,nm) && ~isKey(known,nm) && ~any(strcmp(local,nm))
                entry = sprintf('%s -> %s(',f,nm);
                if ~any(strcmp(out,entry)), out{end+1} = entry; end %#ok<AGROW>
            end
        end
    end
end
end

function code = stripCode(src)
% Blank comments (including %{ ... %} blocks and ... continuations) and the
% contents of string literals, keeping the line structure. A quote after an
% identifier character, a closing bracket, a dot or another quote is a
% transpose, not a string.
lines = regexp(src,'\r?\n','split');
inBlock = false;
for k = 1:numel(lines)
    s = strtrim(lines{k});
    if inBlock
        if strcmp(s,'%}'), inBlock = false; end
        lines{k} = ''; continue
    end
    if strcmp(s,'%{'), inBlock = true; lines{k} = ''; continue; end
    ln = lines{k}; res = blanks(0); i = 1; n = numel(ln); prev = '';
    while i <= n
        c = ln(i);
        if c == '%', break; end
        if c == '.' && i+2 <= n && strcmp(ln(i:i+2),'...'), break; end
        isTranspose = c == '''' && ~isempty(prev) && ...
            (isstrprop(prev,'alphanum') || any(prev == '_)]}.'''));
        if c == '"' || (c == '''' && ~isTranspose)
            j = i + 1;
            while j <= n
                if ln(j) == c
                    if j < n && ln(j+1) == c, j = j + 2; continue; end
                    break
                end
                j = j + 1;
            end
            res = [res c c]; %#ok<AGROW>
            i = j + 1; prev = c;
            continue
        end
        res(end+1) = c; %#ok<AGROW>
        if ~isspace(c), prev = c; end
        i = i + 1;
    end
    lines{k} = res;
end
code = strjoin(lines,newline);
end

function names = definedNames(code)
% Functions, arguments, assigned variables (single and multiple) and loop
% variables defined anywhere in CODE.
names = {};
t = regexp(code,'(?m)^\s*function\s+(?:\[[^\]]*\]\s*=\s*|\w+\s*=\s*)?(\w+)','tokens');
names = [names, cellfun(@(x) x{1},t,'UniformOutput',false)];
t = regexp(code,'(?m)^\s*function[^(\n]*\(([^)]*)\)','tokens');
for k = 1:numel(t)
    names = [names, strtrim(strsplit(t{k}{1},','))]; %#ok<AGROW>
end
t = regexp(code,'(?m)(?:^|[;,\[])\s*([A-Za-z]\w*)\s*(?:\([^=\n]*\))?\s*=(?!=)','tokens');
names = [names, cellfun(@(x) x{1},t,'UniformOutput',false)];
t = regexp(code,'\[([^\]\n]*)\]\s*=(?!=)','tokens');
for k = 1:numel(t)
    parts = regexp(t{k}{1},'[,\s]+','split');
    names = [names, parts(~cellfun(@isempty,parts) & ~strcmp(parts,'~'))]; %#ok<AGROW>
end
t = regexp(code,'\<for\s+(\w+)\s*=','tokens');
names = [names, cellfun(@(x) x{1},t,'UniformOutput',false)];
names = unique(names(~cellfun(@isempty,names)));
end

%-------------------------------------------------------------------------------
function s = rec(check,status,detail)
s = struct('check',check,'status',logical(status),'detail',detail);
end

function tf = fileContains(f,pat)
tf = false;
if exist(f,'file') ~= 2, return; end
try
    t = strrep(fileread(f),char(13),'');   % (2.5.0) CRLF checkouts: a pattern may span lines
    tf = ~isempty(strfind(t,pat)); %#ok<STREMP>
catch
end
end

function dirs = foldersBelow(root,depth)
% ROOT and every folder under it, written out because dir('**') recurses in
% MATLAB but not in Octave.
dirs = {};
if depth > 8 || exist(root,'dir') ~= 7, return; end
dirs = {root};
D = dir(root);
for i = 1:numel(D)
    if ~D(i).isdir || strcmp(D(i).name,'.') || strcmp(D(i).name,'..'), continue; end
    dirs = [dirs foldersBelow(fullfile(root,D(i).name),depth+1)]; %#ok<AGROW>
end
end

function out = grepTree(dirs,pat)
% Files under DIRS containing PAT on a non-comment line.
out = {};
for d = 1:numel(dirs)
    L = dir(fullfile(dirs{d},'*.m'));
    for i = 1:numel(L)
        f = fullfile(L(i).folder,L(i).name);
        if lineHas(f,pat,false), out{end+1} = f; end %#ok<AGROW>
    end
end
end

function out = grepTreeRegex(dirs,pat)
out = {};
for d = 1:numel(dirs)
    L = dir(fullfile(dirs{d},'*.m'));
    for i = 1:numel(L)
        f = fullfile(L(i).folder,L(i).name);
        if lineHas(f,pat,true), out{end+1} = f; end %#ok<AGROW>
    end
end
end

function tf = lineHas(f,pat,isRegex)
tf = false;
try
    fid = fopen(f);
    L = fgetl(fid);
    while ischar(L)
        s = strtrim(L);
        if ~isempty(s) && s(1) ~= '%'
            s = regexprep(s,'%.*$','');   % drop trailing comments
            if isRegex
                hit = ~isempty(regexp(s,pat,'once'));
            else
                hit = ~isempty(strfind(s,pat)); %#ok<STREMP>
            end
            if hit, tf = true; break; end
        end
        L = fgetl(fid);
    end
    fclose(fid);
catch
    try fclose(fid); catch, end
end
end

%-------------------------------------------------------------------------------
function out = findDanglingCalls()
% Identifiers matching hctsa's naming convention that are called but resolve to
% no file anywhere in the tree (called from Operations/, PeripheryFunctions/ or
% the HCTSA2_*.m tooling, at the root and, since 3.0.0, in Tests/).
%
% Three things have to be right or this produces nothing but noise:
%   * the known-file map must be built by RECURSIVE walk. dir('**') recurses in
%     MATLAB but not in Octave, and a nearly-empty map flags every call in the
%     library as dangling.
%   * local subfunctions (SUB_*, and any function defined later in the same
%     file) are not separate files and must not be reported.
%   * struct field indexing looks exactly like a call: TS_DataMat(i,j) is a
%     field, not a function. Names known to be struct fields are excluded.
%   * MEX components have no .m file at all (ML_close_ret ships as .c and
%     .mexw64), so C and C++ sources count as definitions too.
%   * a local VARIABLE indexed with parentheses is indistinguishable from a call
%     by regex -- WL_cwt assigns SC_a then writes SC_a(:). Names assigned
%     anywhere in the same file are therefore excluded.
out = {};
known = containers.Map('KeyType','char','ValueType','logical');
% (HCTSA2 3.0.0: FeatureSets/ too, which holds TS_GiveMeInputFiles)
for d = {'Operations','PeripheryFunctions','Calculation','PlottingAnalysis','Toolboxes','Database','FeatureSets'}
    collectNames(d{1},known,0);
end
% Struct fields that share the naming convention.
for f = {'TS_DataMat','TS_Quality','TS_CalcTime','TS_LoadData'}
    known(f{1}) = true;
end

pat = '(?<![A-Za-z0-9_.])([A-Z]{2}[A-Za-z0-9]*_[A-Za-z0-9_]+)\s*\(';
% The tooling in the root folder calls operations by name too.
L0 = [dir('HCTSA2_*.m'); dir(fullfile('Tests','HCTSA2_*.m'))];
for i = 1:numel(L0), [~,n] = fileparts(L0(i).name); known(n) = true; end
for d = {'Operations','PeripheryFunctions','.','Tests'}
    if strcmp(d{1},'.'), L = dir('HCTSA2_*.m');
    elseif strcmp(d{1},'Tests'), L = dir(fullfile('Tests','HCTSA2_*.m'));
    else, L = dir(fullfile(d{1},'*.m')); end
    for i = 1:numel(L)
        f = fullfile(L(i).folder,L(i).name);
        try
            txt = fileread(f);
        catch
            continue
        end
        lines = strsplit(txt,newline);
        % Names defined inside this file (the main function and any subfunctions)
        localNames = {};
        for k = 1:numel(lines)
            % (HCTSA2 3.0.0: the output part needs its '='; without it, a function with
            % no output, 'function SUB_x(a)', was read as output 'SUB_' and name 'x')
            tok = regexp(lines{k},'^\s*function\s+(?:(?:\[[^\]]*\]|[A-Za-z0-9_]+)\s*=\s*)?([A-Za-z0-9_]+)\s*(?:\(|$)','tokens','once');
            if ~isempty(tok), localNames{end+1} = tok{1}; end %#ok<AGROW>
            % A local VARIABLE indexed with parentheses is indistinguishable from a
            % call by regex: WL_cwt assigns SC_a then writes SC_a(:). Treat anything
            % assigned in this file as local.
            asg = regexp(lines{k},'^\s*([A-Za-z][A-Za-z0-9_]*)\s*=[^=]','tokens','once');
            if ~isempty(asg), localNames{end+1} = asg{1}; end %#ok<AGROW>
        end
        for k = 1:numel(lines)
            ln = strtrim(lines{k});
            if isempty(ln) || ln(1) == '%', continue; end
            ln = regexprep(ln,'%.*$','');
            ln = regexprep(ln,'''[^'']*''','''''');   % drop string literals
            tok = regexp(ln,pat,'tokens');
            for t = 1:numel(tok)
                nm = tok{t}{1};
                if isKey(known,nm) || any(strcmp(localNames,nm)), continue; end
                entry = sprintf('%s calls %s',L(i).name,nm);
                if ~any(strcmp(out,entry)), out{end+1} = entry; end %#ok<AGROW>
            end
        end
    end
end
end

function collectNames(root,map,depth)
if depth > 6 || exist(root,'dir') ~= 7, return; end
for ext = {'*.m','*.c','*.cpp','*.mex*'}
    L = dir(fullfile(root,ext{1}));
    for i = 1:numel(L)
        [~,n] = fileparts(L(i).name);
        map(n) = true;
    end
end
D = dir(root);
for i = 1:numel(D)
    if ~D(i).isdir || strcmp(D(i).name,'.') || strcmp(D(i).name,'..'), continue; end
    collectNames(fullfile(root,D(i).name),map,depth+1);
end
end
