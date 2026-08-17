function report = HCTSA2_CheckInstall(verbose)
% HCTSA2_CheckInstall   Report exactly what this installation can and cannot do.
%
%   HCTSA2_CheckInstall
%   report = HCTSA2_CheckInstall(false)
%
% Run from the hctsa 2.0 root after `startup`. Checks MATLAB version, required
% toolboxes, every compiled component, the TISEAN executables and the Java
% library, then reports how many of the 7,702 library features are available.
%
% Nothing here modifies the installation. It answers one question: if a feature
% comes back NaN, is that a real result, a missing binary, or a broken install?

if nargin < 1, verbose = true; end

report = struct();
root = pwd;
pass = @(b) tern(b,'   OK',' FAIL');

if verbose
    fprintf('\n========== hctsa 2.0 installation check ==========\n');
    fprintf('MATLAB %s on %s\n\n',version,computer);
end

%-------------------------------------------------------------------------------
%% 1. MATLAB version and required toolboxes
%-------------------------------------------------------------------------------
v = sscanf(version('-release'),'%d');
report.matlabRelease = version('-release');
report.matlabOK = ~isempty(v) && v >= 2018;

report.statsToolbox  = license('test','Statistics_Toolbox')  == 1;
report.signalToolbox = license('test','Signal_Toolbox')      == 1;

if verbose
    fprintf('-- MATLAB and required toolboxes --\n');
    fprintf('  %s  MATLAB release %s\n',pass(report.matlabOK),report.matlabRelease);
    fprintf('  %s  Statistics and Machine Learning Toolbox\n',pass(report.statsToolbox));
    fprintf('  %s  Signal Processing Toolbox\n',pass(report.signalToolbox));
    fprintf('\n-- toolboxes hctsa 2.0 no longer needs --\n');
    for t = {'Curve_Fitting_Toolbox','Wavelet_Toolbox','Econometrics_Toolbox', ...
             'Identification_Toolbox','Financial_Toolbox'}
        have = license('test',t{1}) == 1;
        fprintf('     %-26s %s\n',t{1},tern(have,'present (unused)','absent (fine)'));
    end
end

%-------------------------------------------------------------------------------
%% 2. Compiled components
%-------------------------------------------------------------------------------
% {label, probe function, features, hasMatlabFallback}
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
C = {
 'Max_Little (fastdfa/rpde/steps)', 'ML_fastdfa_core',            10, false
 'Michael_Small',                   'MS_complexitybs',            33, false
 'Physionet sampen',                'sampen_mex',                 49, false
 'catch22',                         'catch22_DN_HistogramMode_5', 240, false
 ... % OpenTSTOOL was removed with upstream v2.0.0's TSTOOL replacement. Its 619
 ... % features are now produced by native and TISEAN-based operations, so there
 ... % is no MEX component to check for.
 'gpml (MEX = speed only)',         'solve_chol',                 72, true
};
report.mex = struct('name',{},'available',{},'features',{},'fallback',{});
mexOK = 0; mexBad = 0;
if verbose, fprintf('\n-- compiled components (MEX) --\n'); end
for i = 1:size(C,1)
    nm = C{i,1}; fn = C{i,2}; nf = C{i,3}; fb = C{i,4};
    have = mexIsBuilt(fn);
    if have || fb, mexOK = mexOK + nf; else, mexBad = mexBad + nf; end
    report.mex(end+1) = struct('name',nm,'available',have,'features',nf,'fallback',fb); %#ok<AGROW>
    if verbose
        if have
            fprintf('  %s  %-32s %4d features\n',pass(true),nm,nf);
        elseif fb
            fprintf('  %s  %-32s %4d features (MATLAB fallback in use)\n','  ---',nm,nf);
        else
            fprintf('  %s  %-32s %4d features\n',pass(false),nm,nf);
        end
    end
end

%-------------------------------------------------------------------------------
%% 3. TISEAN executables
%-------------------------------------------------------------------------------
if ispc
    tiseanDir = fullfile(root,'Toolboxes','Tisean_bin');
    probe = 'nstat_z.exe';
else
    tiseanDir = fullfile(getenv('HOME'),'bin');
    probe = 'nstat_z';
end
report.tiseanDir = tiseanDir;
report.tisean = exist(fullfile(tiseanDir,probe),'file') == 2;
report.tiseanFeatures = 276;
if verbose
    fprintf('\n-- TISEAN executables --\n');
    fprintf('  %s  %-32s %4d features\n',pass(report.tisean),tiseanDir,276);
    if ~report.tisean
        fprintf('        (see INSTALL.md section 4 -- these are standalone .exe files,\n');
        fprintf('         not MEX, and are not built by HCTSA2_BuildBinaries)\n');
    end
end

%-------------------------------------------------------------------------------
%% 4. Java library
%-------------------------------------------------------------------------------
jar = fullfile(root,'Toolboxes','infodynamics-dist','infodynamics.jar');
report.infodynamics = exist(jar,'file') == 2;
if verbose
    fprintf('\n-- Java --\n');
    fprintf('  %s  infodynamics.jar\n',pass(report.infodynamics));
end

%-------------------------------------------------------------------------------
%% 5. Live smoke test on the pure-MATLAB core
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
featuresTotal = countFeatureSet();
gated = 0;
for i = 1:numel(report.mex)
    if ~report.mex(i).available && ~report.mex(i).fallback
        gated = gated + report.mex(i).features;
    end
end
if ~report.tisean, gated = gated + report.tiseanFeatures; end
avail = max(featuresTotal - gated,0);
report.featuresAvailable = avail;
report.featuresTotal = featuresTotal;
pureMatlab = featuresTotal - totalCompiledFeatures(report);

if verbose
    fprintf('\n=================================================\n');
    fprintf('  features available : %4d of %u  (%.0f%%)\n',avail,featuresTotal, ...
        100*avail/max(featuresTotal,1));
    fprintf('    pure MATLAB      : %4d  (always)\n',pureMatlab);
    fprintf('    from MEX         : %4d\n',mexOK);
    fprintf('    from TISEAN      : %4d\n',276*report.tisean);
    if mexBad > 0 || ~report.tisean
        fprintf('  unavailable        : %4d  -- run HCTSA2_BuildBinaries, or see INSTALL.md\n', ...
            mexBad + 276*(~report.tisean));
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
function tf = mexIsBuilt(fn)
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
        if exist(fullfile(fileparts(w),[fn ext]),'file') == 2, tf = true; return; end
    end
catch
end

% 3. anywhere under Toolboxes/
tf = findFileBelow(fullfile(pwd,'Toolboxes'),[fn ext],0);
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
if exist(fullfile(root,name),'file') == 2, tf = true; return; end
d = dir(root);
for i = 1:numel(d)
    if ~d(i).isdir || strcmp(d(i).name,'.') || strcmp(d(i).name,'..'), continue; end
    if findFileBelow(fullfile(root,d(i).name),name,depth+1), tf = true; return; end
end
end

%-------------------------------------------------------------------------------
function n = countFeatureSet()
% Number of features the library defines, read from the feature set on disk.
n = 0;
f = fullfile('FeatureSets','INP_ops_hctsa.txt');
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
n = n + report.tiseanFeatures;
end
