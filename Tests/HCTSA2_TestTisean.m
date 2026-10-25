function report = HCTSA2_TestTisean()
% HCTSA2_TestTisean   Check the TISEAN executables and how hctsa calls them.
%
%   report = HCTSA2_TestTisean
%
% Run after startup (from any folder: since HCTSA2 3.0.0 this file is in Tests/,
% which startup puts on the path), with the executables built (HCTSA2_BuildTisean).
%
%---WHAT IT CHECKS, AND WHY EACH CHECK EXISTS
%
% 1. BAD INPUT IS NOT FATAL. Short series and oversized embeddings are fed to
%    every program. Each reaches a TISEAN exit() or a Fortran stop, and each
%    must come back as an ordinary exit status (0 to 127: TISEAN's deliberate
%    exit(n) codes) rather than taking anything with it. A status of 128 or more,
%    or a negative one (a signal on Linux/macOS, an NTSTATUS crash code on
%    Windows), counts as a crash. Five segfaults and two hangs in TISEAN itself
%    were found this way and are fixed in source_c/ and source_f/.
%
% 2. BF_RunTisean RETURNS SOMETHING. Every other kind of check looks at the
%    FILES the programs write, and a backend once passed every one of those
%    while returning an empty string to MATLAB from every call -- the files
%    were correct and hctsa could not compute a single TISEAN feature. Three
%    operations use isempty(res) as their only did-it-run test, so an empty
%    return is indistinguishable from a crash to them. Also: lyap_spec must be
%    built from 2.2.0's source. An older build accepts NL_LyapSpec's -t (Theiler
%    window) without complaint and ignores it, so its values would be wrong
%    with no error anywhere.
%
% 3. TEMPORARY FILES. Static audit that every operation removes its temporary
%    file on every path, then a runtime test that drives the operations down
%    their error branches and counts what is left behind. Both, because the
%    audit cannot see a runtime surprise and the runtime test only covers the
%    paths its inputs happen to reach.
%
% 4. OPERATIONS. The eleven operations that call TISEAN (or ripser, which
%    NL_PersistentHomology runs the same way) are run on a series with real
%    structure and their outputs checked for being present and finite.

fprintf('\n============ TISEAN check ============\n');
if exist('BF_WriteTempFile','file') ~= 2
    error('HCTSA2_TestTisean:noPath','Run startup first.');
end
% (HCTSA2 2.2.1: the checks below re-seed the global random stream with rng(1),
% rng(2) and rng(3); leave it as the caller had it.)
rngState = rng;
restoreRng = onCleanup(@() rng(rngState)); %#ok<NASGU>

report = struct();
report.robustness = checkRobustness();
report.returns    = checkReturns();
report.cleanupOk  = HCTSA2_AuditTiseanCleanup(true);
report.leaks      = checkLeaks();
report.operations = checkOperations();

fprintf('\n================ summary =================\n');
fprintf('  bad input     : %u of %u handled without a crash\n', ...
        report.robustness.nOk,report.robustness.nCase);
fprintf('  returned data : %u of %u calls returned a non-empty string\n', ...
        report.returns.nOk,report.returns.nCase);
if report.cleanupOk
    fprintf('  temp cleanup  : every operation guards every exit\n');
else
    fprintf(2,'  temp cleanup  : AN OPERATION CAN LEAK -- see the audit above\n');
end
fprintf('  temp leaks    : %u file(s) left after %u calls (%u errored)\n', ...
        report.leaks.nLeft,report.leaks.nRun,report.leaks.nErr);
fprintf('  operations    : %u of %u produced usable output\n', ...
        report.operations.nOk,report.operations.nCase);

allGood = report.robustness.nOk == report.robustness.nCase && ...
          report.returns.nOk    == report.returns.nCase && ...
          report.cleanupOk && report.leaks.nLeft == 0 && ...
          report.operations.nOk == report.operations.nCase;
if allGood
    fprintf('  ALL CHECKS PASSED\n');
else
    fprintf(2,'  SOMETHING FAILED -- see above\n');
end
fprintf('==========================================\n');
end

%-------------------------------------------------------------------------------
function r = checkRobustness()
fprintf('\n-- 1. bad input must not be fatal --\n');
rng(2);
short = cumsum(randn(20,1)); short = (short-mean(short))/std(short);
tiny  = cumsum(randn(6,1));  tiny  = (tiny-mean(tiny))/std(tiny);
fShort = BF_WriteTempFile(short);
fTiny  = BF_WriteTempFile(tiny);
% (2.5.0) a constant series: with zero spread, c1's box size was zero and the
% neighbour search indexed out of bounds (a segfault found by fuzzing the merge)
fConst = BF_WriteTempFile(ones(866,1));
c1 = onCleanup(@() BF_CleanTempFiles(fShort)); %#ok<NASGU>
c2 = onCleanup(@() BF_CleanTempFiles(fTiny));  %#ok<NASGU>
c3 = onCleanup(@() BF_CleanTempFiles(fConst)); %#ok<NASGU>

C = {
 'd2',            '-d25 -M1,10 -t10',         fShort
 'd2',            '-d1 -M1,10 -t10',          fTiny
 'nstat_z',       '-# 4 -d5 -m10',            fShort
 'nstat_z',       '-# 8 -d1 -m3',             fTiny
 'false_nearest', '-d10 -m1 -M1,10 -t10 -V0', fShort
 'boxcount',      '-M1,10 -d5 -Q2.0 -#10',    fShort
 'boxcount',      '-M1,10 -d25 -Q2.0 -#10',   fTiny
 'lyap_r',        '-d1 -m10 -t10 -s20',       fShort
 'lyap_r',        '-d5 -m10 -t40 -s10',       fShort
 'c1',            '-d1 -m1 -M3 -t40 -n50',    fShort
 'c1',            '-d5 -m1 -M4 -t10 -n50',    fTiny
 'c1',            '-d10 -m1 -M4 -t58 -n215',  fConst
 'c2d',           '-a2',                      fShort   % not a .c1 file at all
 'c2g',           '',                         fTiny    % not a .c2 file at all
 'lyap_spec',     '-m1,10 -d25 -k30',         fShort
 'lyap_spec',     '-m1,3 -d1 -k30',           fTiny
 'nrlazy',        '-m1,5 -d5 -i2 -v0.3',      fShort
 'nrlazy',        '-m1,5 -d1 -i1 -v0.3',      fTiny
 'ripser',        '--format point-cloud --dim 2', fTiny
};
nOk = 0;
for i = 1:size(C,1)
    % (HCTSA2 2.2.1: an hc###### stem, so that BF_CleanTempFiles removes every
    % file a program writes beside it -- d2's .c2/.d2/.h2/.stat, nrlazy's .1,
    % .2 -- and the stale-file sweep can find what a crash leaves.)
    o = BF_WriteTempFile(0);
    try
        if strcmp(C{i,1},'ripser')           % ripser has no -o: output is stdout
            cmd = sprintf('%s %s "%s"',C{i,1},C{i,2},C{i,3});
        else
            cmd = sprintf('%s %s -o "%s" "%s"',C{i,1},C{i,2},o,C{i,3});
        end
        [st,~] = BF_RunTisean(cmd,true);
        % BF_RunTisean returns whatever status the program ended with, so a
        % crash is a status, not an error: TISEAN's deliberate exit(n) codes are
        % all below 128, a signal death is 128 or more and a Windows crash code
        % is negative.
        if st >= 0 && st < 128
            nOk = nOk + 1;
            tag = sprintf('status %d',st);
        else
            tag = sprintf('CRASHED (status %d)',st);
        end
    catch err
        tag = ['ERROR: ' regexprep(err.message,'\s+',' ')];
    end
    fprintf('  %-14s %-30s %s\n',C{i,1},C{i,2},tag);
    BF_CleanTempFiles(o);
end
fprintf('  %u of %u handled without a crash\n',nOk,size(C,1));
r = struct('nOk',nOk,'nCase',size(C,1));
end

%-------------------------------------------------------------------------------
function r = checkReturns()
fprintf('\n-- 2. what BF_RunTisean returns to the caller --\n');
rng(1);
y = cumsum(randn(800,1)); y = (y-mean(y))/std(y);
inp = BF_WriteTempFile(y);
cleanup = onCleanup(@() BF_CleanTempFiles(inp)); %#ok<NASGU>
inpE = BF_WriteTempFile([y(1:end-2) y(2:end-1) y(3:end)],7,' '); % NL_LyapSpec's embedding
cleanupE = onCleanup(@() BF_CleanTempFiles(inpE)); %#ok<NASGU>

% The command shapes hctsa really builds, with the keepStderr flag each call
% site passes. c2d/c2g/c2t consume d2's and c1's output, so those run first.
BF_RunTisean(sprintf('d2 -d1 -M1,4 -t0 %s',inp),true);
BF_RunTisean(sprintf('c1 -d1 -m1 -M2 -t1 -n50 -o "%s.c1" "%s"',inp,inp),true);
cases = { sprintf('d2 -d1 -M1,4 -t0 %s',inp),                           true
          sprintf('nstat_z -# 2 -d1 -m2 %s',inp),                       false
          sprintf('false_nearest -d1 -m1 -M1,4 -t0 %s',inp),            false
          sprintf('boxcount -M1,4 -d1 -Q2.0 -#2 -o "%s.box" "%s"',inp,inp), true
          sprintf('lyap_r -d1 -m2 -t0 -s10 -o "%s.ros" "%s"',inp,inp),      true
          sprintf('c1 -d1 -m1 -M2 -t1 -n50 -o "%s.c1" "%s"',inp,inp),       true
          sprintf('c2d -a2 %s.c1',inp),                                 false
          sprintf('c2g %s.c2',inp),                                     false
          sprintf('c2t %s.c2',inp),                                     false
          sprintf('lyap_spec -m3,1 -c1,2,3 -k30 -t5 -o "%s.lyap" "%s"',inpE,inpE), true
          sprintf('nrlazy -m1,3 -d1 -i1 -v0.3 -V4 %s',inp),             false
          sprintf('ripser --format point-cloud --dim 1 %s',inp),        false };
nOk = 0;
for i = 1:size(cases,1)
    why = '';
    try
        [st,out] = BF_RunTisean(cases{i,1},cases{i,2});
        if isempty(out), why = sprintf('EMPTY string back (status %d)',st); end
    catch err
        why = ['errored: ' regexprep(err.message,'\s+',' ')];
    end
    if isempty(why)
        nOk = nOk + 1;
    else
        fprintf(2,'  %-14s %s\n',strtok(cases{i,1}),why);
    end
end
% lyap_spec's own usage text names -t only in a build from 2.2.0's source
try
    [~,usage] = BF_RunTisean('lyap_spec -h',true);
    why = '';
    if isempty(strfind(usage,'Theiler'))
        why = ['built from source older than 2.2.0: it ignores -t, so NL_LyapSpec ' ...
               'is wrong. Rebuild with HCTSA2_BuildTisean.'];
    end
catch err
    why = ['errored: ' regexprep(err.message,'\s+',' ')];
end
if isempty(why)
    nOk = nOk + 1;
else
    fprintf(2,'  %-14s %s\n','lyap_spec -h',why);
end
% (HCTSA2 3.0.0) false_nearest writes a fifth column, the number of points behind each
% fraction, only when built from 3.0.0's source; NL_FNN refuses an older build
try
    [~,res] = BF_RunTisean(sprintf('false_nearest -d1 -m1 -M1,4 -t0 -V0 %s',inp));
    line = strtrim(strtok(res,newline));
    why = '';
    if isempty(line)
        why = 'printed nothing (expected 5 columns)';
    elseif numel(strsplit(line)) < 5
        why = sprintf(['built from source older than 3.0.0 (%u columns, not 5): NL_FNN ' ...
               'refuses it, so every ''fnn'' feature fails. Rebuild with HCTSA2_BuildTisean.'], ...
               numel(strsplit(line)));
    end
catch err
    why = ['errored: ' regexprep(err.message,'\s+',' ')];
end
if isempty(why)
    nOk = nOk + 1;
else
    fprintf(2,'  %-14s %s\n','false_nearest',why);
end
nCase = size(cases,1) + 2;
fprintf('  %u of %u calls returned a non-empty string (lyap_spec -t and false_nearest''s fifth column included)\n',nOk,nCase);
r = struct('nOk',nOk,'nCase',nCase);
end

%-------------------------------------------------------------------------------
function r = checkLeaks()
% Drive the operations down their ERROR branches and count what survives. The
% inputs are chosen for that: the degenerate-data fast paths ("out = NaN;
% return") are where three separate attempts at this cleanup leaked.
fprintf('\n-- 3b. temporary files, with errors forced --\n');
rng(3);
t = (1:600)';
good = sin(2*pi*t/17)+0.5*sin(2*pi*t/5)+0.3*randn(600,1);
good = (good-mean(good))/std(good);
inputs = { 'good', good; 'short10', (1:10)'+0.01*randn(10,1); 'short3', [1;2;3]
           'constant', zeros(200,1); 'twovalue', double(rand(300,1)>0.5); 'tiny', [1;2] };
ops = opList();
td = tempdir; nRun = 0; nErr = 0;
before = numel(dir(fullfile(td,'hc*')));
for a = 1:size(inputs,1)
    X = inputs{a,2}; %#ok<NASGU>
    for b = 1:numel(ops)
        nRun = nRun + 1;
        try, v = eval(ops{b}); catch, nErr = nErr + 1; end %#ok<NASGU>
    end
end
nLeft = numel(dir(fullfile(td,'hc*'))) - before;
if nLeft > 0
    fprintf(2,'  %u temporary file(s) LEAKED over %u calls\n',nLeft,nRun);
else
    fprintf('  0 leaked over %u calls (%u of them errored)\n',nRun,nErr);
end
r = struct('nRun',nRun,'nErr',nErr,'nLeft',max(nLeft,0));
end

%-------------------------------------------------------------------------------
function r = checkOperations()
fprintf('\n-- 4. the operations that call TISEAN --\n');
% Two periodic components plus noise: a short correlation time, so BF_Embed
% returns a sane tau and the operations reach TISEAN instead of erroring on
% the embedding. cumsum(randn) would not.
rng(3);
t = (1:600)';
X = sin(2*pi*t/17) + 0.5*sin(2*pi*t/5) + 0.3*randn(600,1);
X = (X-mean(X))/std(X); %#ok<NASGU>
ops = opList();
nOk = 0;
for i = 1:numel(ops)
    name = regexprep(ops{i},'\(.*','');
    try
        v = eval(ops{i});
        if isstruct(v)
            f = fieldnames(v);
            nFinite = sum(cellfun(@(k) isnumeric(v.(k)) && isscalar(v.(k)) && ...
                                       isfinite(v.(k)), f));
            good = nFinite > 0;
            tag = sprintf('%u fields, %u finite',numel(f),nFinite);
        else
            good = isnumeric(v) && isscalar(v) && isfinite(v);
            tag = sprintf('scalar %g',v);
        end
        if good, nOk = nOk + 1; fprintf('  %-22s %s\n',name,tag);
        else,    fprintf(2,'  %-22s no usable output (%s)\n',name,tag); end
    catch err
        fprintf(2,'  %-22s ERROR: %s\n',name,regexprep(err.message,'\s+',' '));
    end
end
fprintf('  %u of %u produced usable output\n',nOk,numel(ops));
r = struct('nOk',nOk,'nCase',numel(ops));
end

%-------------------------------------------------------------------------------
function ops = opList()
% Signatures taken from the master-operation calls in
% FeatureSets/INP_mops_hctsa.txt, so this exercises what hctsa really invokes.
% (NL_TakensEstimator no longer calls TISEAN: upstream made it native.)
ops = { 'NL_d2(X,1,10,{''ac'',1})'
        'NL_c1(X,1,[1,7],0.02,0.5)'
        'NL_FNN(X,''ac'',10,{''ac'',1},0,[],5)'
        'SY_nstat_z(X,4,{1,3})'
        'NL_BoxCountEntropyRate(X,50,{''ac1e'',5})'
        'NL_Dimensions(X,50,{''ac1e'',''fnn''})'
        'NL_LargestLyap(X,-1,{''ac1e'',30},{''ac1e'',1},3,{1,4})'
        'NL_GPCorrSum(X,-1,0.5,{''ac'',1},20,{''ac'',''fnn''},1)'
        'NL_LyapSpec(X,1,3,30,''full'',{''ac'',1})'
        'PP_SchreiberDenoise(X,5,1,1,0.3)'
        'NL_PersistentHomology(X)' };
end
