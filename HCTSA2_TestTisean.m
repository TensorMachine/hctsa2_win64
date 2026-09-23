function report = HCTSA2_TestTisean()
% HCTSA2_TestTisean   Check the TISEAN executables and how hctsa calls them.
%
%   report = HCTSA2_TestTisean
%
% Run from the hctsa root, after startup, with the executables built
% (HCTSA2_BuildTisean).
%
%---WHAT IT CHECKS, AND WHY EACH CHECK EXISTS
%
% 1. BAD INPUT IS NOT FATAL. Short series and oversized embeddings are fed to
%    every program. Each reaches a TISEAN exit() or a Fortran stop, and each
%    must come back as an ordinary non-zero status rather than taking anything
%    with it. Five segfaults and two hangs in TISEAN itself were found this
%    way and are fixed in source_c/ and source_f/.
%
% 2. BF_RunTisean RETURNS SOMETHING. Every other kind of check looks at the
%    FILES the programs write, and a backend once passed every one of those
%    while returning an empty string to MATLAB from every call -- the files
%    were correct and hctsa could not compute a single TISEAN feature. Three
%    operations use isempty(res) as their only did-it-run test, so an empty
%    return is indistinguishable from a crash to them.
%
% 3. TEMPORARY FILES. Static audit that every operation removes its temporary
%    file on every path, then a runtime test that drives the operations down
%    their error branches and counts what is left behind. Both, because the
%    audit cannot see a runtime surprise and the runtime test only covers the
%    paths its inputs happen to reach.
%
% 4. OPERATIONS. The ten operations that call TISEAN are run on a series with
%    real structure and their outputs checked for being present and finite.

fprintf('\n============ TISEAN check ============\n');
if exist('BF_WriteTempFile','file') ~= 2
    error('HCTSA2_TestTisean:noPath','Run startup first.');
end

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
c1 = onCleanup(@() BF_CleanTempFiles(fShort)); %#ok<NASGU>
c2 = onCleanup(@() BF_CleanTempFiles(fTiny));  %#ok<NASGU>

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
 'poincare',      '-d25 -m10 -q1 -C0',        fShort
 'c1',            '-d1 -m1 -M3 -t40 -n50',    fShort
 'c1',            '-d5 -m1 -M4 -t10 -n50',    fTiny
 'c2d',           '-a2',                      fShort   % not a .c1 file at all
 'c2g',           '',                         fTiny    % not a .c2 file at all
};
nOk = 0;
for i = 1:size(C,1)
    o = [tempname '_r'];
    try
        [st,~] = BF_RunTisean(sprintf('%s %s -o %s %s',C{i,1},C{i,2},o,C{i,3}),true);
        nOk = nOk + 1;
        tag = sprintf('status %d',st);
    catch err
        tag = ['ERROR: ' regexprep(err.message,'\s+',' ')];
    end
    fprintf('  %-14s %-28s %s\n',C{i,1},C{i,2},tag);
    for e = {'','.c2','.d2','.h2','.stat','.box','.ros','.poin','.c1'}
        f = [o e{1}];
        if isfile(f), try delete(f); catch, end, end %#ok<CTCH>
    end
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

% The command shapes hctsa really builds, with the keepStderr flag each call
% site passes. c2d/c2g/c2t consume d2's and c1's output, so those run first.
BF_RunTisean(sprintf('d2 -d1 -M1,4 -t0 %s',inp),true);
BF_RunTisean(sprintf('c1 -d1 -m1 -M2 -t1 -n50 -o %s.c1 %s',inp,inp),true);
cases = { sprintf('d2 -d1 -M1,4 -t0 %s',inp),                           true
          sprintf('nstat_z -# 2 -d1 -m2 %s',inp),                       false
          sprintf('false_nearest -d1 -m1 -M1,4 -t0 %s',inp),            false
          sprintf('boxcount -M1,4 -d1 -Q2.0 -#2 -o %s.box %s',inp,inp), true
          sprintf('lyap_r -d1 -m2 -t0 -s10 -o %s.ros %s',inp,inp),      true
          sprintf('poincare -d1 -m2 -q1 -C0 -o %s.poin %s',inp,inp),    true
          sprintf('c1 -d1 -m1 -M2 -t1 -n50 -o %s.c1 %s',inp,inp),       true
          sprintf('c2d -a2 %s.c1',inp),                                 false
          sprintf('c2g %s.c2',inp),                                     false
          sprintf('c2t %s.c2',inp),                                     false };
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
fprintf('  %u of %u calls returned a non-empty string\n',nOk,size(cases,1));
r = struct('nOk',nOk,'nCase',size(cases,1));
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
ops = { 'NL_TISEAN_d2(X,1,10,0)'
        'NL_TISEAN_c1(X,1,[1,7],0.02,0.5)'
        'NL_TISEAN_fnn(X,''ac'',10,0.05,0)'
        'SY_TISEAN_nstat_z(X,4,{1,3})'
        'NL_BoxCorrDim(X,50,{''ac'',5})'
        'NL_dimensions(X,50,{''ac'',''fnn''})'
        'NL_LargestLyap(X,-1,0.1,0.01,3,{1,4})'
        'NL_PoincareSection(X,''max'',{1,3})'
        'NL_GPCorrSum(X,-1,0.5,40,20,{''ac'',''fnn''},1)'
        'NL_TakensEstimator(X,-1,0.05,0.05,{1,8})' };
end
