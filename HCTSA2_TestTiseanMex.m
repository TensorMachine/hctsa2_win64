function report = HCTSA2_TestTiseanMex(nRounds)
% HCTSA2_TestTiseanMex   Prove the TISEAN MEX files match the executables.
%
%   HCTSA2_TestTiseanMex          % 6 rounds
%   report = HCTSA2_TestTiseanMex(20)
%
% Needs BOTH backends present: run HCTSA2_BuildTisean and HCTSA2_BuildTiseanMex
% first. Run from the hctsa 2.0 root, after startup.
%
%---WHAT IT CHECKS, AND WHY EACH CHECK EXISTS
%
% 1. BYTE IDENTITY, INTERLEAVED. Every program is run through the MEX file and
%    through the executable on the same input, and the output files are
%    compared byte for byte. The cases are shuffled and repeated, because the
%    failure this is really looking for only appears on the SECOND call: a
%    TISEAN program mutates its own globals and a MEX file keeps them. Measured
%    during development -- d2 scales EPSMAX by the data interval and never
%    restores it, so call 2 returned a different, entirely plausible-looking
%    correlation sum; and c1 advances a SLATEC random-number seed that a fresh
%    process starts from zero. A single-call test passes both.
%
% 2. MEMORY. RSS is sampled across the rounds. TISEAN frees some of its globals
%    on the normal path and none of them when it exits early, so the gateway
%    tracks every allocation and sweeps what is left. Growth should be flat.
%
% 3. NO HARD FAILURES ON BAD INPUT. Short series and oversized embeddings are
%    fed to every program. In a MEX file a TISEAN exit() or a Fortran stop ends
%    MATLAB, so each of these must come back as an ordinary non-zero status.
%
% 4. OPERATION LEVEL. The ten hctsa operations that call TISEAN are run with
%    the backend forced to 'exe' and then to 'mex', and every output field is
%    compared. This is the check that matters for feature values.

if nargin < 1 || isempty(nRounds), nRounds = 6; end
if exist('BF_WriteTempFile','file') ~= 2
    error('HCTSA2_TestTiseanMex:noPath','Run startup first.');
end

progs = {'d2','nstat_z','false_nearest','boxcount','lyap_r','poincare', ...
         'c1','c2d','c2g','c2t'};
haveMex = false(1,numel(progs));
for i = 1:numel(progs), haveMex(i) = exist(['tisean_' progs{i}],'file') == 3; end
fprintf('\n============ TISEAN MEX verification ============\n');
fprintf('MEX files present : %u of 10\n',sum(haveMex));
if ~any(haveMex)
    error('HCTSA2_TestTiseanMex:noMex','No TISEAN MEX files found. Build them first.');
end

rng(1);
y = cumsum(randn(800,1)); y = (y-mean(y))/std(y);
report = struct('identity',[],'memory',[],'robustness',[],'operations',[]);

report.identity   = checkIdentity(y,nRounds,haveMex,progs);
report.memory     = report.identity.memoryKB;
report.robustness = checkRobustness(haveMex,progs);
report.operations = checkOperations();

fprintf('\n================== summary =====================\n');
fprintf('  byte identity : %u of %u comparisons matched\n', ...
        report.identity.nMatch,report.identity.nCompare);
if ispc
    % memory() reports MATLAB's whole footprint, which moves for reasons
    % unrelated to these MEX files and goes negative routinely. Reported,
    % but not treated as a measurement.
    fprintf('  memory        : %d kB change over %u calls (whole-MATLAB, noisy)\n', ...
            report.memory,report.identity.nCalls);
else
    fprintf('  memory growth : %d kB over %u calls\n', ...
            report.memory,report.identity.nCalls);
end
fprintf('  bad input     : %u of %u handled without a crash\n', ...
        report.robustness.nOk,report.robustness.nCase);
fprintf('  operations    : %u of %u identical across backends', ...
        report.operations.nMatch,report.operations.nCompare);
if report.operations.nBothErrored > 0
    fprintf(' (%u compared nothing)',report.operations.nBothErrored);
end
fprintf('\n');
allGood = report.identity.nMatch == report.identity.nCompare && ...
          report.robustness.nOk == report.robustness.nCase && ...
          report.operations.nMatch == report.operations.nCompare;
if allGood
    fprintf('  ALL CHECKS PASSED\n');
else
    fprintf(2,'  FAILURES ABOVE -- do not ship this build\n');
end
fprintf('===============================================\n\n');
end

%-------------------------------------------------------------------------------
function r = checkIdentity(y,nRounds,haveMex,progs)
fprintf('\n-- 1. byte identity, interleaved over %u rounds --\n',nRounds);
fn = BF_WriteTempFile(y);
fnShort = BF_WriteTempFile(y(1:300));

% {program, arguments, input file, output extensions ('' = the -o file itself)}
C = {
 'd2',            '-d1 -M1,3 -t10',                  fn,      {'.c2','.d2','.h2'}
 'd2',            '-d5 -M1,4 -t25',                  fn,      {'.c2','.d2','.h2'}
 'd2',            '-d1 -M1,5 -t0 -N0 -r0.01 -R2 -#30', fnShort, {'.c2'}
 'nstat_z',       '-# 4 -d1 -m3',                    fn,      {''}
 'nstat_z',       '-# 3 -d2 -m2',                    fnShort, {''}
 'false_nearest', '-d1 -m1 -M1,4 -t10 -V0',          fn,      {''}
 'false_nearest', '-d5 -m1 -M1,6 -t25 -V0',          fn,      {''}
 'boxcount',      '-M1,4 -d1 -Q2.0 -#10',            fn,      {''}
 'boxcount',      '-M1,5 -d2 -Q0.0 -#50',            fnShort, {''}
 'lyap_r',        '-d1 -m3 -t10 -s20',               fn,      {''}
 'lyap_r',        '-d2 -m6 -t40 -s50',               fn,      {''}
 'poincare',      '-d1 -m3 -q1 -C0',                 fn,      {''}
 'poincare',      '-d3 -m4 -q2 -C1',                 fnShort, {''}
 'c1',            '-d1 -m1 -M3 -t10 -n50',           fn,      {''}
 'c1',            '-d3 -m1 -M4 -t25 -n40',           fnShort, {''}
 % -M7 with a large -n, which is what NL_TISEAN_c1's own master operations use.
 % This case exists because the -M3 and -M4 cases above did NOT catch d1.f's
 % out-of-bounds shuffle: the index only goes negative for m > 1, and the higher
 % the embedding the more often it does. A test that only exercises small m
 % passes while c1 returns different numbers on every call after the first.
 'c1',            '-d1 -m1 -M7 -t12 -n300',          fn,      {''}
 'c1',            '-d1 -m2 -M6 -t12 -n300',          fn,      {''}
};
keep = ismember(C(:,1),progs(haveMex));
C = C(keep,:);

% c2d needs a .c1 and c2g/c2t need a .c2, so make them once up front.
base = [fn '.base'];
runBoth('c1',['-d1 -m1 -M3 -t10 -n50 -o ' base '.c1'],fn,'exe');
runBoth('d2',['-d1 -M1,3 -t10 -o ' base],fn,'exe');
D = {};
if isfile([base '.c1']) && any(strcmp(progs(haveMex),'c2d'))
    D(end+1,:) = {'c2d','-a2',[base '.c1'],{''}};
end
if isfile([base '.c2'])
    if any(strcmp(progs(haveMex),'c2g')), D(end+1,:) = {'c2g','',[base '.c2'],{''}}; end
    if any(strcmp(progs(haveMex),'c2t')), D(end+1,:) = {'c2t','',[base '.c2'],{''}}; end
end
C = [C; D];

nCompare = 0; nMatch = 0; nCalls = 0; mem0 = 0; failures = {};
for round = 1:nRounds
    order = randperm(size(C,1));
    for oi = order
        prog = C{oi,1}; args = C{oi,2}; inp = C{oi,3}; exts = C{oi,4};
        % Short output paths derived from the data file, NOT tempname. The
        % Fortran programs declare their filenames character*72 and silently
        % truncate anything longer -- measured under MinGW, a 73-character path
        % is written as its first 72 characters. On Windows tempname alone is
        % about 70, so this test sat one character away from reporting
        % mismatches that were nothing but a truncated name.
        oExe = [fn '.oe']; oMex = [fn '.om'];
        delTree(oExe); delTree(oMex);
        runBoth(prog,[args ' -o ' oExe],inp,'exe');
        runBoth(prog,[args ' -o ' oMex],inp,'mex');
        nCalls = nCalls + 2;
        for k = 1:numel(exts)
            a = [oExe exts{k}]; b = [oMex exts{k}];
            nCompare = nCompare + 1;
            [eq,why] = filesEqual(a,b);
            if eq
                nMatch = nMatch + 1;
            else
                failures{end+1} = sprintf('%s %s%s  [%s]',prog,args,exts{k},why); %#ok<AGROW>
            end
        end
        delTree(oExe); delTree(oMex);
    end
    if round == 2, mem0 = rssKB(); end
    fprintf('  round %2u: %u/%u matched\n',round,nMatch,nCompare);
end
% A deliberately long output path, for every Fortran program. TISEAN shipped
% these with character*72 filename buffers and truncated anything longer in
% silence; they are character*512 now. If that is ever undone, this is what
% notices -- the rest of this test uses short paths and would sail past it.
longDir = fullfile(tempdir,'hctsa_tisean_longpath_check');
if ~isfolder(longDir), mkdir(longDir); end
stem = fullfile(longDir,repmat('p',1,150));
fprintf('  long output paths (%u characters):\n',numel(stem)+2);
for p = {'c1','c2d','c2g','c2t'}
    if ~any(strcmp(progs(haveMex),p{1})), continue; end
    switch p{1}
        case 'c1',  a = '-d1 -m1 -M3 -t10 -n50'; in = fn;
        case 'c2d', a = '-a2';                   in = [base '.c1'];
        otherwise,  a = '';                      in = [base '.c2'];
    end
    if ~isfile(in), continue; end
    le = [stem '_e']; lm = [stem '_m'];
    runBoth(p{1},[a ' -o ' le],in,'exe');
    runBoth(p{1},[a ' -o ' lm],in,'mex');
    nCalls = nCalls + 2; nCompare = nCompare + 1;
    [eq,why] = filesEqual(le,lm);
    if eq
        nMatch = nMatch + 1;
        fprintf('    %-4s ok\n',p{1});
    else
        failures{end+1} = sprintf('%s, long path  [%s]',p{1},why); %#ok<AGROW>
        fprintf(2,'    %-4s FAILED: %s\n',p{1},why);
    end
    delTree(le); delTree(lm);
end
try rmdir(longDir,'s'); catch, end %#ok<CTCH>

if ~isempty(failures)
    fprintf(2,'  MISMATCHES:\n');
    u = unique(failures);
    for i = 1:numel(u), fprintf(2,'    %s\n',u{i}); end
end
delTree(base); delTree(fn); delTree(fnShort);
r = struct('nCompare',nCompare,'nMatch',nMatch,'nCalls',nCalls, ...
           'memoryKB',rssKB()-mem0,'failures',{unique(failures)});
if ispc
    fprintf('  memory: %d kB change after the first two rounds (whole-MATLAB, noisy)\n',r.memoryKB);
else
    fprintf('  memory: %d kB growth after the first two rounds\n',r.memoryKB);
end
end

%-------------------------------------------------------------------------------
function r = checkRobustness(haveMex,progs)
% Every one of these reaches a TISEAN exit() or a Fortran stop. On the
% executable path that ends a child process; in a MEX file it would end MATLAB
% unless the shim caught it. Reaching this function's end at all is most of the
% result -- the counts just say each call also came back with a status.
fprintf('\n-- 3. bad input must not be fatal --\n');
rng(2);
short = cumsum(randn(20,1));  short = (short-mean(short))/std(short);
tiny  = cumsum(randn(6,1));   tiny  = (tiny-mean(tiny))/std(tiny);
fShort = BF_WriteTempFile(short);
fTiny  = BF_WriteTempFile(tiny);

C = {
 'd2',            '-d25 -M1,10 -t10',        fShort
 'd2',            '-d1 -M1,10 -t10',         fTiny
 'nstat_z',       '-# 4 -d5 -m10',           fShort
 'nstat_z',       '-# 8 -d1 -m3',            fTiny
 'false_nearest', '-d10 -m1 -M1,10 -t10 -V0', fShort
 'boxcount',      '-M1,10 -d5 -Q2.0 -#10',   fShort
 'boxcount',      '-M1,10 -d25 -Q2.0 -#10',  fTiny
 'lyap_r',        '-d1 -m10 -t10 -s20',      fShort
 'lyap_r',        '-d5 -m10 -t40 -s10',      fShort
 'poincare',      '-d25 -m10 -q1 -C0',       fShort
 'c1',            '-d1 -m1 -M3 -t40 -n50',   fShort
 'c1',            '-d5 -m1 -M4 -t10 -n50',   fTiny
 'c2d',           '-a2',                     fShort   % not a .c1 file at all
 'c2g',           '',                        fTiny    % not a .c2 file at all
};
keep = ismember(C(:,1),progs(haveMex));
C = C(keep,:);

nOk = 0;
for i = 1:size(C,1)
    o = [tempname '_r'];
    st = NaN;
    try
        [st,~] = runBoth(C{i,1},[C{i,2} ' -o ' o],C{i,3},'mex');
        nOk = nOk + 1;
        tag = sprintf('status %d',st);
    catch err
        tag = ['ERROR: ' err.message];
    end
    fprintf('  %-14s %-28s %s\n',C{i,1},C{i,2},tag);
    delTree(o);
end
delete(fShort); delete(fTiny);
r = struct('nCase',size(C,1),'nOk',nOk);
end

%-------------------------------------------------------------------------------
function r = checkOperations()
% The check that decides whether feature values move. Each operation is run
% twice with the same input and the same seed, once per backend, and every
% output field is compared exactly.
fprintf('\n-- 4. operations, exe backend vs mex backend --\n');
% Taken verbatim from the master-operation calls in
% FeatureSets/INP_mops_hctsa.txt, so this exercises what hctsa actually invokes
% rather than a plausible-looking signature written from memory.
ops = { @(x) NL_TISEAN_d2(x,1,10,0)
        @(x) NL_TISEAN_c1(x,1,[1,7],0.02,0.5)
        @(x) NL_TISEAN_fnn(x,'ac',10,0.05,0)
        @(x) SY_TISEAN_nstat_z(x,4,{1,3})
        @(x) NL_BoxCorrDim(x,50,{'ac',5})
        @(x) NL_dimensions(x,50,{'ac','fnn'})
        @(x) NL_LargestLyap(x,-1,0.1,0.01,3,{1,4})
        @(x) NL_PoincareSection(x,'max',{1,3})
        @(x) NL_GPCorrSum(x,-1,0.5,40,20,{'ac','fnn'},1)
        @(x) NL_TakensEstimator(x,-1,0.05,0.05,{1,8}) };
names = {'NL_TISEAN_d2','NL_TISEAN_c1','NL_TISEAN_fnn','SY_TISEAN_nstat_z', ...
         'NL_BoxCorrDim','NL_dimensions','NL_LargestLyap','NL_PoincareSection', ...
         'NL_GPCorrSum','NL_TakensEstimator'};

% Deliberately NOT a random walk. cumsum(randn) has an autocorrelation time of
% the order of its own length, so BF_Embed('ac') returns a huge tau, the
% embedding window exceeds the series, and half these operations error before
% they reach TISEAN at all -- which would leave this check comparing two error
% messages and calling it a pass. Two periodic components plus noise gives a
% short correlation time and real numbers to compare.
rng(3);
t = (1:600)';
x = sin(2*pi*t/17) + 0.5*sin(2*pi*t/5) + 0.3*randn(600,1);
x = (x-mean(x))/std(x);

nCompare = 0; nMatch = 0; nBothErrored = 0;
for i = 1:numel(ops)
    a = runOp(ops{i},x,'exe');
    b = runOp(ops{i},x,'mex');
    nCompare = nCompare + 1;
    [verdict,why] = sameResult(a,b);
    switch verdict
        case 'same'
            nMatch = nMatch + 1;
            fprintf('  %-20s identical\n',names{i});
        case 'botherrored'
            % The backends agree, so it counts as a match -- but it is reported
            % separately because no feature value was compared and it therefore
            % verifies nothing about the numbers.
            nMatch = nMatch + 1; nBothErrored = nBothErrored + 1;
            fprintf('  %-20s both backends errored alike -- NOTHING COMPARED: %s\n', ...
                    names{i},why);
        otherwise
            fprintf(2,'  %-20s DIFFERS (%s)\n',names{i},why);
    end
end
setBackend('');
if nBothErrored > 0
    fprintf(2,['  %u operation(s) errored under both backends, so their feature\n' ...
               '  values were never compared. Investigate those before trusting\n' ...
               '  this run.\n'],nBothErrored);
end
r = struct('nCompare',nCompare,'nMatch',nMatch,'nBothErrored',nBothErrored);
end

function out = runOp(fh,x,backend)
setBackend(backend);
try
    rng(1);              % the stochastic operations must see the same stream
    out = fh(x);
catch err
    out = err;
end
end

function setBackend(b)
setenv('HCTSA2_TISEAN_BACKEND',b);
clear BF_RunTisean          % drop the cached choice
end

function [verdict,why] = sameResult(a,b)
% VERDICT is 'same', 'botherrored' or 'differ'.
verdict = 'differ'; why = '';
ea = isError(a); eb = isError(b);
if ea || eb
    if ea && eb
        verdict = 'botherrored';
        why = strtrim(errText(a));
        if numel(why) > 70, why = [why(1:67) '...']; end
    else
        why = 'one backend errored and the other did not';
    end
    return
end
if isstruct(a) ~= isstruct(b), why = 'different output types'; return; end
if ~isstruct(a)
    if isequaln(a,b), verdict = 'same'; else, why = 'scalar output differs'; end
    return
end
fa = sort(fieldnames(a)); fb = sort(fieldnames(b));
if ~isequal(fa,fb), why = 'different field names'; return; end
bad = {};
for i = 1:numel(fa)
    if ~isequaln(a.(fa{i}),b.(fa{i})), bad{end+1} = fa{i}; end %#ok<AGROW>
end
if isempty(bad)
    verdict = 'same';
else
    why = sprintf('%u field(s): %s',numel(bad),strjoin(bad(1:min(4,end)),', '));
end
end

function t = isError(v)
% MATLAB hands back an MException; other interpreters hand back a struct with
% the same fields. Both have to be recognised, or two error objects get compared
% field by field and their differing stacks read as a real disagreement.
t = isa(v,'MException') || (isstruct(v) && isfield(v,'message') && isfield(v,'stack'));
end

function s = errText(v)
if isa(v,'MException') || isstruct(v), s = v.message; else, s = ''; end
end

%-------------------------------------------------------------------------------
function [st,out] = runBoth(prog,args,inp,backend)
setenv('HCTSA2_TISEAN_BACKEND',backend);
clear BF_RunTisean
[st,out] = BF_RunTisean(sprintf('%s %s %s',prog,args,inp),true);
end

function [t,why] = filesEqual(a,b)
% WHY says which kind of failure it was. A missing output file and a genuine
% numerical disagreement are completely different problems, and a bare
% pass/fail count cannot tell them apart -- which cost a diagnostic round trip.
t = false; why = '';
ea = isfile(a); eb = isfile(b);
if ~ea || ~eb
    if ~ea && ~eb, why = 'neither backend produced a file';
    elseif ~ea,    why = 'no file from the exe backend';
    else,          why = 'no file from the mex backend';
    end
    return
end
x = readAll(a); y = readAll(b);
if numel(x) ~= numel(y)
    why = sprintf('sizes differ: %u vs %u bytes',numel(x),numel(y));
    return
end
d = find(x ~= y,1);
if isempty(d)
    t = true;
elseif isequal(x(x ~= 13),y(y ~= 13))
    why = 'carriage returns only';
else
    why = sprintf('first difference at byte %u of %u',d,numel(x));
end
end

function b = readAll(f)
fid = fopen(f,'r'); b = fread(fid,Inf,'*uint8'); fclose(fid); b = b(:);
end

function delTree(p)
for e = {'','.c1','.c2','.d2','.h2','.stat','.tso','.oe','.om'}
    f = [p e{1}];
    if isfile(f), try delete(f); catch, end, end %#ok<CTCH>
end
end

function kb = rssKB()
kb = 0;
try
    if ispc
        [~,m] = memory; %#ok<ASGLU>
        u = memory; kb = round(u.MemUsedMATLAB/1024);
    elseif isfile('/proc/self/statm')
        f = fopen('/proc/self/statm','r'); v = fscanf(f,'%d'); fclose(f);
        kb = v(2)*4;
    end
catch
end
end
