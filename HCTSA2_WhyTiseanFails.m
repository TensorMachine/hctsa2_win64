function HCTSA2_WhyTiseanFails()
% HCTSA2_WhyTiseanFails   Find the layer at which TISEAN stops working.
%
%   >> HCTSA2_WhyTiseanFails
%
% Run from the hctsa root, after startup, in the SAME MATLAB session that fails.
%
%---WHY THIS EXISTS
% When TISEAN breaks, every one of the ten programs breaks at once and each
% operation reports its own downstream symptom -- "did not produce a .box output
% file", "not generated?", "Call to TISEAN function failed". None of them names
% the cause. This walks the call stack from the bottom up and stops at the first
% layer that fails, so the answer is the layer number rather than a guess.
%
%   1. every executable is reachable and runs
%   2. BF_RunTisean returns data
%   3. the operation, called directly
%   4. the operation through BF_ParEval  (what hctsa evaluates through)
%   5. the operation through TS_ComputeMasterLoop, serially
%   6. the same on a parallel worker, if a pool is open
%
% Layers 3-6 are the same code path with more wrapping around it. If 3 passes
% and 5 fails, the wrapper is at fault. If 1 fails, nothing above it can work.
%
% Layer 1 RUNS each program rather than just looking for it. A name that
% resolves on the PATH is not the same as a program that executes, and the
% difference stays invisible until something calls it.

progs = {'d2','nstat_z','false_nearest','boxcount','lyap_r','poincare', ...
         'c1','c2d','c2g','c2t'};

fprintf('\n================ why TISEAN fails ================\n');

%-------------------------------------------------------------------------------
%% 0. Environment
%-------------------------------------------------------------------------------
fprintf('\n--- 0. environment ---\n');
fprintf('  MATLAB            : %s (%s)\n',version,computer);
fprintf('  pwd               : %s\n',pwd);
fprintf('  tempdir           : %s (%u chars)\n',tempdir,numel(tempdir));
fprintf('  tempname length   : %u chars\n',numel(tempname));
tmo = getenv('HCTSA2_TISEAN_TIMEOUT');
if isempty(tmo)
    fprintf('  TISEAN timeout    : 300 s (default)\n');
else
    fprintf('  TISEAN timeout    : %s s\n',tmo);
end
% Can we actually write where TISEAN puts its output?
try
    probe = [tempname '.probe']; f = fopen(probe,'w');
    if f < 0
        fprintf(2,'  tempdir writable  : NO -- fopen failed. Nothing can work.\n');
    else
        fwrite(f,'x'); fclose(f); delete(probe);
        fprintf('  tempdir writable  : yes\n');
    end
catch e
    fprintf(2,'  tempdir writable  : NO (%s)\n',e.message);
end
p = currentPool();
if isempty(p)
    fprintf('  parallel pool     : none open\n');
else
    fprintf('  parallel pool     : %s, %u workers\n',class(p),p.NumWorkers);
end

%-------------------------------------------------------------------------------
%% 1. Does each executable run?
%-------------------------------------------------------------------------------
fprintf('\n--- 1. executables: present, and do they actually run? ---\n');
% A real data file, because several TISEAN programs fall back to READING STDIN
% when their input file is missing, which would hang MATLAB rather than fail.
y = cumsum(randn(300,1)); y = (y-mean(y))/std(y);
dataFile = BF_WriteTempFile(y);
cleanupData = onCleanup(@() cleanupFiles(dataFile));
capture = [tempname '.err'];

minimalArgs = struct( ...
    'd2',            {{'-d1','-M1,3','-t1'}}, ...
    'nstat_z',       {{'-#','2','-d1','-m2'}}, ...
    'false_nearest', {{'-d1','-m1','-M1,3','-t1'}}, ...
    'boxcount',      {{'-M1,3','-d1'}}, ...
    'lyap_r',        {{'-d1','-m2','-t1','-s10'}}, ...
    'poincare',      {{'-d1','-m2','-q1','-C0'}}, ...
    'c1',            {{'-d1','-m1','-M2','-t1','-n50'}}, ...
    'c2d',           {{'-a2'}}, ...
    'c2g',           {{}}, ...
    'c2t',           {{}});

nLoad = 0; firstLoadErr = '';
for i = 1:numel(progs)
    args = [minimalArgs.(progs{i}), {'-o',[dataFile '.probe'], dataFile}];
    try
        [st,out] = BF_RunTisean(sprintf('%s %s',progs{i}, ...
                       strjoin(args,' ')),true); %#ok<ASGLU>
        nLoad = nLoad + 1;
        fprintf('  %-16s runs (status %d, %u chars back)\n', ...
                progs{i},st,numel(out));
    catch err
        if isempty(firstLoadErr), firstLoadErr = err.message; end
        fprintf(2,'  %-16s *** WILL NOT RUN ***\n',progs{i});
        fprintf(2,'      %s\n',oneLine(err.message));
    end
end
cleanupFiles(capture);
fprintf('  --> %u of 10 run\n',nLoad);
if nLoad < 10
    fprintf(2,'\n  STOP. Layer 1 failed, so nothing above it can work.\n');
    fprintf(2,['  The executables are build artifacts and do not travel with the\n' ...
               '  repository. Build them on THIS machine with HCTSA2_BuildTisean,\n' ...
               '  then run startup.\n']);
    return
end

%-------------------------------------------------------------------------------
%% 2. BF_RunTisean
%-------------------------------------------------------------------------------
fprintf('\n--- 2. BF_RunTisean (the dispatcher hctsa calls) ---\n');
cmds = { sprintf('d2 -d1 -M1,3 -t1 %s',dataFile),                true
         sprintf('nstat_z -# 2 -d1 -m2 %s',dataFile),            false
         sprintf('false_nearest -d1 -m1 -M1,3 -t1 %s',dataFile), false
         sprintf('boxcount -M1,3 -d1 -o %s.box %s',dataFile,dataFile), true
         sprintf('c2g %s.c2',dataFile),                          false };
nRun = 0;
for i = 1:size(cmds,1)
    c = cmds{i,1};
    try
        [st,out] = BF_RunTisean(c,cmds{i,2});
        if isempty(out)
            fprintf(2,'  %-14s status %-4d returned NOTHING\n',strtok(c),st);
        else
            nRun = nRun + 1;
            fprintf('  %-14s status %-4d %u chars back\n',strtok(c),st,numel(out));
        end
    catch err
        fprintf(2,'  %-14s ERROR: %s\n',strtok(c),oneLine(err.message));
    end
end
fprintf('  --> %u of %u returned data\n',nRun,size(cmds,1));

%-------------------------------------------------------------------------------
%% 3-6. The operation, with progressively more of hctsa wrapped around it
%-------------------------------------------------------------------------------
rng(3);
t = (1:600)';
xr  = sin(2*pi*t/17) + 0.5*sin(2*pi*t/5) + 0.3*randn(600,1);
x_z = (xr-mean(xr))/std(xr);
code = 'SY_TISEAN_nstat_z(x_z,4,{1,3})';

fprintf('\n--- 3. the operation, called directly (what the test script does) ---\n');
r3 = attempt(@() evalIn(code,xr,x_z));

fprintf('\n--- 4. through BF_ParEval, evalc branch (hctsa''s evaluator) ---\n');
r4 = attempt(@() BF_ParEval(xr,x_z,code,0));

fprintf('\n--- 5. through BF_ParEval, plain eval branch (howVocal = ''full'') ---\n');
r5 = attempt(@() BF_ParEval(xr,x_z,code,1));

fprintf('\n--- 6. through TS_ComputeMasterLoop, SERIAL (line 234) ---\n');
r6 = attempt(@() TS_ComputeMasterLoop(xr,x_z,code,1,1,'minimal',1,1));

r7 = [];
if ~isempty(currentPool())
    fprintf('\n--- 7. through TS_ComputeMasterLoop on a WORKER (line 224) ---\n');
    out = cell(1,1);
    try
        parfor k = 1:1
            out{k} = TS_ComputeMasterLoop(xr,x_z,code,1,1,'minimal',1,1); %#ok<PFBNS>
        end
        r7 = ~isempty(out{1});
        if r7, fprintf('  ok\n'); else, fprintf(2,'  returned empty\n'); end
    catch err
        r7 = false;
        fprintf(2,'  ERROR: %s\n',oneLine(err.message));
    end
end

%-------------------------------------------------------------------------------
fprintf('\n===================== verdict =====================\n');
lbl = {'3 direct','4 evalc','5 eval','6 serial loop'};
val = [r3 r4 r5 r6];
if ~isempty(r7), lbl{end+1} = '7 worker'; val(end+1) = r7; end
for i = 1:numel(val)
    fprintf('  %-16s %s\n',lbl{i},ternary(val(i),'ok','FAILED'));
end
if all(val)
    fprintf(['\n  Every layer works in this session. The failure is not in how\n' ...
             '  hctsa calls TISEAN -- it is specific to your data or to a\n' ...
             '  parameter combination these defaults do not reach. Run your real\n' ...
             '  computation and send the first "---Error evaluating ..." line.\n']);
elseif val(1) && ~all(val)
    fprintf(2,['\n  The operation works alone but fails once hctsa wraps it. The\n' ...
               '  first FAILED row above is the wrapper at fault.\n']);
else
    fprintf(2,'\n  The operation fails even called directly -- layers 1 and 2 above say why.\n');
end
fprintf('===================================================\n');
end

%-------------------------------------------------------------------------------
function p = currentPool()
% gcp does not exist without the Parallel Computing Toolbox (or in Octave).
p = [];
try
    p = gcp('nocreate');
catch
end
end

%-------------------------------------------------------------------------------
function out = evalIn(code,x,x_z) %#ok<INUSD>
out = eval(code);
end

function ok = attempt(fn)
try
    v = fn();
    ok = ~isempty(v) && ~(iscell(v) && isempty(v));
    if ok
        fprintf('  ok\n');
    else
        fprintf(2,'  returned EMPTY (hctsa records this as a failed operation)\n');
    end
catch err
    ok = false;
    fprintf(2,'  ERROR: %s\n',err.message);
    for k = 1:min(4,numel(err.stack))
        fprintf(2,'      at %s (line %u)\n',err.stack(k).name,err.stack(k).line);
    end
end
end

function s = oneLine(m)
s = regexprep(m,'\s+',' ');
if numel(s) > 200, s = [s(1:197) '...']; end
end

function v = ternary(c,a,b)
if c, v = a; else, v = b; end
end

function cleanupFiles(base)
ws = warning('off','MATLAB:DELETE:FileNotFound');
restore = onCleanup(@() warning(ws)); %#ok<NASGU>
for e = {'','.c2','.d2','.h2','.box','.tso','.err','.c1','.stat','.ros', ...
         '.probe','.probe.tso'}
    f = [base e{1}];
    if isfile(f), try delete(f); catch, end, end %#ok<CTCH>
end
end
