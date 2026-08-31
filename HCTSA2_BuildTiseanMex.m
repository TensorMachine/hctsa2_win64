function report = HCTSA2_BuildTiseanMex(outDir)
% HCTSA2_BuildTiseanMex   Build the TISEAN programs hctsa needs as MEX files.
%
%   HCTSA2_BuildTiseanMex
%   report = HCTSA2_BuildTiseanMex('Toolboxes/Tisean_mex')
%
% Produces tisean_d2, tisean_nstat_z, tisean_false_nearest, tisean_boxcount,
% tisean_lyap_r, tisean_poincare, tisean_c1, tisean_c2d, tisean_c2g and
% tisean_c2t as MEX files instead of standalone executables. Run from the
% hctsa 2.0 root. HCTSA2_BuildTisean (the .exe build) still works and is kept
% as a fallback; BF_RunTisean prefers whichever is present.
%
%---WHY
% The executables are the only .exe files hctsa 2.0 produces. On a machine
% whose policy restricts executables but not the DLLs MATLAB loads, MEX files
% get the same features without one. Removing system() also removes the
% cmd.exe launch per call -- roughly fifty per time series across the 29
% TISEAN master operations.
%
%---WHAT IT TAKES
% A TISEAN program is written to run once and exit. Three things have to change
% before one can be called repeatedly inside MATLAB, and all three are handled
% at build time without touching the shared sources:
%
%   1. main() is renamed with -Dmain=tisean_<name>_main (C) or the program is
%      rewritten as a subroutine (Fortran, generated below).
%   2. exit() and `stop` are routed back to the gateway instead of ending the
%      process. See Toolboxes/Tisean_3.0.1/mex/tisean_mex.h.
%   3. The file-scope globals are restored between calls. Without that, d2
%      returns different numbers on its second call. See tisean_reset.h.
%
%---REQUIREMENTS
% MATLAB R2024a or later with the MinGW-w64 add-on. R2024a is the first release
% in which MathWorks supports gfortran for MEX on Windows, and four of the ten
% programs are Fortran. On an earlier release the six C programs still build
% and the Fortran four are skipped -- use HCTSA2_BuildTisean for those.

if nargin < 1 || isempty(outDir)
    outDir = fullfile('Toolboxes','Tisean_mex');
end

root = pwd;
src  = fullfile(root,'Toolboxes','Tisean_3.0.1');
mexd = fullfile(src,'mex');
if exist(fullfile(src,'source_c'),'dir') ~= 7
    error('HCTSA2_BuildTiseanMex:noSource', ...
        'TISEAN source not found at %s. Run this from the hctsa 2.0 root.',src);
end
if exist(outDir,'dir') ~= 7, mkdir(outDir); end
absOut = fullfile(root,outDir);

fprintf('\n========= building TISEAN as MEX files =========\n');
fprintf('source : %s\n',src);
fprintf('output : %s\n',absOut);
fprintf('MATLAB %s on %s\n',version,computer);

%-------------------------------------------------------------------------------
%% Toolchain
%-------------------------------------------------------------------------------
cc = [];
try
    cc = mex.getCompilerConfigurations('C','Selected');
catch
end
if isempty(cc)
    error('HCTSA2_BuildTiseanMex:noCompiler', ...
        'No C compiler is configured for MEX. Run "mex -setup C" first.');
end
fprintf('C compiler : %s\n',cc(1).Name);

[gfortran,fortranWhy] = findGfortran(cc(1));
if isempty(gfortran)
    fprintf(2,'gfortran   : NOT USABLE (%s)\n',fortranWhy);
    fprintf(2,'  c1, c2d, c2g and c2t will be skipped. NL_TISEAN_c1,\n');
    fprintf(2,'  NL_TISEAN_d2 and NL_TakensEstimator need them.\n');
else
    fprintf('gfortran   : %s\n',gfortran);
end

build = fullfile(tempdir,sprintf('tisean_mexbuild_%06u',randi(999999)));
mkdir(build);
cleanup = onCleanup(@() rmdir(build,'s'));
fprintf('scratch    : %s\n',build);

% Off Windows a MEX file is an ELF shared object and the dynamic loader
% interposes global symbols across them: d2, nstat_z, boxcount, lyap_r,
% false_nearest and poincare all define `series`, `verbosity` and `length` at
% file scope, so whichever loaded first would win for all of them. Measured on
% Linux: the twelfth interleaved call segfaulted. Windows DLLs never share
% symbols this way, so ispc passes nothing extra and the supported build stays
% the simplest one.
%
% The flag has to reach the compiler through mex, and drivers disagree about
% how: MATLAB wants "CFLAGS=$CFLAGS ..." and rejects a bare compiler option,
% while others take the bare option and reject the assignment. Both spellings
% are tried. There is deliberately NO rung that simply drops the flag -- a build
% without it links and then segfaults after about a dozen interleaved calls, and
% shipping something that crashes later is worse than failing here.
if ispc
    visVariants = {{}};
    fFlags = '-O2 -std=legacy';
else
    visVariants = {{'CFLAGS=$CFLAGS -fvisibility=hidden'}, {'-fvisibility=hidden'}};
    fFlags = '-O2 -std=legacy -fPIC -fvisibility=hidden';
end

% mex classifies its inputs by extension and documents object and library files
% as .obj/.lib on Windows and .o/.a elsewhere. gfortran will happily write
% either, and MinGW's linker reads an archive whatever it is called, so name
% them the way mex expects rather than hoping it passes an unknown extension
% through to the link line.
if ispc, objExt = '.obj'; libExt = '.lib'; else, objExt = '.o'; libExt = '.a'; end

report = struct('name',{},'status',{},'detail',{});

%-------------------------------------------------------------------------------
%% 1. C sources, with the shim prepended
%-------------------------------------------------------------------------------
% Prepending #include rather than passing -include keeps this independent of
% how the selected mex configuration spells its compiler-flag variable.
fprintf('\n--- 1. C programs ---\n');
cRoutines = dir(fullfile(src,'source_c','routines','*.c'));
routineGen = cell(1,numel(cRoutines));
for i = 1:numel(cRoutines)
    routineGen{i} = shimCopy(fullfile(cRoutines(i).folder,cRoutines(i).name), ...
                             build,mexd);
end

cProgs = {'d2','nstat_z','false_nearest','boxcount','lyap_r','poincare'};
for i = 1:numel(cProgs)
    name = cProgs{i};
    progGen = shimCopy(fullfile(src,'source_c',[name '.c']),build,mexd);
    % source_c/routines is on the include path because the support routines say
    % #include "tsa.h" and #include "tisean_cec.h", and a quoted include is
    % resolved relative to the directory of the file doing the including. The
    % originals sit next to those headers; the shimmed copies compiled here do
    % not, so without this they fail to find them.
    base = [{'-output',fullfile(absOut,['tisean_' name]), ...
             ['-I' fullfile(src,'source_c')], ...
             ['-I' fullfile(src,'source_c','routines')], ['-I' mexd], ...
             ['-Dmain=tisean_' name '_main']}, ...
            {fullfile(mexd,'tisean_mex.c'), fullfile(mexd,['tisean_' name '.c']), ...
             progGen}, routineGen];
    v = cell(1,numel(visVariants)); lbl = cell(1,numel(visVariants));
    for k = 1:numel(visVariants)
        v{k} = [visVariants{k} base];  lbl{k} = '';
    end
    report(end+1) = tryBuild(name,absOut,v,lbl); %#ok<AGROW>
end

%-------------------------------------------------------------------------------
%% 2. Fortran sources, rewritten as subroutines
%-------------------------------------------------------------------------------
fProgs = {'c1','c2d','c2g','c2t'};
if isempty(gfortran)
    for i = 1:numel(fProgs)
        report(end+1) = mk(fProgs{i},'skipped','no usable gfortran'); %#ok<AGROW>
    end
else
    fprintf('\n--- 2. Fortran support library ---\n');
    [libF,libS,fobj,ok] = buildFortranLibs(gfortran,src,mexd,build,fFlags,objExt,libExt);
    if ~ok
        for i = 1:numel(fProgs)
            report(end+1) = mk(fProgs{i},'failed','support library did not build'); %#ok<AGROW>
        end
    else
        [gfLibs,libWhy] = gfortranLibs(gfortran,build,libExt);
        if isempty(gfLibs)
            fprintf(2,'  could not locate libgfortran (%s)\n',libWhy);
        end
        fprintf('\n--- 3. Fortran programs ---\n');
        for i = 1:numel(fProgs)
            name = fProgs{i};
            gen = genFortran(fullfile(src,'source_f',[name '.f']),build,name);
            obj = fullfile(build,[name objExt]);
            [st,msg] = runTool(fileparts(gfortran), ...
                sprintf('"%s" %s -c "%s" -o "%s"',gfortran,fFlags,gen,obj));
            if st ~= 0
                fprintf(2,'  %-14s FAILED (gfortran)\n',name);
                showError(msg);
                report(end+1) = mk(name,'failed',firstLine(msg)); %#ok<AGROW>
                continue
            end
            base = [{'-output',fullfile(absOut,['tisean_' name]), ...
                     ['-I' fullfile(src,'source_c')], ...
                     ['-I' fullfile(src,'source_c','routines')], ['-I' mexd]}, ...
                    {fullfile(mexd,'tisean_mex.c'), fullfile(mexd,['tisean_' name '.c']), ...
                     obj, fobj, libF, libS}];
            % Static libgfortran keeps the MEX file self-contained -- the same
            % reasoning as -static in the executable build. Some ELF
            % distributions ship a non-PIC libgfortran.a that cannot go into a
            % shared object at all, so linking it dynamically is the last rung.
            v = {}; lbl = {};
            for k = 1:numel(visVariants)
                v{end+1} = [visVariants{k} base gfLibs];  lbl{end+1} = ''; %#ok<AGROW>
            end
            if ~ispc
                % Some ELF distributions ship a non-PIC libgfortran.a that
                % cannot go into a shared object at all. Linking it dynamically
                % is a real loss of self-containment, so it is reported.
                for k = 1:numel(visVariants)
                    v{end+1} = [visVariants{k} base {'-lgfortran'}]; %#ok<AGROW>
                    lbl{end+1} = 'against the shared libgfortran: it must be present at run time'; %#ok<AGROW>
                end
            end
            report(end+1) = tryBuild(name,absOut,v,lbl); %#ok<AGROW>
        end
    end
end

%-------------------------------------------------------------------------------
%% 4. Smoke test
%-------------------------------------------------------------------------------
fprintf('\n--- 4. smoke test ---\n');
addpath(absOut);
smoke(report,absOut);

nBuilt = sum(strcmp({report.status},'built'));
fprintf('\n===============================================\n');
fprintf('  built %u of 10 TISEAN MEX files into %s\n',nBuilt,outDir);
if nBuilt == 10
    fprintf('  Re-run startup, then HCTSA2_CheckInstall to confirm.\n');
    fprintf('  HCTSA2_TestTiseanMex checks them against the executables.\n');
elseif nBuilt >= 6
    fprintf(2,'  The Fortran programs are missing: NL_TISEAN_c1, NL_TISEAN_d2\n');
    fprintf(2,'  and NL_TakensEstimator will not run. R2024a or later is\n');
    fprintf(2,'  required for gfortran MEX support on Windows.\n');
end
fprintf('===============================================\n\n');
end

%-------------------------------------------------------------------------------
function out = shimCopy(srcFile,build,mexd)
% Copy a TISEAN C source into the build directory with the MEX shim prepended.
% #line keeps any compiler diagnostic pointing at the original file.
[~,base,ext] = fileparts(srcFile);
out = fullfile(build,[base ext]);
txt = fileread(srcFile);
hdr = sprintf('#include "%s"\n#line 1 "%s"\n', ...
              strrep(fullfile(mexd,'tisean_mex.h'),'\','/'), ...
              strrep(srcFile,'\','/'));
fid = fopen(out,'w');
fwrite(fid,[hdr txt]);
fclose(fid);
end

%-------------------------------------------------------------------------------
function out = genFortran(srcFile,build,prog)
% Rewrite a TISEAN Fortran source so it can be called more than once, in
% process, without ending MATLAB.
%
%   * An unnamed main program becomes `subroutine <prog>_main`. gfortran
%     compiles an unnamed main program to a LOCAL MAIN__ plus a global main,
%     so there is nothing callable from outside until this is done.
%   * `stop` in the main program becomes `return`; anywhere else it becomes
%     `call tsfatal(1)`, which longjmps back to the gateway.
%   * `stop "message"` becomes `call tsdie("message")` so the diagnostic
%     survives. The message is often on a continuation line, which is why the
%     opening paren and the closing paren are placed on different lines.
%
% PROG empty means a shared helper (readfile, xreadfile, d1, normal): no
% subroutine header, and every `stop` becomes fatal rather than a return.
[~,base,ext] = fileparts(srcFile);
out = fullfile(build,[base ext]);
L = regexp(fileread(srcFile),'\r\n|\n|\r','split');
if ~isempty(prog)
    L = [{sprintf('      subroutine %s_main',prog)}, L];
end
inMain = ~isempty(prog);
changed = [];
i = 1;
while i <= numel(L)
    l = L{i};
    if ~isempty(regexpi(l,'^\s{6,}subroutine\b','once')) && ~(inMain && i == 1)
        inMain = false;
    end
    if isComment(l) || isCont(l), i = i + 1; continue; end
    % last continuation line of this statement
    j = i; k = i + 1;
    while k <= numel(L)
        if isComment(L{k}), k = k + 1; continue; end
        if isCont(L{k}), j = k; k = k + 1; continue; end
        break
    end
    if ~isempty(regexp(l,'\<stop[ \t]+("[^"]*"|''[^'']*'')[ \t]*$','once'))
        L{i} = regexprep(l,'\<stop[ \t]+("[^"]*"|''[^'']*'')[ \t]*$','call tsdie($1)');
        changed(end+1) = i; %#ok<AGROW>
    elseif ~isempty(regexp(l,'\<stop[ \t]*$','once'))
        if j > i                       % message continues on the next line
            L{i} = regexprep(l,'\<stop[ \t]*$','call tsdie(');
            L{j} = [deblank(L{j}) ')'];
            changed(end+1:end+2) = [i j]; %#ok<AGROW>
        elseif inMain
            L{i} = regexprep(l,'\<stop[ \t]*$','return');
            changed(end+1) = i; %#ok<AGROW>
        else
            L{i} = regexprep(l,'\<stop[ \t]*$','call tsfatal(1)');
            changed(end+1) = i; %#ok<AGROW>
        end
    end
    i = j + 1;
end
% Fixed-form Fortran ignores anything past column 72. A rewritten line that
% overran would lose its tail silently, so fail loudly instead. Only lines this
% function touched are checked -- several TISEAN sources carry trailing "!"
% comments past column 72 already, which the compiler discards as intended.
overlong = changed(cellfun(@(x) numel(x) > 72, L(changed)));
if ~isempty(overlong)
    error('HCTSA2_BuildTiseanMex:column72', ...
        '%s: rewritten line(s) %s exceed column 72.',srcFile,mat2str(overlong));
end
fid = fopen(out,'w');
fprintf(fid,'%s\n',L{:});
fclose(fid);
end

function t = isComment(l)
t = isempty(strtrim(l)) || any(l(1) == 'cC*!');
end

function t = isCont(l)
t = numel(l) > 5 && ~isComment(l) && l(6) ~= ' ' && l(6) ~= '0';
end

%-------------------------------------------------------------------------------
function [libF,libS,fobj,ok] = buildFortranLibs(gfortran,src,mexd,build,fFlags,objExt,libExt)
% libtsa (with the rewritten helpers), libsla, and tsreset.o.
ok = false; libF = ''; libS = ''; fobj = '';
binDir = fileparts(gfortran);
[ar,ranlib] = findAr(binDir);
if isempty(ar)
    fprintf(2,'  ar not found next to gfortran; cannot build static libraries\n');
    return
end

% istdio.f is what configure would generate: the Fortran unit for stderr.
tmpl = fullfile(src,'source_f','istdio_temp.f');
gen  = fullfile(build,'istdio.f');
fid = fopen(gen,'w'); fwrite(fid,strrep(fileread(tmpl),'ERRUNIT','0')); fclose(fid);

plain   = {'arguments','commandline','any_s','help','verbose','neigh', ...
           'rank','nmore','store_spec','tospec'};
rewrite = {'readfile','xreadfile','d1','normal'};   % these contain `stop`
objs = {gen};
for i = 1:numel(plain)
    f = fullfile(src,'source_f',[plain{i} '.f']);
    if isfile(f), objs{end+1} = f; end %#ok<AGROW>
end
for i = 1:numel(rewrite)
    objs{end+1} = genFortran(fullfile(src,'source_f',[rewrite{i} '.f']),build,''); %#ok<AGROW>
end

good = {};
for i = 1:numel(objs)
    [~,b] = fileparts(objs{i});
    o = fullfile(build,['tsa_' b objExt]);
    if compileF(gfortran,fFlags,objs{i},o), good{end+1} = o; end %#ok<AGROW>
end
fprintf('  compiled %u of %u libtsa sources\n',numel(good),numel(objs));

slaSrc = dir(fullfile(src,'source_f','slatec','*.f'));
sgood = {};
for i = 1:numel(slaSrc)
    o = fullfile(build,['sla_' slaSrc(i).name(1:end-2) objExt]);
    if compileF(gfortran,fFlags,fullfile(slaSrc(i).folder,slaSrc(i).name),o)
        sgood{end+1} = o; %#ok<AGROW>
    end
end
fprintf('  compiled %u of %u SLATEC sources\n',numel(sgood),numel(slaSrc));

fobj = fullfile(build,['tsreset' objExt]);
if ~compileF(gfortran,fFlags,fullfile(mexd,'tsreset.f'),fobj)
    fprintf(2,'  tsreset.f did not compile\n'); return
end

libF = fullfile(build,['libtsa' libExt]);
libS = fullfile(build,['libsla' libExt]);
archive(ar,ranlib,libF,good);
archive(ar,ranlib,libS,sgood);
ok = isfile(libF) && isfile(libS) && ~isempty(good) && ~isempty(sgood);
end

%-------------------------------------------------------------------------------
function [p,why] = findGfortran(cc)
% MathWorks' MinGW-w64 add-on ships gfortran, and R2024a is the first release
% in which it is a supported Fortran compiler for MEX on Windows. Earlier
% releases can still run the C half of this build.
p = ''; why = '';
if ispc && verLessThan('matlab','24.1')
    why = sprintf('MATLAB %s predates R2024a, which added gfortran MEX support on Windows',version('-release'));
    return
end
root = '';
if ~isempty(cc), root = cc.Location; end
if isempty(root), root = getenv('MW_MINGW64_LOC'); end
if ~isempty(root)
    c = fullfile(root,'bin',['gfortran' exeExt()]);
    if isfile(c), p = c; return; end
end
p = onPath('gfortran');
if isempty(p)
    why = 'gfortran not found beside the C compiler or on PATH';
end
end

function [ar,ranlib] = findAr(binDir)
ar = fullfile(binDir,['ar' exeExt()]);
if ~isfile(ar), ar = onPath('ar'); end
ranlib = fullfile(binDir,['ranlib' exeExt()]);
if ~isfile(ranlib), ranlib = onPath('ranlib'); end
end

function [libs,why] = gfortranLibs(gfortran,build,libExt)
% Ask the compiler where its own runtime lives and hand the archives to mex by
% full path. Linking libgfortran.a rather than -lgfortran keeps the MEX files
% self-contained -- the same reasoning as -static in the executable build, and
% it means MinGW does not have to be on PATH to run them. Copied into the build
% directory under the extension mex expects.
libs = {}; why = '';
for want = {'libgfortran.a','libquadmath.a'}
    [st,out] = runTool(fileparts(gfortran), ...
        sprintf('"%s" -print-file-name=%s',gfortran,want{1}));
    out = strtrim(out);
    if st == 0 && isfile(out)
        dest = fullfile(build,[want{1}(1:end-2) libExt]);
        if strcmp(libExt,'.a')
            dest = out;
        else
            copyfile(out,dest);
        end
        libs{end+1} = dest; %#ok<AGROW>
    elseif strcmp(want{1},'libgfortran.a')
        why = sprintf('gfortran -print-file-name=%s did not resolve',want{1});
        libs = {}; return
    end
end
end

function e = exeExt()
if ispc, e = '.exe'; else, e = ''; end
end

function p = onPath(tool)
p = '';
if ispc, q = ['where ' tool]; else, q = ['which ' tool]; end
[st,w] = system(q);
if st == 0
    w = strtrim(w);
    nl = strfind(w,newline);
    if ~isempty(nl), w = w(1:nl(1)-1); end
    p = strtrim(w);
end
end

function ok = compileF(fc,fFlags,in,out)
cmd = sprintf('"%s" %s -fallow-argument-mismatch -c "%s" -o "%s"',fc,fFlags,in,out);
[st,~] = runTool(fileparts(fc),cmd);
if st ~= 0
    % -fallow-argument-mismatch is gfortran 10+; MathWorks ships 8.1.0, which
    % rejects the flag outright, so retry without it.
    [st,~] = runTool(fileparts(fc), ...
        sprintf('"%s" %s -c "%s" -o "%s"',fc,fFlags,in,out));
end
ok = (st == 0);
end

function archive(ar,ranlib,lib,objs)
if isempty(objs), return; end
if isfile(lib), delete(lib); end
% A single command line with sixty-odd paths can exceed the shell's limit.
for i = 1:40:numel(objs)
    chunk = objs(i:min(i+39,numel(objs)));
    q = sprintf('"%s" ',chunk{:});
    [st,msg] = runTool(fileparts(ar),sprintf('"%s" rcs "%s" %s',ar,lib,q));
    if st ~= 0, fprintf(2,'    ar failed: %s\n',firstLine(msg)); end
end
if ~isempty(ranlib)
    runTool(fileparts(ranlib),sprintf('"%s" "%s"',ranlib,lib));
end
end

function [st,msg] = runTool(binDir,cmd)
% Invoke a MinGW tool with its own bin directory FIRST on PATH.
%
% gcc.exe and gfortran.exe are drivers: they launch f951.exe, as.exe and
% collect2.exe from the same installation and find them via PATH. MATLAB does
% not put the support-package MinGW on PATH, so calling by absolute path alone
% makes them fail with an EMPTY error message. MATLAB's own bin directory also
% carries a libstdc++/libgcc that can shadow MinGW's; prepending wins both.
old = getenv('PATH');
c = onCleanup(@() setenv('PATH',old)); %#ok<NASGU>
setenv('PATH',[binDir pathsep old]);
[st,msg] = system([cmd ' 2>&1']);
end

%-------------------------------------------------------------------------------
function r = tryBuild(name,absOut,variants,labels)
% Try each link variant in order and take the first that produces a file.
%
% On Windows there is exactly one variant, so this is a plain attempt. Off
% Windows the later rungs drop things that are wanted but not essential, and
% each one that is used is reported -- a build that quietly gave up a guarantee
% is worse than one that failed.
%
% On failure the compiler's whole message is printed. An earlier version
% truncated it to 90 characters, which turned six different failures into six
% identical-looking lines and cut off the one naming the missing header.
msg = '';
for k = 1:numel(variants)
    [ok,msg] = attempt(name,variants{k},absOut);
    if ok
        if isempty(labels{k})
            fprintf('  %-14s built\n',name);
            r = mk(name,'built','');
        else
            fprintf(2,'  %-14s built, %s\n',name,labels{k});
            r = mk(name,'built',labels{k});
        end
        return
    end
end
fprintf(2,'  %-14s FAILED\n',name);
showError(msg);
r = mk(name,'failed',firstLine(msg));
end

function [ok,msg] = attempt(name,args,absOut)
% isfile, not exist(...,'file')==2: the 'file' qualifier limits what exist
% searches, but the value it returns still classifies by type, and for a MEX
% file that is 3. An == 2 test could never be true, so every target built
% successfully and was then reported as having produced no file.
target = fullfile(absOut,['tisean_' name '.' mexext]);
if isfile(target), delete(target); end
try
    mex(args{:});
    ok = isfile(target);
    msg = '';
    if ~ok, msg = 'mex reported success but produced no file'; end
catch err
    ok = false; msg = err.message;
end
end


function showError(msg)
lines = regexp(strtrim(msg),'\r\n|\n|\r','split');
for i = 1:numel(lines)
    if ~isempty(strtrim(lines{i}))
        fprintf(2,'      %s\n',lines{i});
    end
end
end

%-------------------------------------------------------------------------------
function smoke(report,absOut)
% Run each MEX file once on a real series. A build that links but cannot run --
% a missing libgfortran, an export that did not take -- would otherwise only
% show up much later as a NaN in a feature vector.
%
% Each program names the ONE file it is expected to produce, and that file is
% checked and removed before the next program runs. Checking "did any output
% appear" would pass on a leftover: d2's .c2 has to survive for c2g and c2t to
% have an input, so a later program producing nothing would look fine.
if exist('BF_WriteTempFile','file') ~= 2
    fprintf(2,'  skipped: run startup first (BF_WriteTempFile is not on the path)\n');
    return
end
addpath(absOut);
rng(1);
y = cumsum(randn(500,1)); y = (y-mean(y))/std(y);
fn = BF_WriteTempFile(y);
capture = [tempname '.err'];
out = [fn '.out'];
built = {report(strcmp({report.status},'built')).name};

% {program, arguments, the file it must produce}
% Order matters: d2 leaves the .c2 that c2g and c2t read, c1 leaves the .c1
% that c2d reads.
cases = {
 'd2',            {'-d1','-M1,3','-t10','-o',fn,fn},                        [fn '.c2']
 'nstat_z',       {'-#','3','-d1','-m2','-o',out,fn},                       out
 'false_nearest', {'-d1','-m1','-M1,4','-t10','-V0','-o',out,fn},           out
 'boxcount',      {'-M1,4','-d1','-Q2.0','-#10','-o',out,fn},               out
 'lyap_r',        {'-d1','-m3','-t10','-s20','-o',out,fn},                  out
 'poincare',      {'-d1','-m3','-q1','-C0','-o',out,fn},                    out
 'c1',            {'-d1','-m1','-M3','-t10','-n50','-o',[fn '.c1'],fn},     [fn '.c1']
 'c2d',           {'-a2','-o',out,[fn '.c1']},                              out
 'c2g',           {'-o',out,[fn '.c2']},                                    out
 'c2t',           {'-o',out,[fn '.c2']},                                    out
};

keep = {[fn '.c2'],[fn '.c1']};      % inputs for later programs
for i = 1:size(cases,1)
    n = cases{i,1};
    if ~any(strcmp(n,built)), continue; end
    want = cases{i,3};
    if isfile(want) && ~any(strcmp(want,keep)), delete(want); end
    why = '';
    try
        [st,msg] = feval(['tisean_' n],cases{i,2},capture);
        good = (st == 0) && isfile(want);
        if ~good
            why = sprintf(' (status %d)',st);
            if ~isempty(strtrim(msg)), why = [why ' ' firstLine(msg)]; end
        end
    catch err
        good = false; why = [' ' firstLine(err.message)];
    end
    fprintf('  %-14s %s%s\n',n,ternS(good,'runs','DID NOT RUN'),why);
    if isfile(out), delete(out); end
end
for e = {'','.c1','.c2','.d2','.h2','.stat','.out'}
    if isfile([fn e{1}]), delete([fn e{1}]); end
end
if isfile(capture), delete(capture); end
end

%-------------------------------------------------------------------------------
function s = mk(name,status,detail)
s = struct('name',name,'status',status,'detail',detail);
end

function s = firstLine(m)
s = strtrim(regexprep(m,'\n.*','','once'));
if numel(s) > 90, s = [s(1:87) '...']; end
end

function s = ternS(c,a,b)
if c, s = a; else, s = b; end
end
