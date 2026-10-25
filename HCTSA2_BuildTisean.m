function report = HCTSA2_BuildTisean(outDir)
% HCTSA2_BuildTisean   Build the TISEAN programs hctsa needs, from source.
%
%   HCTSA2_BuildTisean
%   report = HCTSA2_BuildTisean('Toolboxes/Tisean_bin')
%
% Compiles the eleven TISEAN executables hctsa calls, and ripser, using the
% compilers that come with MATLAB's MinGW-w64 add-on. No autotools, no shell, no
% prebuilt binaries. Run from the hctsa2 root (the folder holding startup.m).
%
%   d2, nstat_z, false_nearest,  C        -> gcc
%   boxcount, lyap_r,
%   lyap_spec, nrlazy
%   c1, c2d, c2g, c2t            Fortran  -> gfortran
%   ripser                       C++      -> g++   (Toolboxes/ripser, not TISEAN;
%                                                   NL_PersistentHomology)
%
% The bundled false_nearest.c carries the official TISEAN 3.0.1 patch of
% 2009-03-08: the univariate branch built its delay embedding as vemb[i]=i,
% omitting the delay, while the multivariate branch alongside it correctly used
% (i/comp)*delay. hctsa always calls it univariate (-m1) with -d tau, so every
% call with tau > 1 embedded at the wrong lags. See KNOWN_ISSUES.md.
%
% The random number generators in d2.c and routines/rand.c (used by lyap_spec)
% are 64-bit on every platform. The originals chose their generator by
% sizeof(long), which is 4 on Windows, so a Windows build gave different d2 and
% lyap_spec output from Linux/macOS. With the change the Windows programs
% reproduce the Linux output; on Linux nothing changes.
%
% Run HCTSA2_ProbeToolchain first: four of the eleven are Fortran, and whether
% MathWorks' MinGW package includes gfortran decides whether they can be built.
%
%---WHAT THIS REPLACES
% TISEAN normally builds through ./configure && make, which needs a POSIX shell.
% Everything configure actually does for these programs is small enough to do
% directly:
%
%   * istdio.f is generated from istdio_temp.f by substituting the Fortran unit
%     number for stderr. gfortran uses unit 0.
%   * libtsa.a is built from a fixed list of utility sources.
%   * libsla.a is built from the bundled SLATEC routines (c2g needs dqk15).
%   * routines/libddtsa.a is built from the C support routines.
%
%---OUTPUT
% Executables are written to OUTDIR (default Toolboxes/Tisean_bin), which is
% exactly where startup.m looks for them. Each is smoke-tested after building.

if nargin < 1 || isempty(outDir)
    outDir = fullfile('Toolboxes','Tisean_bin');
end

root = pwd;
src = fullfile(root,'Toolboxes','Tisean_3.0.1');
if exist(src,'dir') ~= 7
    error('HCTSA2_BuildTisean:noSource', ...
        'TISEAN source not found at %s. Run this from the hctsa2 root (the folder holding startup.m).',src);
end
if exist(outDir,'dir') ~= 7, mkdir(outDir); end

fprintf('\n=========== building TISEAN from source ===========\n');
fprintf('source : %s\n',src);
fprintf('output : %s\n',outDir);

%-------------------------------------------------------------------------------
%% Locate the toolchain
%-------------------------------------------------------------------------------
[gcc,gfortran,ar,ranlib,gxx] = findTools();
fprintf('\ngcc      : %s\n',emptyAs(gcc,'NOT FOUND'));
fprintf('g++      : %s\n',emptyAs(gxx,'NOT FOUND'));
fprintf('gfortran : %s\n',emptyAs(gfortran,'NOT FOUND'));
fprintf('ar       : %s\n',emptyAs(ar,'NOT FOUND'));
if isempty(gcc)
    error('HCTSA2_BuildTisean:noGcc','No C compiler found. Run "mex -setup C".');
end
if isempty(ar)
    error('HCTSA2_BuildTisean:noAr','ar not found; cannot build static libraries.');
end
if isempty(gfortran)
    fprintf(2,'\ngfortran is missing: c1, c2d, c2g and c2t cannot be built.\n');
    fprintf(2,'The seven C programs will still be produced.\n');
end

binDir = fileparts(gcc);          % MinGW bin: needed on PATH for gcc to run
% (HCTSA2 2.2.1: the scratch name and the smoke test below draw random numbers
% and re-seed the global stream; put it back as the caller left it.)
rngState = rng;
restoreRng = onCleanup(@() rng(rngState)); %#ok<NASGU>
build = fullfile(tempdir,sprintf('tisean_build_%06u',randi(999999)));
mkdir(build);
cleanupObj = onCleanup(@() rmdir(build,'s'));
fprintf('scratch  : %s\n',build);

report = struct('name',{},'status',{},'detail',{});

%-------------------------------------------------------------------------------
%% 1. C support library
%-------------------------------------------------------------------------------
fprintf('\n--- 1. C support routines -> libddtsa.a ---\n');
cRoutines = dir(fullfile(src,'source_c','routines','*.c'));
objs = {};
nFail = 0;
for i = 1:numel(cRoutines)
    in = fullfile(cRoutines(i).folder,cRoutines(i).name);
    out = fullfile(build,[cRoutines(i).name(1:end-2) '.o']);
    cmd = sprintf('"%s" -O2 -c "%s" -o "%s" -I"%s"',gcc,in,out, ...
        fullfile(src,'source_c'));
    [st,msg] = runTool(binDir,cmd);
    if st == 0
        objs{end+1} = out; %#ok<AGROW>
    else
        nFail = nFail + 1;
        if nFail <= 3, fprintf(2,'    %s: %s\n',cRoutines(i).name,firstLine(msg)); end
    end
end
fprintf('  compiled %u of %u C routines\n',numel(objs),numel(cRoutines));
libC = fullfile(build,'libddtsa.a');
archive(ar,ranlib,libC,objs);

%-------------------------------------------------------------------------------
%% 2. Fortran support libraries
%-------------------------------------------------------------------------------
libF = ''; libS = '';
if ~isempty(gfortran)
    fprintf('\n--- 2. Fortran support routines -> libtsa.a, libsla.a ---\n');

    % istdio.f: configure substitutes the stderr unit number. gfortran uses 0.
    tmpl = fullfile(src,'source_f','istdio_temp.f');
    gen  = fullfile(build,'istdio.f');
    t = fileread(tmpl);
    t = strrep(t,'ERRUNIT','0');
    fid = fopen(gen,'w'); fwrite(fid,t); fclose(fid);
    fprintf('  generated istdio.f (stderr unit 0)\n');

    incNames = {'readfile','xreadfile','arguments','commandline','any_s', ...
                'help','verbose','d1','neigh','normal','rank','nmore', ...
                'store_spec','tospec'};
    fobjs = {};
    o = fullfile(build,'istdio.o');
    if compileF(gfortran,gen,o), fobjs{end+1} = o; end
    for i = 1:numel(incNames)
        in = fullfile(src,'source_f',[incNames{i} '.f']);
        if exist(in,'file') ~= 2, continue; end
        o = fullfile(build,[incNames{i} '.o']);
        if compileF(gfortran,in,o), fobjs{end+1} = o; end %#ok<AGROW>
    end
    fprintf('  compiled %u of %u libtsa sources\n',numel(fobjs),numel(incNames)+1);
    libF = fullfile(build,'libtsa.a');
    archive(ar,ranlib,libF,fobjs);

    slaSrc = dir(fullfile(src,'source_f','slatec','*.f'));
    sobjs = {};
    for i = 1:numel(slaSrc)
        in = fullfile(slaSrc(i).folder,slaSrc(i).name);
        o = fullfile(build,['sla_' slaSrc(i).name(1:end-2) '.o']);
        if compileF(gfortran,in,o), sobjs{end+1} = o; end %#ok<AGROW>
    end
    fprintf('  compiled %u of %u SLATEC sources\n',numel(sobjs),numel(slaSrc));
    libS = fullfile(build,'libsla.a');
    archive(ar,ranlib,libS,sobjs);
end

%-------------------------------------------------------------------------------
%% 3. The programs
%-------------------------------------------------------------------------------
fprintf('\n--- 3. programs ---\n');
% Seven programs, not six. NL_FNN (then NL_TISEAN_fnn) builds its command into a
% variable before calling system(), so an earlier enumeration that scanned for
% system(sprintf('<name> ... missed false_nearest entirely. The reliable check
% is to test every TISEAN program name against the Operations source, which is
% how this list was rebuilt.
% Ten programs. The last three are needed by upstream v2.0.0's TSTOOL
% replacements -- NL_BoxCorrDim (NL_BoxCountEntropyRate since upstream fc43329)
% and NL_Dimensions call boxcount, NL_LargestLyap calls lyap_r, NL_PoincareSection
% calls poincare. All three are C, so they add no toolchain requirement beyond
% what was already needed.
% Twelve programs: upstream 3.0.0-alpha added NL_LyapSpec (lyap_spec) and
% PP_SchreiberDenoise (nrlazy), both C.
% Eleven (HCTSA2 2.5.0): upstream f54e4d1 computes NL_PoincareSection's section in
% MATLAB, so poincare is no longer called and is not built (its source stays).
progs = {
 'd2',            'C'
 'nstat_z',       'C'
 'false_nearest', 'C'
 'boxcount',      'C'
 'lyap_r',        'C'
 'lyap_spec',     'C'
 'nrlazy',        'C'
 'c1',            'F'
 'c2d',           'F'
 'c2g',           'F'
 'c2t',           'F'
 'ripser',        'CXX'
};
for i = 1:size(progs,1)
    name = progs{i,1}; lang = progs{i,2};
    exe = fullfile(root,outDir,[name '.exe']);
    if strcmp(lang,'CXX')
        % ripser: upstream's compile_ripser flags, plus -static as below. It is
        % written to the same folder as the TISEAN programs, where BF_RunTisean
        % looks first.
        if isempty(gxx)
            fprintf('  %-9s SKIPPED (no g++)\n',name);
            report(end+1) = mk(name,'skipped','no g++'); %#ok<AGROW>
            continue
        end
        in = fullfile(root,'Toolboxes','ripser','ripser.cpp');
        cmd = sprintf('"%s" -std=c++11 -O3 -D NDEBUG -static "%s" -o "%s"', ...
            gxx,in,exe);
    elseif strcmp(lang,'C')
        in = fullfile(src,'source_c',[name '.c']);
        % -static: link the MinGW runtime in. Without it the executables need
        % libgcc_s_seh-1.dll (and for Fortran, libgfortran-5.dll,
        % libquadmath-0.dll, libwinpthread-1.dll) beside them or on PATH.
        % Verified: a dynamically linked build depends on libgfortran and
        % libgcc_s; the static one has no such dependencies and behaves
        % identically. Self-contained matters more than size here -- an offline
        % machine should not need MinGW on PATH to run these.
        cmd = sprintf('"%s" -O2 -static "%s" -o "%s" -I"%s" "%s" -lm', ...
            gcc,in,exe,fullfile(src,'source_c'),libC);
    else
        if isempty(gfortran)
            fprintf('  %-9s SKIPPED (no gfortran)\n',name);
            report(end+1) = mk(name,'skipped','no gfortran'); %#ok<AGROW>
            continue
        end
        in = fullfile(src,'source_f',[name '.f']);
        cmd = sprintf('"%s" -O2 -static -std=legacy -fallow-argument-mismatch "%s" -o "%s" "%s" "%s"', ...
            gfortran,in,exe,libF,libS);
    end
    % (HCTSA2 2.2.1: remove the previous executable first. A failed build used to
    % leave the OLD one here, where HCTSA2_CheckInstall counts it and BF_RunTisean
    % runs it -- an old lyap_spec.exe ignores 2.2.0's -t option and gives
    % different numbers. On Windows it cannot be removed while a pool or another
    % MATLAB is running it.)
    if exist(exe,'file') == 2
        wsDel = warning('off','MATLAB:DELETE:Permission');
        try delete(exe); catch, end %#ok<CTCH,NOSEM>
        warning(wsDel);
        if exist(exe,'file') == 2
            fprintf(2,'  %-9s FAILED: cannot replace the old %s.exe (still running?)\n',name,name);
            report(end+1) = mk(name,'failed','old executable in use, not replaced'); %#ok<AGROW>
            continue
        end
    end
    [st,msg] = runTool(binDir,cmd);
    if st ~= 0 && strcmp(lang,'F')
        % -fallow-argument-mismatch is gfortran 10+. MathWorks ships 8.1.0, which
        % rejects the flag outright, so retry without it.
        [st,msg] = runTool(binDir,strrep(cmd,' -fallow-argument-mismatch',''));
    end
    if st == 0 && exist(exe,'file') == 2
        fprintf('  %-9s built\n',name);
        report(end+1) = mk(name,'built',''); %#ok<AGROW>
    else
        fprintf(2,'  %-9s FAILED: %s\n',name,firstLine(msg));
        report(end+1) = mk(name,'failed',firstLine(msg)); %#ok<AGROW>
    end
end

%-------------------------------------------------------------------------------
%% 4. Smoke test
%-------------------------------------------------------------------------------
fprintf('\n--- 4. smoke test ---\n');
oldPath = getenv('PATH');
% The executables are statically linked, so MinGW should not be needed here.
% Adding its bin directory anyway means a partially dynamic build still gets a
% fair test rather than failing for a missing DLL.
setenv('PATH',[fullfile(root,outDir) pathsep binDir pathsep oldPath]);
restorePath = onCleanup(@() setenv('PATH',oldPath));
rng(1);
y = cumsum(randn(500,1)); y = (y-mean(y))/std(y);
% (HCTSA2 2.2.1: an hc###### file from BF_WriteTempFile, removed with every file
% the programs write beside it, on every way out. The fixed suffix list this
% replaces missed .box, .ros, .stat and .laz.1.)
if exist('BF_WriteTempFile','file') ~= 2
    addpath(fullfile(root,'PeripheryFunctions'));
end
fn = BF_WriteTempFile(y,7);
cleanSmoke = onCleanup(@() BF_CleanTempFiles(fn)); %#ok<NASGU>
for i = 1:numel(report)
    if ~strcmp(report(i).status,'built'), continue; end
    n = report(i).name;
    switch n
        case 'd2',      c = sprintf('d2 -d1 -M1,3 -t10 "%s"',fn);
        case 'nstat_z', c = sprintf('nstat_z -# 3 -d1 -m2 "%s"',fn);
        case 'false_nearest', c = sprintf('false_nearest -d1 -m1 -M1,4 -t10 -V0 "%s"',fn);
        case 'boxcount', c = sprintf('boxcount -d1 -M1,4 -Q10 "%s"',fn);
        case 'lyap_r',   c = sprintf('lyap_r -d1 -m3 -t10 -s20 "%s"',fn);
        case 'lyap_spec', c = sprintf('lyap_spec -m1,3 -d1 -k30 "%s"',fn);
        case 'nrlazy',   c = sprintf('nrlazy -m1,3 -d1 -i1 -v0.5 "%s"',fn);
        case 'ripser',   c = sprintf('ripser --format point-cloud --dim 1 "%s"',fn);
        case 'c1',      c = sprintf('c1 -d1 -m2 -M3 -t10 -n50 -o "%s.c1" "%s"',fn,fn);
        otherwise,      c = sprintf('%s -h',n);
    end
    [st,out] = system([c ' 2>&1']);
    ok = ~isempty(regexpi(out,'TISEAN','once')) || st == 0;
    fprintf('  %-9s %s\n',n,ternS(ok,'runs','DID NOT RUN'));
    if ~ok, report(i).detail = firstLine(out); end
end

nBuilt = sum(strcmp({report.status},'built'));
fprintf('\n=================================================\n');
fprintf('  built %u of %u programs into %s\n',nBuilt,size(progs,1),outDir);
if nBuilt == size(progs,1)
    fprintf('  Re-run startup, then HCTSA2_CheckInstall to confirm.\n');
else
    % (HCTSA2 2.2.1: name what is actually missing and what that costs. This used
    % to blame the Fortran programs for any shortfall, whatever failed.)
    notBuilt = report(~strcmp({report.status},'built'));
    fprintf(2,'  Not built, so these operations will not run:\n');
    for i = 1:numel(notBuilt)
        fprintf(2,'    %-14s %-8s %s\n',notBuilt(i).name,notBuilt(i).status,usedBy(notBuilt(i).name));
    end
end
fprintf('=================================================\n\n');
end

%-------------------------------------------------------------------------------
function [gcc,gfortran,ar,ranlib,gxx] = findTools()
root = '';
try
    cc = mex.getCompilerConfigurations('C','Selected');
    if ~isempty(cc), root = cc(1).Location; end
catch
end
if isempty(root), root = getenv('MW_MINGW64_LOC'); end
gcc = pick(root,'gcc'); gfortran = pick(root,'gfortran'); gxx = pick(root,'g++');
ar = pick(root,'ar');   ranlib = pick(root,'ranlib');
end

function p = pick(root,tool)
p = '';
if ~isempty(root)
    c = fullfile(root,'bin',[tool '.exe']);
    if exist(c,'file') == 2, p = c; return; end
end
[st,w] = system(['where ' tool]);
if st == 0
    w = strtrim(w); nl = strfind(w,newline);
    if ~isempty(nl), w = w(1:nl(1)-1); end
    p = strtrim(w);
end
end

function ok = compileF(fc,in,out)
binDir = fileparts(fc);
cmd = sprintf('"%s" -O2 -std=legacy -fallow-argument-mismatch -c "%s" -o "%s"',fc,in,out);
[st,~] = runTool(binDir,cmd);
if st ~= 0
    % gfortran 8.1 (what MathWorks ships) does not know that flag.
    cmd = sprintf('"%s" -O2 -std=legacy -c "%s" -o "%s"',fc,in,out);
    [st,~] = runTool(binDir,cmd);
end
ok = (st == 0);
end

function archive(ar,ranlib,lib,objs)
if isempty(objs), return; end
if exist(lib,'file') == 2, delete(lib); end
% Batch the objects: a single command line with 60+ paths can exceed limits.
for i = 1:40:numel(objs)
    chunk = objs(i:min(i+39,numel(objs)));
    q = sprintf('"%s" ',chunk{:});
    [st,msg] = runTool(fileparts(ar),sprintf('"%s" rcs "%s" %s',ar,lib,q));
    if st ~= 0, fprintf(2,'    ar failed: %s\n',firstLine(msg)); end
end
if ~isempty(ranlib), runTool(fileparts(ranlib),sprintf('"%s" "%s"',ranlib,lib)); end
end


%-------------------------------------------------------------------------------
function [st,msg] = runTool(binDir,cmd)
% Invoke a MinGW tool with its own bin directory FIRST on PATH.
%
% gcc.exe is a driver: it launches cc1.exe, as.exe, ld.exe and collect2.exe from
% the same installation, and it locates them via PATH. MATLAB does not put the
% support-package MinGW on PATH, so calling gcc.exe by absolute path alone makes
% it fail with an EMPTY error message -- the driver cannot start its own
% subprocesses and reports nothing useful.
%
% MATLAB also puts its own bin directory (with its own libstdc++/libgcc) ahead of
% everything else, which can shadow MinGW's runtime. Prepending MinGW wins both.
old = getenv('PATH');
c = onCleanup(@() setenv('PATH',old));
setenv('PATH',[binDir pathsep old]);
[st,msg] = system([cmd ' 2>&1']);
end

function s = usedBy(prog)
% The operations that cannot run without PROG. NL_Dimensions, NL_GPCorrSum and
% NL_LargestLyap also need false_nearest for their 'fnn' embeddings (see the
% table in HCTSA2_CheckInstall).
switch prog
    case 'd2',            s = 'NL_d2, NL_Dimensions, NL_GPCorrSum';
    case 'nstat_z',       s = 'SY_nstat_z';
    case 'false_nearest', s = 'NL_FNN, NL_Dimensions, NL_GPCorrSum, NL_LargestLyap and every operation with an ''fnn'' embedding';
    case 'boxcount',      s = 'NL_BoxCountEntropyRate, NL_Dimensions';
    case 'lyap_r',        s = 'NL_LargestLyap';
    case 'lyap_spec',     s = 'NL_LyapSpec';
    case 'nrlazy',        s = 'PP_SchreiberDenoise';
    case {'c1','c2d'},    s = 'NL_c1';
    case {'c2g','c2t'},   s = 'NL_d2';
    case 'ripser',        s = 'NL_PersistentHomology';
    otherwise,            s = '';
end
end

function s = mk(name,status,detail)
s = struct('name',name,'status',status,'detail',detail);
end

function s = firstLine(m)
s = strtrim(regexprep(m,'\n.*','','once'));
if numel(s) > 80, s = [s(1:77) '...']; end
end

function s = emptyAs(v,alt)
if isempty(v), s = alt; else, s = v; end
end

function s = ternS(c,a,b)
if c, s = a; else, s = b; end
end
