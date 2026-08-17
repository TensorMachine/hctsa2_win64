function HCTSA2_ProbeToolchain
% HCTSA2_ProbeToolchain   Can this machine build TISEAN from source?
%
%   HCTSA2_ProbeToolchain
%
% hctsa calls six TISEAN programs. Four of them are FORTRAN:
%
%   d2, nstat_z        C
%   c1, c2d, c2g, c2t  Fortran 77
%
% MATLAB's MinGW-w64 add-on is documented as a C/C++ compiler. Whether the
% package MathWorks ships also contains gfortran varies by version, and that
% single fact decides whether TISEAN can be built here at all. This checks,
% rather than assuming either way.
%
% Run it and send the whole console output.

fprintf('\n=========== toolchain probe ===========\n');
fprintf('MATLAB %s on %s\n',version,computer);

%-------------------------------------------------------------------------------
%% Where is MinGW?
%-------------------------------------------------------------------------------
fprintf('\n--- MinGW location ---\n');
root = '';
try
    cc = mex.getCompilerConfigurations('C','Selected');
    if ~isempty(cc)
        fprintf('  selected C compiler : %s\n',cc(1).Name);
        fprintf('  location            : %s\n',cc(1).Location);
        root = cc(1).Location;
    end
catch err
    fprintf('  mex.getCompilerConfigurations failed: %s\n',err.message);
end
if isempty(root)
    root = getenv('MW_MINGW64_LOC');
    if ~isempty(root)
        fprintf('  from MW_MINGW64_LOC : %s\n',root);
    end
end
if isempty(root)
    fprintf(2,'  Could not locate MinGW. Run "mex -setup C" first.\n');
end

%-------------------------------------------------------------------------------
%% Which compilers are actually present?
%-------------------------------------------------------------------------------
fprintf('\n--- compilers in that installation ---\n');
tools = {'gcc','g++','gfortran','ar','ranlib','make','mingw32-make'};
found = struct();
for i = 1:numel(tools)
    t = tools{i};
    p = '';
    if ~isempty(root)
        cand = fullfile(root,'bin',[t '.exe']);
        if exist(cand,'file') == 2, p = cand; end
    end
    if isempty(p)
        [st,w] = system(['where ' t]);      % fall back to PATH
        if st == 0
            w = strtrim(w);
            nl = strfind(w,newline);
            if ~isempty(nl), w = w(1:nl(1)-1); end
            p = strtrim(w);
        end
    end
    found.(matlab.lang.makeValidName(t)) = p;
    if isempty(p)
        fprintf('  %-14s NOT FOUND\n',t);
    else
        [~,v] = system(['"' p '" --version']);
        v = strtrim(v); nl = strfind(v,newline);
        if ~isempty(nl), v = v(1:nl(1)-1); end
        fprintf('  %-14s %s\n',t,p);
        fprintf('  %-14s   %s\n','',v);
    end
end

%-------------------------------------------------------------------------------
%% The decisive question
%-------------------------------------------------------------------------------
fprintf('\n--- verdict ---\n');
haveC = ~isempty(found.gcc);
haveF = ~isempty(found.gfortran);
if haveC && haveF
    fprintf('  gcc AND gfortran are both present.\n');
    fprintf('  All six TISEAN programs can be built from source on this machine.\n');
elseif haveC
    fprintf(2,'  gcc is present but gfortran is NOT.\n');
    fprintf(2,'  d2 and nstat_z (C) can be built; c1, c2d, c2g and c2t (Fortran)\n');
    fprintf(2,'  cannot. Those four feed NL_TISEAN_c1 -- worth 132 features.\n');
    fprintf(2,'  NL_TISEAN_d2 and SY_TISEAN_nstat_z would still work.\n');
else
    fprintf(2,'  No C compiler found. Run "mex -setup C" first.\n');
end

% A real compile test is worth more than a version string.
fprintf('\n--- actual compile test ---\n');
td = tempname; mkdir(td);
cu = onCleanup(@() rmdir(td,'s'));
if haveC
    f = fullfile(td,'t.c');
    fid = fopen(f,'w');
    fprintf(fid,'#include <stdio.h>\n#include <math.h>\nint main(void){printf("%%g\\n",sqrt(2.0));return 0;}\n');
    fclose(fid);
    [st,msg] = runTool(fileparts(found.gcc), ...
        sprintf('"%s" "%s" -o "%s" -lm',found.gcc,f,fullfile(td,'t.exe')));
    fprintf('  C   : %s\n',verdict(st,msg));
end
if haveF
    f = fullfile(td,'t.f');
    fid = fopen(f,'w');
    fprintf(fid,'      program t\n      write(*,*) sqrt(2.0)\n      end\n');
    fclose(fid);
    [st,msg] = runTool(fileparts(found.gfortran), ...
        sprintf('"%s" -std=legacy "%s" -o "%s"',found.gfortran,f,fullfile(td,'tf.exe')));
    fprintf('  F77 : %s\n',verdict(st,msg));
end
fprintf('\n======================================\n\n');
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

function s = verdict(st,msg)
if st == 0
    s = 'compiles OK';
else
    m = strtrim(msg);
    nl = strfind(m,newline);
    if ~isempty(nl), m = m(1:nl(1)-1); end
    if isempty(m), m = sprintf('(no output, exit code %d)',st); end
    s = ['FAILED -- ' m];
end
end
