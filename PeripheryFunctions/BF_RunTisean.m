function [status,out] = BF_RunTisean(cmd,keepStderr)
% BF_RunTisean   Run a TISEAN program and return its data.
%
%   [status,out] = BF_RunTisean('nstat_z -# 4 -d1 -m3 /tmp/hc123456')
%   [status,out] = BF_RunTisean('d2 ...',true)   % keep stderr
%
% Runs one of the twelve TISEAN executables, or ripser, and returns its exit
% status and its text. Every operation that calls an external program goes
% through here rather than calling system() directly, so the timeout, the
% not-installed check and the temporary-file cleanup apply uniformly.
% (NL_PersistentHomology's ripser is not part of TISEAN; it shares the runner
% because it needs exactly the same protection.)
%
% KEEPSTDERR (default false) controls whether the program's stderr is appended.
%
%   false -- for programs whose DATA is what the caller parses (nstat_z, c2d,
%            c2g, c2t, false_nearest). Mixing the banner in corrupts the parse.
%   true  -- for programs that write their results to FILES and use stdout for
%            nothing at all (d2, c1, boxcount, lyap_r, poincare). Their only
%            text is the stderr banner, so discarding it returns an empty
%            string even on success, and callers that test isempty(out) as a
%            did-it-run check then fail.
%
%---WHY STDERR IS SEPARATED AT ALL
% MATLAB's system() merges stdout and stderr into one string. TISEAN writes its
% data to stdout and its banner and progress messages to stderr, and the two
% streams are buffered differently -- stderr unbuffered, stdout block-buffered
% when redirected. On Windows this produces two failures that are invisible on
% Linux and macOS:
%
%   1. ORDER. The messages can arrive AFTER the data, so the "Writing to
%      stdout" line that hctsa's parsers use as a marker ends up at the END of
%      the capture. Slicing the text after that marker then discards every data
%      row. This is why nstat_z reported "0 data rows" while its output plainly
%      contained all 16.
%
%   2. SPLIT LINES. A stderr message can land in the middle of a stdout line,
%      breaking one data row into two fragments -- as in a c2g capture that
%      began "-005  5.14077711", the tail of a number severed from its head.
%
% stderr is therefore redirected away unless the caller asks for it, and
% callers parse by content rather than by position relative to a stderr marker,
% which is more robust on every platform.

% ------------------------------------------------------------------------------
% Copyright (C) 2020, Ben D. Fulcher <ben.d.fulcher@gmail.com>
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version. See <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

if nargin < 2 || isempty(keepStderr), keepStderr = false; end

args = localTokenize(cmd);
prog = args{1};
args = args(2:end);

% Anything that leaves this function by raising takes the caller's operation
% down with it, so the operation's own cleanup never runs and its temporary
% files survive. Measured: an aborted d2 left five -- the input file plus .c2,
% .d2, .h2 and .stat. Sweep them here, where the abort actually happens.
%
% Only on the RAISING paths. On a normal return the caller often still needs
% what TISEAN just wrote: NL_d2 runs d2, then feeds the .c2 it produced
% to c2g and then to c2t. Deleting outputs on success would break the operation
% one call later.
try
    checkPathLengths(prog,args);

    % Check the executable is reachable before shelling out to it. Without
    % this the failure is silent and misdiagnosed: system() runs a command
    % that does not exist, the shell writes "not recognized" to stderr, the
    % redirect sends stderr to nul, and BF_RunTisean returns an EMPTY string
    % with no sign that nothing ran. Each operation then invents its own
    % downstream complaint, so a machine with no TISEAN installed produces ten
    % unrelated-looking errors and no clue the cause is common to them.
    exePath = assertProgramAvailable(prog);
    [status,out] = runViaExe(exePath,args,cmd,keepStderr,timeoutSeconds());
catch err
    if ~isempty(args)
        BF_CleanTempFiles(args{end});
    end
    rethrow(err);
end
end

%-------------------------------------------------------------------------------
function exePath = assertProgramAvailable(prog)
% Resolve PROG to an executable, or raise a single actionable error.
%
% Returns the absolute path when the program is where startup puts it, and the
% bare name when it was found somewhere else on the PATH. Callers run whatever
% comes back, so the file that is checked is the file that runs.
%
%---WHY THIS DOES NOT ASK THE SHELL FIRST
% It used to: "where %s" on Windows, "command -v %s" elsewhere. That reported
% the program missing on a machine where it was installed and working, and the
% resulting error fired on every TISEAN call.
%
% The reason the shell probe failed there was never established -- it could not
% be reproduced off Windows -- but it does not need to be, because the probe was
% answering a harder question than the one being asked. The executables are
% BUILT INTO A KNOWN DIRECTORY by HCTSA2_BuildTisean, and that is the same
% directory startup adds to the PATH and the same one HCTSA2_CheckInstall looks
% in with isfile, which reports correctly on the machine where the probe did
% not. Looking straight at the file needs no shell, no PATH lookup, no PATHEXT
% and no process, and it is exact about the extension: d2.exe on Windows, d2
% everywhere else.
%
% The shell probe is kept as a fallback for an installation somewhere else
% entirely, with the extension made explicit on Windows.
%---WHY ONLY A SUCCESS IS REMEMBERED
% A failure is deliberately NOT cached. Caching one makes a recoverable
% situation permanent: probe before the executables are built, or before
% startup has put their directory on the PATH, and every later call in that
% MATLAB session keeps the stale answer even after the cause is fixed. The only
% way out is "clear functions", which nobody thinks to try -- the symptom is a
% program reported missing while it sits in the folder the error names.
% Re-probing a failure costs one isfile() on a path that is about to raise an
% error anyway, which is nothing beside making the error impossible to clear.
persistent resolved
if isempty(resolved), resolved = struct(); end
field = ['p_' prog];
if isfield(resolved,field)
    exePath = resolved.(field);
    return                       % only successes are ever stored
end
exePath = locateProgram(prog);
if ~isempty(exePath)
    resolved.(field) = exePath;
    return
end

error('hctsa:tisean:notInstalled', ...
    ['TISEAN program ''%s'' could not be found, so no TISEAN feature can be ' ...
     'computed on this machine.\n  Looked for: %s\n  and on the PATH as %s.' ...
     '\n  The executables are BUILD ARTIFACTS: they are machine-specific, ' ...
     'they are listed in .gitignore, and they do not travel with the ' ...
     'repository.\n  Build them on this machine with HCTSA2_BuildTisean, then ' ...
     'run startup.\n  HCTSA2_CheckInstall reports what is present.'], ...
     prog,fullfile(tiseanBinDir(),exeName(prog)),exeName(prog));
end

%-------------------------------------------------------------------------------
function p = locateProgram(prog)
% The install directory first, then the rest of the PATH. HCTSA2_BuildTisean
% puts ripser in the same install directory as the TISEAN programs; ripser built
% with upstream's compile_ripser lands in Toolboxes/ripser/bin instead.
p = fullfile(tiseanBinDir(),exeName(prog));
if isfile(p), return; end
if strcmp(prog,'ripser')
    p = fullfile(fileparts(fileparts(mfilename('fullpath'))),'Toolboxes','ripser','bin',exeName(prog));
    if isfile(p), return; end
end
if ispc
    probe = sprintf('where %s >nul 2>&1',exeName(prog));
else
    probe = sprintf('command -v %s >/dev/null 2>&1',exeName(prog));
end
if system(probe) == 0
    p = prog;          % somewhere else on the PATH; let the OS resolve it
else
    p = '';
end
end

%-------------------------------------------------------------------------------
function n = exeName(prog)
if ispc, n = [prog '.exe']; else, n = prog; end
end

%-------------------------------------------------------------------------------
function d = tiseanBinDir()
% Where the executables are built, and where startup points the PATH. This file
% lives in PeripheryFunctions/, so the hctsa root is one level up. Windows:
% Toolboxes/Tisean_bin (HCTSA2_BuildTisean). Elsewhere: upstream's
% Toolboxes/Tisean_3.0.1/bin (compile_tisean) when it exists, else ~/bin.
root = fileparts(fileparts(mfilename('fullpath')));
if ispc
    d = fullfile(root,'Toolboxes','Tisean_bin');
else
    d = fullfile(root,'Toolboxes','Tisean_3.0.1','bin');
    if ~isfolder(d)
        d = fullfile(getenv('HOME'),'bin');
    end
end
end


%-------------------------------------------------------------------------------
function sec = timeoutSeconds()
% Wall-clock budget for one TISEAN run, in seconds. 0 disables the limit.
%
% Read on every call so it can be changed without clearing this function, and
% settable per worker with
%
%   pctRunOnAll setenv('HCTSA2_TISEAN_TIMEOUT','900')
%
% 300 s is far above any healthy call -- the overwhelming majority finish in
% well under a second -- and far below the runs that make hctsa look hung.
%
% It is a real trade-off, not a free win. d2 called with -N0 on a 100k-point
% series needs roughly 25 minutes to finish LEGITIMATELY, and this kills it at
% five. Those features become NaN rather than wrong, which is the right way
% round, but if you want them, raise this.
sec = 300;
v = str2double(strtrim(getenv('HCTSA2_TISEAN_TIMEOUT')));
if ~isnan(v) && v >= 0
    sec = v;
end
end

%-------------------------------------------------------------------------------
function [status,out] = runViaExe(exePath,args,cmd,keepStderr,timeoutSec)
% Run the executable, abandoning it if it outstays TIMEOUTSEC.
%
%---WHY A TIMEOUT IS NEEDED AT ALL
% Neither d2 nor false_nearest can loop forever -- every loop in both is
% bounded, and the length guards that would otherwise let an unsigned
% subtraction wrap are intact. What they can do is take a quadratic number of
% steps with no way to stop early:
%
%   d2            -N0 sets the pair limit to ULONG_MAX (d2.c:121), which makes
%                 the only early exit, "if (imin > howoften1) exit(0)" at
%                 d2.c:579, unreachable. d2 then always runs its full
%                 nmax-iteration loop with an O(n) box scan inside each one.
%                 Measured at -M1,10: 4k 2.3s, 8k 9.5s, 16k 38.4s -- four times
%                 the work for twice the data. Extrapolating, 500k points is
%                 about ten hours. NL_Dimensions passes -N0 outright;
%                 NL_GPCorrSum reaches it through Nref = -1. NL_d2 passes no
%                 -N, gets the default 1000, and is unaffected -- which is why
%                 only some features hang.
%
%   false_nearest has no pair limit at all. A point whose nearest neighbour
%                 lies farther than varianz/rt can never be marked done
%                 (false_nearest.c:196) and exact-duplicate vectors (since
%                 HCTSA2 3.0.0, any within round-off) are rejected outright
%                 (:185), so a constant or near-constant
%                 stretch leaves points that are rescanned at every epsilon
%                 step while all sitting in one search box. Measured at
%                 N = 16000: 4.1s on clean data, 41.8s when half the series is
%                 flat.
%
% Because both terminate eventually, killing the process is sound: nothing is
% deadlocked, the run is merely too expensive to wait for.
%
%---WHY THIS RAISES AN ERROR RATHER THAN RETURNING WHAT IT GOT
% d2 rewrites its .c2/.d2/.h2 files every 120 seconds (d2.c:32), so a killed
% run usually leaves a COMPLETE-LOOKING but partial set of output files on
% disk. Returning them would hand the caller plausible numbers computed from a
% fraction of the data -- silently wrong features, which is worse than none.
% An error becomes a NaN in the feature matrix, which is honest.
%
%---WHY JAVA AND NOT system()
% system() blocks until the child exits; there is no way to take it back.
% MATLAB's own Java (no jar of ours: 2.4.0 removed JIDT) gives both a timed
% wait and a kill. It also takes
% argv as a list, so nothing passes through a shell and destroyForcibly kills
% the TISEAN process itself rather than a shell that owns it. When there is no
% JVM, or Java cannot start the program, this falls back to the old blocking
% system() call rather than failing.
persistent warned
[~,prog] = fileparts(exePath);       % short name, for messages only
if timeoutSec <= 0 || ~usejava('jvm')
    [status,out] = runViaSystem(exePath,cmd,keepStderr);
    return
end

[outFile,errFile] = captureFiles(args);
% (HCTSA2 2.2.1: procRef is filled in once the child exists, so that the one
% cleanup below can kill it first and delete its capture files second. An
% interrupt or an error while waiting used to leave the child running.)
procRef = javaObject('java.util.concurrent.atomic.AtomicReference');
cleanup = onCleanup(@() cleanRun(procRef,{outFile,errFile})); %#ok<NASGU>

try
    jargs = javaArray('java.lang.String',numel(args)+1);
    jargs(1) = javaObject('java.lang.String',exePath);
    for i = 1:numel(args)
        jargs(i+1) = javaObject('java.lang.String',args{i});
    end
    pb = javaObject('java.lang.ProcessBuilder',jargs);
    if isunix && ~ismac
        % MATLAB on Linux puts its own runtime libraries (libstdc++, libgfortran)
        % on LD_LIBRARY_PATH, and they can shadow the ones the programs were
        % built against -- ripser (C++) and the Fortran programs fail to start.
        % Upstream's runner clears it for the child too. Windows and macOS do
        % not use this variable.
        env = pb.environment();      % (no chained call: MATLAB cannot index a
        env.remove('LD_LIBRARY_PATH');  %  method call's result directly)
    end
    pb.redirectOutput(javaObject('java.io.File',outFile));
    if keepStderr
        % Merged into stdout, matching what system() hands back.
        pb.redirectErrorStream(true);
    else
        pb.redirectError(javaObject('java.io.File',errFile));
    end
    proc = pb.start();
    procRef.set(proc);
catch err
    % Java could not launch it at all -- a PATH resolution difference, a
    % missing JVM class, a sandbox. Fall back rather than turn a working
    % configuration into a broken one, but say so once so a silently
    % unprotected run is not mistaken for a protected one.
    if isempty(warned)
        warned = true;
        warning('hctsa:tisean:noTimeout', ...
            ['Could not launch %s through Java, so the %g-second TISEAN ' ...
             'timeout is NOT in force for this session: %s'], ...
            prog,timeoutSec,err.message);
    end
    [status,out] = runViaSystem(exePath,cmd,keepStderr);
    return
end

% (HCTSA2 2.2.1: no program reads standard input here, but when no argument names
% an existing file every TISEAN program falls back to reading it, so a missing
% input file waited out the whole timeout. Closing it makes the program see end
% of file and stop at once.)
try
    procStdin = proc.getOutputStream();
    procStdin.close();
catch
end

finished = proc.waitFor(int64(timeoutSec), ...
    javaMethod('valueOf','java.util.concurrent.TimeUnit','SECONDS'));

if ~finished
    proc.destroyForcibly();
    try proc.waitFor(); catch, end %#ok<CTCH>
    error('hctsa:tisean:timeout', ...
        ['TISEAN program ''%s'' was still running after %g seconds and was ' ...
         'stopped.\n  Command: %s\n  This is cost, not a deadlock: d2 with ' ...
         '-N0 and false_nearest on near-constant data both grow roughly with ' ...
         'the SQUARE of the series length.\n  Shorten the series, or raise ' ...
         'the budget with setenv(''HCTSA2_TISEAN_TIMEOUT'',''1800'') ' ...
         '(pctRunOnAll on a parallel pool).'], ...
        prog,timeoutSec,cmd);
end

status = double(proc.exitValue());
out = readTextFile(outFile);
out = strrep(out,char(13),'');
end

%-------------------------------------------------------------------------------
function cleanRun(procRef,files)
% Kill the child if it is still running, then remove its capture files. A normal
% return and the timeout branch have already finished it, so only an interrupt
% (Ctrl+C, a cancelled parfeval) or an error raised while waiting gets here with
% a live process. Killing first matters on Windows, where a running program
% holds the files open and delete would fail.
try
    p = procRef.get();
    if ~isempty(p) && p.isAlive()
        p.destroyForcibly();
        p.waitFor(int64(5),javaMethod('valueOf','java.util.concurrent.TimeUnit','SECONDS'));
    end
catch
end
cleanFiles(files);
end

%-------------------------------------------------------------------------------
function [outFile,errFile] = captureFiles(args)
% Where the program's stdout and stderr are captured. Named after the hc######
% temporary file in the command (hc123456.tout, hc123456.terr), so that anything
% left behind -- a killed MATLAB never runs the onCleanup -- belongs to that file
% set: BF_CleanTempFiles and TS_Compute's stale-file sweep remove it with the
% rest. Two runs never share a stem at the same time: the stem was claimed
% atomically by BF_WriteTempFile, and one operation runs its programs in turn.
% A command with no hc###### argument falls back to tempname.
outFile = ''; errFile = '';
for i = 1:numel(args)
    [folder,name,ext] = fileparts(args{i});
    tok = regexp([name ext],'^(hc\d{6})($|[._])','tokens','once');
    if ~isempty(tok)
        stem = fullfile(folder,tok{1});
        outFile = [stem '.tout'];
        errFile = [stem '.terr'];
        return
    end
end
outFile = [tempname '.tout'];
errFile = [tempname '.terr'];
end

%-------------------------------------------------------------------------------
function [status,out] = runViaSystem(exePath,cmd,keepStderr)
% Run the RESOLVED executable, not the bare name the caller wrote. This path is
% only reached without a JVM, or when Java could not launch the program, and in
% both cases leaving the name for the shell to resolve would reintroduce the
% PATH lookup that assertProgramAvailable exists to avoid.
if ~isempty(exePath) && ~strcmp(exePath,strtok(cmd))
    rest = cmd(numel(strtok(cmd))+1:end);
    if any(exePath == ' '), cmd = ['"' exePath '"' rest]; else, cmd = [exePath rest]; end
end
if keepStderr
    % Merge explicitly rather than relying on the host's default. MATLAB's
    % system() folds stderr into the captured string; Octave's does not. An
    % explicit 2>&1 behaves identically on both.
    redirect = ' 2>&1';
elseif ispc
    redirect = ' 2>nul';
else
    redirect = ' 2>/dev/null';
end
if isunix && ~ismac
    % See runViaExe: MATLAB's LD_LIBRARY_PATH can shadow the programs' own
    % runtime libraries. Cleared for this call only.
    oldLd = getenv('LD_LIBRARY_PATH');
    setenv('LD_LIBRARY_PATH','');
    restoreLd = onCleanup(@() setenv('LD_LIBRARY_PATH',oldLd)); %#ok<NASGU>
end
% (HCTSA2 2.2.1: stdin from the null device, as in runViaExe, so a missing input
% file ends the program at once instead of blocking on the keyboard.)
if ispc, nullIn = ' <nul'; else, nullIn = ' </dev/null'; end
[status,out] = system([cmd redirect nullIn]);
out = strrep(out,char(13),'');
end

%-------------------------------------------------------------------------------
function s = readTextFile(f)
s = '';
if ~isfile(f), return; end
fid = fopen(f,'r');
if fid < 0, return; end
s = fread(fid,Inf,'*char')';
fclose(fid);
end

%-------------------------------------------------------------------------------
function checkPathLengths(prog,args)
% TISEAN's Fortran programs hold filenames in fixed-length buffers and gfortran
% truncates a longer name SILENTLY -- the program writes somewhere the caller
% never asked for, and the caller then looks for a file that does not exist.
% That is a nasty failure to diagnose: it surfaces as a missing output far from
% its cause.
%
% The buffers were character*72 as TISEAN shipped them, which a Windows tempname
% reaches on its own. They are now character*512 (see the note in c1.f). This
% check is the backstop for that: it should never fire in practice, and if the
% declarations are ever changed back it says so instead of losing a file.
FORTRAN_NAME_LIMIT = 512;
if ~any(strcmp(prog,{'c1','c2d','c2g','c2t'})), return; end
for i = 1:numel(args)
    a = args{i};
    if numel(a) > FORTRAN_NAME_LIMIT && (any(a == '/') || any(a == '\'))
        error('hctsa:tisean:pathTooLong', ...
            ['TISEAN''s Fortran programs hold filenames in character*%u ' ...
             'buffers and truncate anything longer without warning. %s was ' ...
             'given a path of %u characters:\n  %s'], ...
            FORTRAN_NAME_LIMIT,prog,numel(a),a);
    end
end
end

%-------------------------------------------------------------------------------
function t = localTokenize(cmd)
% Split a command line into argv. Double-quoted runs are kept together.
t = regexp(cmd,'"[^"]*"|\S+','match');
t = strrep(t,'"','');
% Nothing hctsa passes contains a space except, possibly, the temporary
% directory in a data path. A path that is not the last argument (-o <file>
% <data>) must be quoted by the caller, which the split above keeps whole.
% (HCTSA2 2.2.1: the call sites that write -o files quote them; this used to
% rely on the data path being last, which left the -o path split in two.)
% Call sites that pass an unquoted data path LAST still work: if the naive split
% left a last token that is not a file, glue the trailing tokens back together
% until it is. On the usual space-free path the first test succeeds and this
% costs one exist() call.
if numel(t) > 2 && ~isfile(t{end})
    for k = numel(t)-1:-1:2
        cand = strjoin(t(k:end),' ');
        if isfile(cand)
            t = [t(1:k-1), {cand}];
            break
        end
    end
end
end

%-------------------------------------------------------------------------------
function cleanFiles(files)
% delete() WARNS rather than errors when a file is already gone, so try/catch
% does not silence it -- the warning has to be turned off around the call. The
% On Windows a name can still be visible to isfile for a moment after its
% handle closes, so this fires routinely and means nothing. Only the two specific identifiers are
% suppressed, and the previous state is restored.
ws = warning('off','MATLAB:DELETE:FileNotFound');
ws(end+1) = warning('off','MATLAB:DELETE:Permission');
restore = onCleanup(@() warning(ws)); %#ok<NASGU>
for i = 1:numel(files)
    if isempty(files{i}), continue; end
    try %#ok<TRYNC>
        if isfile(files{i}), delete(files{i}); end
    catch
    end
end
end
