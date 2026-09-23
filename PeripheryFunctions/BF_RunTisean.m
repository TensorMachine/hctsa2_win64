function [status,out] = BF_RunTisean(cmd,keepStderr)
% BF_RunTisean   Run a TISEAN program and return its data.
%
%   [status,out] = BF_RunTisean('nstat_z -# 4 -d1 -m3 /tmp/hc123456')
%   [status,out] = BF_RunTisean('d2 ...',true)   % keep stderr
%
% Runs one of the ten TISEAN executables and returns its exit status and its
% text. Every TISEAN operation goes through here rather than calling system()
% directly, so the timeout, the not-installed check and the temporary-file
% cleanup apply uniformly.
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
% what TISEAN just wrote: NL_TISEAN_d2 runs d2, then feeds the .c2 it produced
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
    assertProgramAvailable(prog);
    [status,out] = runViaExe(prog,args,cmd,keepStderr,timeoutSeconds());
catch err
    if ~isempty(args)
        BF_CleanTempFiles(args{end});
    end
    rethrow(err);
end
end

%-------------------------------------------------------------------------------
function assertProgramAvailable(prog)
% Raise a single, actionable error when neither backend can supply PROG.
%
% The probe runs once per program per MATLAB session, so the installed case --
% the normal one -- costs a single 'where'/'command -v' and nothing after that.
persistent haveExe
if isempty(haveExe), haveExe = struct(); end
field = ['p_' prog];
if ~isfield(haveExe,field)
    if ispc
        probe = sprintf('where %s >nul 2>&1',prog);
    else
        probe = sprintf('command -v %s >/dev/null 2>&1',prog);
    end
    haveExe.(field) = (system(probe) == 0);
end
if haveExe.(field), return; end

error('hctsa:tisean:notInstalled', ...
    ['TISEAN program ''%s'' is not on the PATH, so no TISEAN feature can be ' ...
     'computed on this machine.\n  The executables are BUILD ARTIFACTS: they ' ...
     'are machine-specific, they are listed in .gitignore, and they do not ' ...
     'travel with the repository.\n  Build them on this machine with ' ...
     'HCTSA2_BuildTisean, then run startup.\n  HCTSA2_CheckInstall reports ' ...
     'what is present.'],prog);
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
function [status,out] = runViaExe(prog,args,cmd,keepStderr,timeoutSec)
% Run the executable, abandoning it if it outstays TIMEOUTSEC.
%
%---WHY A TIMEOUT IS NEEDED AT ALL
% Neither d2 nor false_nearest can loop forever -- every loop in both is
% bounded, and the length guards that would otherwise let an unsigned
% subtraction wrap are intact. What they can do is take a quadratic number of
% steps with no way to stop early:
%
%   d2            -N0 sets the pair limit to ULONG_MAX (d2.c:120), which makes
%                 the only early exit, "if (imin > howoften1) exit(0)" at
%                 d2.c:561, unreachable. d2 then always runs its full
%                 nmax-iteration loop with an O(n) box scan inside each one.
%                 Measured at -M1,10: 4k 2.3s, 8k 9.5s, 16k 38.4s -- four times
%                 the work for twice the data. Extrapolating, 500k points is
%                 about ten hours. NL_dimensions passes -N0 outright;
%                 NL_TakensEstimator and NL_GPCorrSum reach it through
%                 Nref = -1. NL_TISEAN_d2 passes no -N, gets the default 1000,
%                 and is unaffected -- which is why only some features hang.
%
%   false_nearest has no pair limit at all. A point whose nearest neighbour
%                 lies farther than varianz/rt can never be marked done
%                 (false_nearest.c:168) and exact-duplicate vectors are
%                 rejected outright (:158), so a constant or near-constant
%                 stretch leaves points that are rescanned at every epsilon
%                 step while all sitting in one search box. Measured at
%                 N = 16000: 4.1s on clean data, 41.8s when half the series is
%                 flat.
%
% Because both terminate eventually, killing the process is sound: nothing is
% deadlocked, the run is merely too expensive to wait for.
%
%---WHY THIS RAISES AN ERROR RATHER THAN RETURNING WHAT IT GOT
% d2 rewrites its .c2/.d2/.h2 files every 120 seconds (d2.c:31), so a killed
% run usually leaves a COMPLETE-LOOKING but partial set of output files on
% disk. Returning them would hand the caller plausible numbers computed from a
% fraction of the data -- silently wrong features, which is worse than none.
% An error becomes a NaN in the feature matrix, which is honest.
%
%---WHY JAVA AND NOT system()
% system() blocks until the child exits; there is no way to take it back. Java
% is present in any MATLAB that runs hctsa (startup javaaddpath's
% infodynamics.jar), and it gives both a timed wait and a kill. It also takes
% argv as a list, so nothing passes through a shell and destroyForcibly kills
% the TISEAN process itself rather than a shell that owns it. When there is no
% JVM, or Java cannot start the program, this falls back to the old blocking
% system() call rather than failing.
persistent warned
if timeoutSec <= 0 || ~usejava('jvm')
    [status,out] = runViaSystem(cmd,keepStderr);
    return
end

outFile = [tempname '.tout'];
errFile = [tempname '.terr'];
cleanup = onCleanup(@() cleanFiles({outFile,errFile})); %#ok<NASGU>

try
    jargs = javaArray('java.lang.String',numel(args)+1);
    jargs(1) = javaObject('java.lang.String',prog);
    for i = 1:numel(args)
        jargs(i+1) = javaObject('java.lang.String',args{i});
    end
    pb = javaObject('java.lang.ProcessBuilder',jargs);
    pb.redirectOutput(javaObject('java.io.File',outFile));
    if keepStderr
        % Merged into stdout, matching what system() hands back.
        pb.redirectErrorStream(true);
    else
        pb.redirectError(javaObject('java.io.File',errFile));
    end
    proc = pb.start();
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
    [status,out] = runViaSystem(cmd,keepStderr);
    return
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
function [status,out] = runViaSystem(cmd,keepStderr)
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
[status,out] = system([cmd redirect]);
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
% directory in the data path -- which is the last argument in every call site.
% If the naive split left a last token that is not a file, glue the trailing
% tokens back together until it is. On the usual space-free path the first test
% succeeds and this costs one exist() call.
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
