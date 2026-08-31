function [status,out] = BF_RunTisean(cmd,keepStderr)
% BF_RunTisean   Run a TISEAN program and return its data.
%
%   [status,out] = BF_RunTisean('nstat_z -# 4 -d1 -m3 /tmp/hc123456')
%   [status,out] = BF_RunTisean('d2 ...',true)   % keep stderr
%
% Dispatches to the TISEAN MEX files when they are present and falls back to
% the standalone executables otherwise. Callers see the same two outputs either
% way, so no operation needs to know which is in use.
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
%---CHOOSING THE BACKEND
% Set the environment variable HCTSA2_TISEAN_BACKEND to 'mex' or 'exe' to force
% one; anything else (or unset) means "MEX if built, executables if not". It is
% read once per MATLAB process, so on parallel workers set it with
%
%   pctRunOnAll setenv('HCTSA2_TISEAN_BACKEND','exe')
%
% This exists so the two paths can be diffed against each other on the target
% machine. HCTSA2_TestTiseanMex does exactly that.
%
%---WHY THE MEX PATH REDIRECTS DATA TO A FILE
% A MEX file cannot hand its stdout back the way system() does. The five
% programs that write their data to stdout all accept -o, and their -o output
% is byte-for-byte identical to their stdout output (verified for all five), so
% this function supplies a temporary -o and reads it back. The command strings
% the operations build are untouched.
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
% Neither can happen on the MEX path, where the two streams go to different
% places by construction. On the executable path stderr is redirected away.
% Callers parse by content rather than by position relative to a stderr marker,
% which is more robust on every platform.

% ------------------------------------------------------------------------------
% Copyright (C) 2020, Ben D. Fulcher <ben.d.fulcher@gmail.com>
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version. See <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

persistent haveMex
if isempty(haveMex), haveMex = struct(); end
if nargin < 2 || isempty(keepStderr), keepStderr = false; end

% Read the override on every call rather than caching it. Caching meant the
% choice could only be changed by clearing this function, and `clear <name>`
% from inside another function is not something to depend on -- if it silently
% did nothing, a test that switches backends would compare one backend with
% itself and report everything as matching. getenv is nothing beside writing a
% temp file and running the program.
backendMode = lower(strtrim(getenv('HCTSA2_TISEAN_BACKEND')));
if ~any(strcmp(backendMode,{'mex','exe'})), backendMode = 'auto'; end

args = localTokenize(cmd);
prog = args{1};
args = args(2:end);
checkPathLengths(prog,args);

useMex = false;
if ~strcmp(backendMode,'exe')
    field = ['p_' prog];        % TISEAN program names are already identifiers
    if ~isfield(haveMex,field)
        haveMex.(field) = (exist(['tisean_' prog],'file') == 3);
    end
    useMex = haveMex.(field);
end

if useMex
    [status,out] = runViaMex(prog,args,keepStderr);
else
    [status,out] = runViaExe(cmd,keepStderr);
end
end

%-------------------------------------------------------------------------------
function [status,out] = runViaMex(prog,args,keepStderr)

% The five programs whose data goes to stdout when -o is absent.
stdoutProgs = {'nstat_z','false_nearest','c2d','c2g','c2t'};

dataFile = '';
if any(strcmp(prog,stdoutProgs)) && ~any(strncmp(args,'-o',2))
    % Derived from the input path, which BF_WriteTempFile already keeps short:
    % the Fortran programs declare their filenames character*72 and silently
    % truncate anything longer.
    dataFile = [args{end} '.tso'];
    args = [args(1:end-1), {'-o',dataFile}, args(end)];
end

capture = [tempname '.tiseanerr'];
cleanup = onCleanup(@() cleanFiles({capture,dataFile}));

[status,msg] = feval(['tisean_' prog],args,capture);

data = '';
if ~isempty(dataFile) && isfile(dataFile)
    data = fileread(dataFile);
end

% Windows executables emit CRLF; strip the CR so numeric parsing works. The MEX
% files inherit the same C runtime, so they do too.
data = strrep(data,char(13),'');
msg  = strrep(msg,char(13),'');

if keepStderr
    out = [data msg];
else
    out = data;
end
end

%-------------------------------------------------------------------------------
function [status,out] = runViaExe(cmd,keepStderr)
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
% MEX file removes the capture itself when it reads it back, and on Windows the
% name can still be visible to isfile for a moment after the handle closes, so
% this fires routinely and means nothing. Only the two specific identifiers are
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
