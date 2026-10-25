function filePath = BF_WriteTempFile(dataVector,whatPrecision,delimiter)
% BF_WriteTempFile   Write a temporary file in the system temporary directory
%
%---INPUTS:
% dataVector, a vector of data to write to file, or a matrix (one row per line).
% whatPrecision, the number of significant digits to write (default 7).
% delimiter, the column delimiter for a matrix (default ',', as dlmwrite's
%       default, so every existing call writes exactly what it always has). TISEAN
%       programs that read several columns need whitespace: pass ' '. ripser, the
%       other reader of a matrix (NL_PersistentHomology), accepts either.

% ------------------------------------------------------------------------------
% Copyright (C) 2020, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite the following two papers:
%
% (1) B.D. Fulcher and N.S. Jones, "hctsa: A Computational Framework for Automated
% Time-Series Phenotyping Using Massive Feature Extraction, Cell Systems 5: 527 (2017).
% DOI: 10.1016/j.cels.2017.10.001
%
% (2) B.D. Fulcher, M.A. Little, N.S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013).
% DOI: 10.1098/rsif.2013.0048
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program. If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

if nargin < 2
    % Precision (number of significant digits) to write: 17, enough to represent
    % every double exactly, so the TISEAN program reads the same numbers that MATLAB
    % holds. (With 7 digits, quantized or near-tied data gave different neighbor
    % sets from the exact values, e.g. in lyap_r and false_nearest.)
    whatPrecision = 17;
end
if nargin < 3
    delimiter = ',';
end

% Short temp filename, forward slashes.
%
% (HCTSA2 2.2.1: only the LEAF name is short and space-free. The folder is
% tempdir as it is, and it can contain a space -- a user-set TEMP such as
% D:\My Temp, or a user name with a space where 8.3 short names are off. Callers
% that put this path on a TISEAN command line must quote it: "%s".)
%
% The TISEAN executables receive this path on a command line. MATLAB's tempname
% produces something like
%   C:\Users\<you>\AppData\Local\Temp\tp9b5f2030_a76e_4847_bb58_fb718fba2436
% -- about 70 characters, with backslashes. The Windows TISEAN 3.x binaries fail
% on it: nstat_z returns nothing at all, d2 reports "Cannot open input file",
% and c1 never creates its output. hctsa's own error text already suspected
% this ("Could be due to an overly long filename").
%
% A short name in the temp directory avoids it, and forward slashes work on
% Windows while surviving the shell without escaping.
%
%---WHY THE NAME COMES FROM A PRIVATE GENERATOR AND IS CLAIMED ATOMICALLY
% Parallel workers share one temporary directory. The name used to be drawn
% with randi() from the GLOBAL random stream, and operations reset that stream
% with BF_ResetSeed, so two workers running the same operations reached this
% point in the same RNG state and drew the SAME name. TS_Compute's
% series-parallel mode does exactly that. The check that followed (exist, then
% write) is not atomic across processes, so both could take the name, and one
% worker's cleanup then deleted the other's files mid-run.
%
% So: the name comes from a generator private to this MATLAB process, seeded
% from its process id and the time (drawName), which also stops this function
% from moving the caller's random stream; and the name is claimed with an
% atomic create-if-absent (claimName). A claim can only fail if the name is
% taken, and then the next name is tried.
filePath = '';
for attempt = 1:20
    candidate = fullfile(tempdir,sprintf('hc%06u',drawName()));
    if ispc
        candidate = strrep(candidate,'\','/');
    end
    if claimName(candidate)
        filePath = candidate;
        break
    end
end
if isempty(filePath)
    error('hctsa:tempfile:noName', ...
        'Could not claim an unused temporary file name in %s after 20 attempts.', ...
        tempdir);
end

% Write the file. The name is already claimed (an empty file exists), so if the
% write fails, remove it here: the caller's cleanup only starts after this returns.
try
    dlmwrite(filePath,dataVector,'delimiter',delimiter,'precision',whatPrecision);
catch err
    if exist(filePath,'file'), delete(filePath); end
    rethrow(err);
end

end

%-------------------------------------------------------------------------------
function n = drawName()
% Six digits from a random stream private to this process: the Park-Miller
% minimal standard generator (multiplier 48271, modulus 2^31-1). Its products
% stay below 2^53, so it is exact in double arithmetic, and it behaves the same
% in MATLAB and Octave. Seeded once per process from the process id and the
% time of day in microseconds, so parallel workers draw different sequences.
persistent state
if isempty(state)
    usOfDay = floor(mod(now,1)*864e8);
    state = mod(processId()*104729 + usOfDay,2147483646) + 1;   % 1 .. 2^31-2
end
state = mod(state*48271,2147483647);
n = mod(state,1000000);
end

%-------------------------------------------------------------------------------
function pid = processId()
pid = 0;
try
    if exist('OCTAVE_VERSION','builtin')
        pid = getpid();
    else
        pid = feature('getpid');
    end
catch
end
end

%-------------------------------------------------------------------------------
function ok = claimName(f)
% Create F if and only if it does not exist yet, as one atomic operation.
% java.io.File.createNewFile guarantees that across processes. Without a JVM,
% fall back to the old check-then-write.
if usejava('jvm')
    try
        jf = javaObject('java.io.File',f);
        ok = logical(jf.createNewFile());
        return
    catch
    end
end
ok = ~exist(f,'file');
end
