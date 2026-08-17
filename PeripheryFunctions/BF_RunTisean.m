function [status,out] = BF_RunTisean(cmd,keepStderr)
% BF_RunTisean   Run a TISEAN executable and capture ONLY its stdout.
%
%   [status,out] = BF_RunTisean('nstat_z -# 4 -d1 -m3 /tmp/hc123456')
%   [status,out] = BF_RunTisean('d2 ...',true)   % keep stderr
%
% KEEPSTDERR (default false) controls whether stderr is merged in.
%
%   false -- for programs whose DATA is on stdout (nstat_z, c2g, c2t). Merging
%            stderr corrupts the parse; see below.
%   true  -- for programs that write their results to FILES and use stdout for
%            nothing at all (d2, c1). Their only text is the stderr banner, so
%            discarding it returns an empty string even on success, and callers
%            that test isempty(res) as a did-it-run check then fail.
%
% MATLAB's system() merges stdout and stderr into one string. The TISEAN
% programs write their data to stdout and their banner and progress messages to
% stderr, and the two streams are buffered differently -- stderr unbuffered,
% stdout block-buffered when redirected. On Windows this produces two failures
% that are invisible on Linux and macOS:
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
% Redirecting stderr away fixes both: what comes back is exactly the data, in
% order, with no interleaving. Callers must then parse by content rather than
% by position relative to a stderr marker, which is more robust anyway.

if nargin < 2 || isempty(keepStderr), keepStderr = false; end

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

% Windows executables emit CRLF; strip the CR so numeric parsing works.
out = strrep(out,char(13),'');
end
