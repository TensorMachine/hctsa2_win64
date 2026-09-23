function BF_CleanTempFiles(stem)
% BF_CleanTempFiles   Remove one TISEAN temporary file set, silently.
%
%   BF_CleanTempFiles(filePath)
%
% Removes the temporary time-series file BF_WriteTempFile created and every
% file TISEAN derived from it -- .c2, .d2, .h2, .stat, .box, .ros, .poin, .c1,
% .fnn, .tso and any other suffix a program invents. Safe to call more than
% once, and safe to call when nothing is there.
%
% STEM may be the temporary file itself or any path derived from it, so
% BF_CleanTempFiles('/tmp/hc123456.c2') cleans the whole set.
%
%---WHY THIS MATCHES SO NARROWLY
% It only ever touches a name of the form hc###### -- exactly six digits, the
% shape BF_WriteTempFile produces -- optionally followed by "." or "_" and a
% suffix. Anything else is left strictly alone.
%
% That precision is the entire point. The obvious implementation, deleting
% tempdir/hc*, is unsafe: a parallel pool shares one temporary directory, so a
% blanket wildcard would delete the INPUT FILES OF OTHER WORKERS mid-run. The
% failures that follow -- a program reporting a missing file it was handed
% moments ago, on a random subset of time series, never reproducibly -- would
% be very hard to trace back to the cleanup that caused them.
%
% Matching exactly six digits also means one call's stem can never be a prefix
% of another's: hc000123 and hc123456 are the same length, so hc12345 cannot
% swallow hc123456. The "." or "_" requirement covers the rest, since c1 names
% its default output file_c1 with an underscore rather than a dot.

% ------------------------------------------------------------------------------
% Copyright (C) 2020, Ben D. Fulcher <ben.d.fulcher@gmail.com>
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version. See <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

if nargin < 1 || isempty(stem) || ~ischar(stem), return; end

[folder,name,ext] = fileparts(stem);
tok = regexp([name ext],'^(hc\d{6})','tokens','once');
if isempty(tok), return; end        % not one of ours; never touch it
base = tok{1};
if isempty(folder), folder = tempdir; end

% delete() WARNS rather than errors on a file that is already gone, so try/catch
% cannot silence it -- the warning has to be turned off around the call. This
% fires routinely and harmlessly: the operations delete some of these files
% themselves, and on Windows a name can stay visible to dir() for a moment
% after its handle closes.
ws = warning('off','MATLAB:DELETE:FileNotFound');
ws(end+1) = warning('off','MATLAB:DELETE:Permission');
restore = onCleanup(@() warning(ws)); %#ok<NASGU>

listing = dir(fullfile(folder,[base '*']));
for i = 1:numel(listing)
    if listing(i).isdir, continue; end
    % "hc123456", "hc123456.c2", "hc123456_c1" -- and nothing else.
    if isempty(regexp(listing(i).name,['^' base '($|[._])'],'once'))
        continue
    end
    f = fullfile(folder,listing(i).name);
    try %#ok<TRYNC>
        if isfile(f), delete(f); end
    catch
    end
end
