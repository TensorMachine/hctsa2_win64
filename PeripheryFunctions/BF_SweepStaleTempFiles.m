function nRemoved = BF_SweepStaleTempFiles(maxAgeHours)
% BF_SweepStaleTempFiles   Remove abandoned TISEAN temporary file sets.
%
%   nRemoved = BF_SweepStaleTempFiles        % sets untouched for 24 hours
%   nRemoved = BF_SweepStaleTempFiles(48)
%
% Every TISEAN operation removes its own temporary files on every way out
% (HCTSA2_AuditTiseanCleanup proves that). What in-process cleanup cannot
% catch is a MATLAB that never gets to run it: Ctrl-C at the wrong moment, a
% killed worker, a crash. Their files stay in the temporary directory. TS_Compute
% calls this once before it starts, so they do not accumulate.
%
% A file set is hc###### plus everything derived from it (hc######.c2,
% hc######_c1, hc######.tout, ...), the names BF_WriteTempFile and the TISEAN
% programs produce. A set is removed only when its NEWEST file is older than
% MAXAGEHOURS. A run in progress -- in this MATLAB, on a parallel worker, or in
% another MATLAB entirely -- writes its files seconds before using them, so it
% is never touched. Anything not named hc###### is never touched either
% (BF_CleanTempFiles does the deleting and matches as narrowly).
%
% Returns the number of files removed.

% ------------------------------------------------------------------------------
% Copyright (C) 2026, hctsa2 contributors. This function is free software: you
% can redistribute it and/or modify it under the terms of the GNU General Public
% License as published by the Free Software Foundation, either version 3 of the
% License, or (at your option) any later version.
% See <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

if nargin < 1 || isempty(maxAgeHours)
    maxAgeHours = 24;
end
nRemoved = 0;

listing = dir(fullfile(tempdir,'hc*'));
listing = listing(~[listing.isdir]);
if isempty(listing), return; end

tok = regexp({listing.name},'^(hc\d{6})($|[._])','tokens','once');
isOurs = ~cellfun(@isempty,tok);
if ~any(isOurs), return; end
stems = cellfun(@(t) t{1},tok(isOurs),'UniformOutput',false);
dates = [listing(isOurs).datenum];

[uStems,~,g] = unique(stems);
newest = accumarray(g(:),dates(:),[],@max);
cutoff = now - maxAgeHours/24;
for k = find(newest(:)' < cutoff)
    BF_CleanTempFiles(fullfile(tempdir,uStems{k}));
    nRemoved = nRemoved + sum(g == k);
end

end
