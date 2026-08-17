function report = HCTSA2_ReportSkipped(whatData,verbose)
% HCTSA2_ReportSkipped   Explain which features did not compute, and why.
%
%   HCTSA2_ReportSkipped                    % uses HCTSA.mat
%   HCTSA2_ReportSkipped('HCTSA_myrun.mat')
%   report = HCTSA2_ReportSkipped('HCTSA.mat',false)
%
% hctsa never aborts a run because one operation failed: each master operation
% is wrapped in try/catch and the outcome is recorded per feature in
% TS_Quality. That is good for robustness but it means failures scroll past
% during a long run and are easy to miss.
%
% This groups the outcome by cause after the fact, and — where the failure is a
% missing compiled component — names the component and points at the fix.
%
%---TS_Quality codes (hctsa's own convention)
%   0  computed normally
%   1  fatal error in the operation
%   2  output was NaN
%   3  output was +Inf
%   4  output was -Inf
%   5  output was complex
%
% Codes 2-5 are usually legitimate results for a particular time series, not
% installation problems. Code 1 is what this function is really for.

if nargin < 1 || isempty(whatData), whatData = 'HCTSA.mat'; end
if nargin < 2, verbose = true; end

if ~exist(whatData,'file')
    error('HCTSA2_ReportSkipped:noFile','Could not find %s',whatData);
end
D = load(whatData,'TS_Quality','Operations','TimeSeries');
if ~isfield(D,'TS_Quality')
    error('HCTSA2_ReportSkipped:noQuality','%s has no TS_Quality.',whatData);
end
Q = D.TS_Quality;
[nTS,nOps] = size(Q);

% Map each compiled component to the operation-name prefixes it gates, so a
% wall of "undefined function" errors becomes one actionable line.
comp = {
 'catch22 (MEX)',           {'catch22'},                                      240
 'TISEAN (executables)',    {'NL_TISEAN','SY_Tisean','MF_tisean'},            276
 'Physionet sampen (MEX)',  {'EN_SampEn','EN_Randomize','EN_mse'},             49
 'Michael_Small (MEX)',     {'MS_','NL_MS'},                                   33
 'Max_Little (MEX)',        {'SC_fastdfa','PH_Walker','ML_'},                  10
};

report = struct();
report.file = whatData;
report.nTimeSeries = nTS;
report.nOperations = nOps;
report.total = numel(Q);

codes = [0 1 2 3 4 5];
labels = {'computed normally','fatal error','NaN output','+Inf output', ...
          '-Inf output','complex output'};
counts = zeros(1,numel(codes));
for i = 1:numel(codes)
    counts(i) = sum(Q(:) == codes(i));
end
report.counts = counts;
report.nanQuality = sum(isnan(Q(:)));

if verbose
    fprintf('\n========== hctsa 2.0: computation outcome ==========\n');
    fprintf('%s  --  %u time series x %u operations\n\n',whatData,nTS,nOps);
    for i = 1:numel(codes)
        if counts(i) > 0
            fprintf('  %-20s %9u  (%.1f%%)\n',labels{i},counts(i), ...
                100*counts(i)/numel(Q));
        end
    end
    if report.nanQuality > 0
        fprintf('  %-20s %9u  (never attempted)\n','not computed',report.nanQuality);
    end
end

%-------------------------------------------------------------------------------
%% Which operations failed outright, and can we name the cause?
%-------------------------------------------------------------------------------
failedPerOp = sum(Q == 1,1);
badOps = find(failedPerOp > 0);
report.failedOps = {};
report.blamed = {};

if isempty(badOps)
    if verbose
        fprintf('\nNo fatal errors. Every operation ran.\n\n');
    end
    return
end

% Operations is a table in MATLAB; tolerate a struct or cell array too so the
% report still works on partially-loaded or older files.
opNames = repmat({''},1,nOps);
if isfield(D,'Operations')
    O = D.Operations;
    nm = {};
    % Order matters: check the plain types first. `istable` does not exist in
    % older MATLAB or in Octave, so testing it first aborts the whole block.
    if iscell(O)
        nm = O;
    elseif isstruct(O) && isfield(O,'Name')
        nm = {O.Name};
    else
        try
            nm = O.Name;                 % table
        catch
            nm = {};
        end
    end
    if ~iscell(nm), nm = cellstr(nm); end
    if numel(nm) == nOps, opNames = nm(:)'; end
end

% Attribute failures to a compiled component where the name matches
blame = zeros(1,numel(badOps));
for k = 1:numel(badOps)
    nm = opNames{badOps(k)};
    for c = 1:size(comp,1)
        if any(cellfun(@(p) strncmp(nm,p,numel(p)),comp{c,2}))
            blame(k) = c; break
        end
    end
end

if verbose
    fprintf('\n-- fatal errors, grouped by likely cause --\n');
    for c = 1:size(comp,1)
        sel = (blame == c);
        if any(sel)
            nf = sum(failedPerOp(badOps(sel)));
            fprintf('  %-26s %5u operations, %7u feature values\n', ...
                comp{c,1},sum(sel),nf);
        end
    end
    sel = (blame == 0);
    if any(sel)
        nf = sum(failedPerOp(badOps(sel)));
        fprintf('  %-26s %5u operations, %7u feature values\n', ...
            'not attributable',sum(sel),nf);
        show = badOps(sel);
        fprintf('      e.g. ');
        for k = 1:min(5,numel(show))
            fprintf('%s ',opNames{show(k)});
        end
        fprintf('\n');
    end
    if any(blame > 0)
        fprintf('\n  Missing compiled components are the likely cause above.\n');
        fprintf('  Run HCTSA2_CheckInstall to confirm, then HCTSA2_BuildBinaries\n');
        fprintf('  (or see INSTALL.md section 4 for TISEAN).\n');
    end
    fprintf('\n  These failures did not stop the run: every other feature was\n');
    fprintf('  computed and stored normally.\n\n');
end

report.failedOps = opNames(badOps);
report.blamed = comp(unique(blame(blame>0)),1);
end
