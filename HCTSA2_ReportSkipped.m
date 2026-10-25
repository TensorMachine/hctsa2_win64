function report = HCTSA2_ReportSkipped(whatData,verbose,installReport)
% HCTSA2_ReportSkipped   Explain which features did not compute, and why.
%
%   HCTSA2_ReportSkipped                    % uses HCTSA.mat
%   HCTSA2_ReportSkipped('HCTSA_myrun.mat')
%   report = HCTSA2_ReportSkipped('HCTSA.mat',false)
%   report = HCTSA2_ReportSkipped('HCTSA.mat',true,HCTSA2_CheckInstall(false))
%
% With the install report from HCTSA2_CheckInstall, a failure is blamed on a
% compiled component only when that component is actually missing; without it,
% the component is named as the likely cause from the operation name alone.
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
%   6  output was empty
%   7  link error (the registered output field does not exist)
%
% Codes 2-5 are usually legitimate results for a particular time series, not
% installation problems. Code 1 is what this function is really for.

if nargin < 1 || isempty(whatData), whatData = 'HCTSA.mat'; end
if nargin < 2, verbose = true; end
if nargin < 3, installReport = []; end

if ~exist(whatData,'file')
    error('HCTSA2_ReportSkipped:noFile','Could not find %s',whatData);
end
D = load(whatData,'TS_Quality','Operations','MasterOperations');
if ~isfield(D,'TS_Quality')
    error('HCTSA2_ReportSkipped:noQuality','%s has no TS_Quality.',whatData);
end
Q = D.TS_Quality;
[nTS,nOps] = size(Q);

% Map each compiled component to the operation-name prefixes it gates (counts
% as in HCTSA2_CheckInstall, the 3.0.0 feature set of upstream v3.0.0), so a
% wall of "undefined function" errors becomes one actionable line. The prefixes
% are matched against the feature name and its master operation's label (the
% catch22 feature sets name their features 'mode_5' etc. but their masters
% 'catch22_...'); operations reaching TISEAN's false_nearest through an 'fnn'
% embedding argument (50 features) are matched by their master code, and so are the
% users of the KSG MEX (2.4.0): a 'kraskov' estimator or a 'mi' time delay.
comp = {
 'TISEAN/ripser (executables)', {'NL_d2','NL_c1','NL_FNN','SY_nstat_z','NL_BoxCountEntropyRate', ...
                                'NL_Dimensions','NL_LargestLyap', ...
                                'NL_GPCorrSum','NL_LyapSpec','PP_SchreiberDenoise', ...
                                'NL_PersistentHomology'},                        538
 'Physionet sampen (MEX)',  {'EN_SampEn','EN_MSE'},                             91
 'Michael_Small (MEX)',     {'EN_LZComplexity','EX_ExtremeEventOrder','EN_Shannon'}, 47
 'Max_Little (MEX)',        {'CP_ML_StepDetect','CP_l1pwc_SweepLambda','EN_rpde'}, 36
 'KSG mutual information (MEX)', {'NL_PersistentHomology'},                440
 'catch22 (MEX)',           {'catch22'},                                         0
};

report = struct();
report.file = whatData;
report.nTimeSeries = nTS;
report.nOperations = nOps;
report.total = numel(Q);

codes = [0 1 2 3 4 5 6 7];
labels = {'computed normally','fatal error','NaN output','+Inf output', ...
          '-Inf output','complex output','empty output','link error'};
counts = zeros(1,numel(codes));
for i = 1:numel(codes)
    counts(i) = sum(Q(:) == codes(i));
end
report.counts = counts;
report.nanQuality = sum(isnan(Q(:)));

if verbose
    fprintf('\n========== hctsa2: computation outcome ==========\n');
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

% Master label and code of each operation (when the file has them)
mLabel = repmat({''},1,nOps); mCode = repmat({''},1,nOps);
try
    [tf,loc] = ismember(D.Operations.MasterID,D.MasterOperations.ID);
    mLabel(tf) = D.MasterOperations.Label(loc(tf));
    mCode(tf) = D.MasterOperations.Code(loc(tf));
catch
end

% Which components are missing, when the install report says so
missing = true(1,size(comp,1));
if ~isempty(installReport)
    missing = cellfun(@(c) localMissing(c,installReport),comp(:,1))';
end

% Attribute failures to a (missing) compiled component where the name matches
blame = zeros(1,numel(badOps));
for k = 1:numel(badOps)
    j = badOps(k);
    for c = 1:size(comp,1)
        hit = any(cellfun(@(p) strncmp(opNames{j},p,numel(p)) || ...
                               strncmp(mLabel{j},p,numel(p)),comp{c,2}));
        if c == 1, hit = hit || ~isempty(strfind(mCode{j},'''fnn''')); end
        if strncmp(comp{c,1},'KSG',3)   % (NL_PersistentHomology's delay defaults to 'mi')
            hit = hit || ~isempty(regexp(mCode{j},'kraskov|''mi''','once'));
        end
        if hit && missing(c)
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
        if isempty(installReport)
            whyText = 'not attributable';
        else
            whyText = 'data-dependent';  % components present: errors on these series
        end
        fprintf('  %-26s %5u operations, %7u feature values\n', ...
            whyText,sum(sel),nf);
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

%===============================================================================
function tf = localMissing(compName,R)
% Is this component unavailable according to HCTSA2_CheckInstall's report?
tf = false;
if strncmp(compName,'TISEAN',6)
    tf = isfield(R,'tisean') && ~R.tisean;
    return
end
key = strtok(compName,' ');   % 'Physionet', 'Michael_Small', 'Max_Little', 'catch22'
if isfield(R,'mex')
    for i = 1:numel(R.mex)
        if strncmp(R.mex(i).name,key,numel(key)) && ~R.mex(i).available && ~R.mex(i).fallback
            tf = true;
        end
    end
end
end
