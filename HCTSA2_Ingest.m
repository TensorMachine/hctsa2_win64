function [tsData,labels,groups,dropped] = HCTSA2_Ingest(data,opt)
% HCTSA2_Ingest   Read time series out of a .mat, a .csv, or a variable.
%
%   [tsData,labels,groups,dropped] = HCTSA2_Ingest('mydata.mat',opt)
%   [tsData,labels,groups,dropped] = HCTSA2_Ingest(X,opt)     % matrix, cell or table
%
% OPT is a struct that may carry: DataVar, Labels, Groups, IdVar, TimeVar,
% ValueVar, ClassVar, Format ('auto' | 'long' | 'wide'). Missing fields are
% filled with defaults, so a bare struct() is fine.
%
% Handles four layouts:
%   - LONG-FORMAT TABLE  (id, time, value, class) -- one row per observation
%   - wide table         (one column per series)
%   - numeric matrix     (rows = series; long format is detected too)
%   - cell array         (one vector per series)
%
% Missing values: NaN/Inf at either end of a series is padding and is trimmed.
% A series with NaN/Inf inside it, nothing left after trimming, or more than
% 50,000 samples (hctsa's limit) is dropped and listed in DROPPED (a struct of
% two cell columns, Name and Reason; struct2table makes it a table), rather than
% silently spliced or skipped later.
%
% Shared by HCTSA2_Analyse and HCTSA2_SamplingSafeFeatures so both read data
% identically.
%
% (Type tests use isa(x,'categorical') and isa(x,'string') rather than
% iscategorical/isstring: equivalent in MATLAB, and also available in Octave,
% where parts of this code are smoke-tested.)

defaults = struct('DataVar','','Labels',{{}},'Groups',[], ...
                  'IdVar',[],'TimeVar',[],'ValueVar',[],'ClassVar',[],'Format','auto');
fn = fieldnames(defaults);
for i = 1:numel(fn)
    if ~isfield(opt,fn{i}), opt.(fn{i}) = defaults.(fn{i}); end
end
opt.Format = validatestring(char(opt.Format),{'auto','long','wide'});

[tsData,labels,groups] = localIngest(data,opt);
[tsData,labels,groups,dropped] = localMissing(tsData,labels,groups);
if isempty(tsData)
    error('HCTSA2_Analyse:noUsableSeries', ...
        'No usable time series: all %u were dropped (see the list above).',numel(dropped.Name));
end
end

%===============================================================================
function [tsData,labels,groups] = localIngest(data,opt)

if ischar(data) || (isa(data,'string') && isscalar(data))
    dataFile = char(data);
    if exist(dataFile,'file') ~= 2
        error('HCTSA2_Analyse:noFile','Could not find %s',dataFile);
    end
    [~,~,ext] = fileparts(dataFile);
    if strcmpi(ext,'.csv')
        D = struct('data',readtable(dataFile));
    else
        D = load(dataFile);
    end
else
    % in-memory data: groups and labels must be given as values
    D = struct('data',{data});
    dataFile = '(variable)';
end
fn = fieldnames(D);

raw = [];
if ~isempty(opt.DataVar)
    if ~isfield(D,opt.DataVar)
        error('HCTSA2_Analyse:noVar','%s has no variable %s',dataFile,opt.DataVar);
    end
    raw = D.(opt.DataVar);
else
    for i = 1:numel(fn)                      % a table wins over anything else
        if isa(D.(fn{i}),'table') || isa(D.(fn{i}),'timetable')
            raw = D.(fn{i});
            fprintf('    using table variable %s\n',fn{i});
            break
        end
    end
    if isempty(raw)
        best = 0;
        for i = 1:numel(fn)
            v = D.(fn{i});
            if iscell(v) && ~isempty(v) && isnumeric(v{1})
                raw = v; break
            elseif isnumeric(v) && ismatrix(v) && numel(v) > best && min(size(v)) > 1
                raw = v; best = numel(v);
            end
        end
    end
    if isempty(raw)
        error('HCTSA2_Analyse:noData', ...
            ['Could not find time-series data in %s. Variables: %s. ' ...
             'Pass DataVar to name it explicitly.'],dataFile,strjoin(fn,', '));
    end
end

groups = [];
if isa(raw,'timetable'), raw = timetable2table(raw); end

if isa(raw,'table')
    [tsData,labels,groups] = localFromTable(raw,opt);
elseif isnumeric(raw) && ~strcmp(opt.Format,'wide') && ...
        (strcmp(opt.Format,'long') || isLongFormatMatrix(raw))
    % A numeric matrix can be long format too: (id, time, value[, class]) with
    % one row per observation. Without this test such a matrix is transposed by
    % the branch below and read as a handful of very long series -- a 5030-by-4
    % table of 12 series came back as 4 series of length 5030, with no error and
    % entirely plausible-looking output. Silent misreading is worse than a
    % failure, so it is detected explicitly.
    fprintf('    numeric matrix looks LONG format (%u rows, %u columns)\n', ...
        size(raw,1),size(raw,2));
    vn = {'id','time','value','class'};
    T = struct();
    for c = 1:min(4,size(raw,2)), T.(vn{c}) = raw(:,c); end
    [tsData,labels,groups] = localFromStructCols(T,opt);
elseif iscell(raw)
    tsData = cellfun(@(x) double(x(:)),raw(:),'UniformOutput',false);
    labels = {};
else
    if size(raw,1) > size(raw,2)
        raw = raw';
        fprintf('    (transposed: assuming series-by-time)\n');
    end
    tsData = cell(size(raw,1),1);
    for i = 1:size(raw,1)
        tsData{i} = double(raw(i,:))';
    end
    labels = {};
end
n = numel(tsData);

% The variable named for groups is never a candidate for the series names
groupVarName = '';
if ischar(opt.Groups) || (isa(opt.Groups,'string') && isscalar(opt.Groups))
    groupVarName = char(opt.Groups);
end

if ~isempty(opt.Labels), labels = opt.Labels; end
if isempty(labels)
    % Series names: the first cellstr with one UNIQUE entry per series that is
    % not a keyword or class-label variable (a cellstr of group names has
    % repeats, and would make TS_Init fail on duplicate names)
    skip = [{'keywords',groupVarName},{'groups','y','classes','class','group'}];
    for i = 1:numel(fn)
        v = D.(fn{i});
        if iscellstr(v) && numel(v) == n && ~any(strcmpi(fn{i},skip)) && ... %#ok<ISCLSTR>
                numel(unique(v)) == n
            labels = v; break
        end
    end
end
if isempty(labels)
    labels = arrayfun(@(i) sprintf('Series_%u',i),1:n,'UniformOutput',false);
end
labels = cellstr(labels);

if ~isempty(opt.Groups)
    groups = opt.Groups;
    if ~isempty(groupVarName)
        if isfield(D,groupVarName), groups = D.(groupVarName);
        else, error('HCTSA2_Analyse:noGroupVar','No variable %s in %s',groupVarName,dataFile);
        end
    end
end
if isempty(groups)
    for cand = {'groups','y','classes','class','group'}
        if isfield(D,cand{1}) && numel(D.(cand{1})) == n
            groups = D.(cand{1});
            fprintf('    (using %s as class labels)\n',cand{1});
            break
        end
    end
end
if isempty(groups)
    warning('HCTSA2_Analyse:noGroups', ...
        ['No class labels found -- separability analysis will be skipped. ' ...
         'Pass Groups to enable it.']);
    groups = ones(n,1);
end
if numel(groups) ~= n
    error('HCTSA2_Analyse:groupCount','%u group labels for %u time series.',numel(groups),n);
end
groups = groups(:);
end

%===============================================================================
function [tsData,labels,groups,dropped] = localMissing(tsData,labels,groups)
% Trim non-finite padding at both ends; drop series with non-finite values
% inside, nothing left, or more than hctsa's 50,000 samples; list what was dropped.
maxL = 50000; % TS_ReadInputFile's limit
n = numel(tsData);
labels = labels(:);
reason = repmat({''},n,1);
nTrimmed = 0;
for i = 1:n
    x = tsData{i};
    ok = isfinite(x);
    if ~any(ok)
        reason{i} = 'no finite values';
        continue
    end
    first = find(ok,1,'first'); last = find(ok,1,'last');
    if first > 1 || last < numel(x)
        nTrimmed = nTrimmed + 1;
        x = x(first:last);
        tsData{i} = x;
    end
    if ~all(isfinite(x))
        reason{i} = sprintf('%u NaN/Inf values inside the series',sum(~isfinite(x)));
    elseif numel(x) > maxL
        reason{i} = sprintf('%u samples (hctsa accepts up to %u)',numel(x),maxL);
    end
end
if nTrimmed > 0
    fprintf('    %u series had NaN/Inf padding at an end, which was trimmed.\n',nTrimmed);
end
isDropped = ~cellfun(@isempty,reason);
dropped = struct('Name',{labels(isDropped)},'Reason',{reason(isDropped)});
if any(isDropped)
    warning('HCTSA2_Analyse:droppedSeries','%u of %u series dropped at ingest:',sum(isDropped),n);
    for i = find(isDropped)'
        fprintf('      %s: %s\n',labels{i},reason{i});
    end
end
tsData = tsData(~isDropped);
labels = labels(~isDropped);
groups = groups(~isDropped);
end

%===============================================================================
function [tsData,labels,groups] = localFromTable(T,opt)
% Reshape a table into one vector per series.
%
% Columns are taken POSITIONALLY by default -- (id, time, value, class) --
% because column names vary between datasets; override with
% IdVar/TimeVar/ValueVar/ClassVar, each of which accepts a name or a column index.

vn = T.Properties.VariableNames;
idVar    = localPickVar(opt.IdVar,1,vn);
timeVar  = localPickVar(opt.TimeVar,2,vn);
valVar   = localPickVar(opt.ValueVar,3,vn);
classVar = localPickVar(opt.ClassVar,4,vn);

switch opt.Format
case 'long'
    isLong = true;
case 'wide'
    isLong = false;
otherwise
    % Long format needs an id column of few, repeated values (the same rule as
    % for numeric matrices). Any repeat at all was enough before, so a wide table
    % whose first column held integer or quantised data was read as long format.
    userCols = ~isempty(opt.IdVar) || ~isempty(opt.TimeVar) || ~isempty(opt.ValueVar);
    isLong = userCols || (~isempty(idVar) && numel(vn) >= 3 && isLongIdColumn(T.(idVar)));
end

if ~isLong
    fprintf('    table looks WIDE (one column per series)\n');
    num = varfun(@isnumeric,T,'OutputFormat','uniform');
    cols = vn(num);
    tsData = cell(numel(cols),1);
    for i = 1:numel(cols)
        tsData{i} = double(T.(cols{i}));   % missing values: see localMissing
    end
    labels = cols;
    groups = [];
    return
end

fprintf('    table is LONG format: id=%s, time=%s, value=%s',idVar,timeVar,valVar);
if ~isempty(classVar), fprintf(', class=%s',classVar); end
fprintf('\n');

[tsData,labels,groups] = reshapeLong(T.(idVar),T.(timeVar),T.(valVar), ...
    localClassCol(T,classVar));
end

%===============================================================================
function c = localClassCol(T,classVar)
if isempty(classVar), c = []; else, c = T.(classVar); end
end

%===============================================================================
function [tsData,labels,groups] = reshapeLong(ids,times,vals,classes)
% Reshape long-format columns into one vector per series.
%
% Shared by the table path and the numeric-matrix path so the two cannot
% drift apart -- the grouping, time-ordering and validation warnings below
% are subtle enough that maintaining them twice would be a mistake.
if isa(ids,'categorical'), ids = string(ids); end

% Group the rows by series id WITHOUT relying on the third output of unique():
% sorting makes each series contiguous, which is O(n log n) rather than a scan
% per series, and avoids a portability trap (Octave returns an empty inverse
% index for unique(...,'stable') where MATLAB returns a usable one).
[sortedIds,ord] = sort(ids);
if isnumeric(sortedIds) || islogical(sortedIds)
    isNew = [true; sortedIds(2:end) ~= sortedIds(1:end-1)];
else
    sortedIds = string(sortedIds);
    isNew = [true; sortedIds(2:end) ~= sortedIds(1:end-1)];
end
blockStart = find(isNew);
blockStop  = [blockStart(2:end)-1; numel(ord)];
nSeries = numel(blockStart);
uid = sortedIds(blockStart);

% Series come out in sorted id order, not order of first appearance. Long-format
% files are often shuffled, which makes first-appearance arbitrary; sorted order
% is reproducible and matches what a reader expects from the id column.

tsData = cell(nSeries,1);
labels = cell(1,nSeries);
groups = cell(nSeries,1);
nUnsorted = 0; nDup = 0; nIrregular = 0;

for i = 1:nSeries
    sel = ord(blockStart(i):blockStop(i));
    tt = double(times(sel));
    vv = double(vals(sel));

    if any(diff(tt) < 0)                 % long tables are not always in order
        % NOTE the local name. `ord` is the outer grouping index over every row;
        % reusing it here clobbered it, so the second series indexed past the end
        % of a 420-element vector. Present since this function was written, and
        % only reachable when a series has out-of-order timepoints -- which is
        % exactly what a shuffled long-format export looks like.
        nUnsorted = nUnsorted + 1;
        [tt,tOrder] = sort(tt);
        vv = vv(tOrder);
    end
    if any(diff(tt) == 0), nDup = nDup + 1; end
    if numel(tt) > 2
        d = diff(tt);
        if max(abs(d - d(1))) > 1e-9*max(1,abs(d(1))), nIrregular = nIrregular + 1; end
    end

    tsData{i} = vv(:);
    if isnumeric(uid)
        labels{i} = sprintf('series_%g',uid(i));
    else
        labels{i} = char(string(uid(i)));
    end

    if ~isempty(classes)
        cv = classes(sel);
        if isa(cv,'categorical'), cv = string(cv); end
        if isnumeric(cv)
            u = unique(cv);
            groups{i} = cv(1);
        else
            u = unique(string(cv));
            groups{i} = char(string(cv(1)));
        end
        if numel(u) > 1
            warning('HCTSA2_Analyse:mixedClass', ...
                'Series %s carries %u different class labels; using the first.', ...
                labels{i},numel(u));
        end
    end
end

if isempty(classes)
    groups = [];
elseif all(cellfun(@isnumeric,groups))
    groups = cell2mat(groups);
end

len = cellfun(@numel,tsData);
fprintf('    %u series, lengths %u..%u\n',nSeries,min(len),max(len));
if nUnsorted > 0
    fprintf('    NOTE: %u series were not in time order and have been sorted.\n',nUnsorted);
end
if nDup > 0
    warning('HCTSA2_Analyse:dupTime', ...
        '%u series contain duplicate timepoints -- check your data.',nDup);
end
if nIrregular > 0
    warning('HCTSA2_Analyse:irregular', ...
        ['%u series are irregularly sampled. hctsa assumes uniform sampling, so ' ...
         'many features will be misleading on unevenly spaced data.'],nIrregular);
end
if sum(len < 20) > 0
    warning('HCTSA2_Analyse:shortSeries', ...
        '%u series are shorter than 20 points; many features will return NaN.',sum(len < 20));
end
end

%===============================================================================
function tf = isLongFormatMatrix(M)
% Is this numeric matrix one row per observation, rather than one series per
% row or column?
%
% Long format has a narrow shape (3 or 4 columns: id, time, value, class) and a
% first column of repeated integer series identifiers -- few distinct values,
% each appearing many times. A genuine wide matrix of 4 series has 5030 mostly
% distinct floating-point values in its first column, so the two are easy to
% tell apart. The test is deliberately conservative: when in doubt it declines,
% leaving the previous behaviour.
tf = false;
if ~ismatrix(M) || size(M,1) < 20, return; end
nc = size(M,2);
if nc < 3 || nc > 4, return; end
tf = isLongIdColumn(M(:,1));
end

%===============================================================================
function tf = isLongIdColumn(c1)
% A series-id column: integer (or text/categorical) values, at least 2 distinct,
% each repeated -- at most one distinct value per 10 rows.
tf = false;
n = numel(c1);
if n < 20, return; end
if isnumeric(c1) || islogical(c1)
    c1 = double(c1);
    if any(~isfinite(c1)), return; end
    if any(abs(c1 - round(c1)) > 0), return; end      % ids are integers
    u = unique(c1);
else
    u = unique(string(c1));
end
tf = numel(u) >= 2 && numel(u) <= n/10;
end

%===============================================================================
function [tsData,labels,groups] = localFromStructCols(T,opt)
% Reshape a struct of equal-length columns using the same logic as a table.
% Kept separate so the table path is untouched.
vn = fieldnames(T);
idVar    = localPickVar(opt.IdVar,1,vn);
timeVar  = localPickVar(opt.TimeVar,2,vn);
valVar   = localPickVar(opt.ValueVar,3,vn);
classVar = localPickVar(opt.ClassVar,4,vn);
[tsData,labels,groups] = reshapeLong(T.(idVar),T.(timeVar),T.(valVar), ...
    localClassCol(T,classVar));
end

%===============================================================================
function name = localPickVar(userVal,defaultIdx,vn)
% Resolve a column given by name, index, or positional default.
name = '';
if ~isempty(userVal)
    if isnumeric(userVal)
        if userVal >= 1 && userVal <= numel(vn), name = vn{userVal}; end
    else
        userVal = char(userVal);
        if any(strcmp(vn,userVal)), name = userVal;
        else, error('HCTSA2_Analyse:badVar','No column named %s',userVal);
        end
    end
elseif defaultIdx >= 1 && defaultIdx <= numel(vn)
    name = vn{defaultIdx};
end
end
