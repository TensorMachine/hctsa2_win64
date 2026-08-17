function [tsData,labels,groups] = HCTSA2_Ingest(dataFile,opt)
% (Type tests use isa(x,'categorical') and isa(x,'string') rather than
% iscategorical/isstring. The two forms are equivalent in MATLAB, but only isa
% exists in Octave -- which is where this code is developed and smoke-tested, so
% the equivalent form keeps that possible without changing MATLAB behaviour.)
% HCTSA2_Ingest   Read time series out of a .mat in whatever shape it arrives.
%
%   [tsData,labels,groups] = HCTSA2_Ingest('mydata.mat',opt)
%
% OPT is a struct that may carry: DataVar, Labels, Groups, IdVar, TimeVar,
% ValueVar, ClassVar. Missing fields are filled with defaults, so a bare
% struct() is fine.
%
% Handles four layouts:
%   - LONG-FORMAT TABLE  (id, time, value, class) -- one row per observation
%   - wide table         (one column per series)
%   - numeric matrix     (rows = series)
%   - cell array         (one vector per series)
%
% Shared by HCTSA2_Analyse and HCTSA2_SamplingSafeFeatures so both read data
% identically.

defaults = struct('DataVar','','Labels',{{}},'Groups',[], ...
                  'IdVar',[],'TimeVar',[],'ValueVar',[],'ClassVar',[]);
fn = fieldnames(defaults);
for i = 1:numel(fn)
    if ~isfield(opt,fn{i}) || (isempty(opt.(fn{i})) && ~ischar(defaults.(fn{i})))
        if ~isfield(opt,fn{i}), opt.(fn{i}) = defaults.(fn{i}); end
    end
end

[tsData,labels,groups] = localIngest(dataFile,opt);
end

%===============================================================================
function [tsData,labels,groups] = localIngest(dataFile,opt)
% Find the time series in an arbitrary .mat. Handles four layouts:
%   - LONG-FORMAT TABLE  (id, time, value, class) -- one row per observation
%   - wide table         (one column per series)
%   - numeric matrix     (rows = series)
%   - cell array         (one vector per series)

if exist(dataFile,'file') ~= 2
    error('HCTSA2_Analyse:noFile','Could not find %s',dataFile);
end
D = load(dataFile);
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
elseif isnumeric(raw) && isLongFormatMatrix(raw)
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

if ~isempty(opt.Labels), labels = opt.Labels; end
if isempty(labels)
    for i = 1:numel(fn)
        v = D.(fn{i});
        if iscellstr(v) && numel(v) == n && ~strcmpi(fn{i},'keywords') %#ok<ISCLSTR>
            labels = v; break
        end
    end
end
if isempty(labels)
    labels = arrayfun(@(i) sprintf('Series_%u',i),1:n,'UniformOutput',false);
end

if ~isempty(opt.Groups)
    groups = opt.Groups;
    if ischar(groups)
        if isfield(D,groups), groups = D.(groups);
        else, error('HCTSA2_Analyse:noGroupVar','No variable %s in %s',groups,dataFile);
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
end

%===============================================================================
function [tsData,labels,groups] = localFromTable(T,opt)
% Reshape a table into one vector per series.
%
% Long format is detected by an ID column that repeats. Columns are taken
% POSITIONALLY by default -- (id, time, value, class) -- because column names
% vary between datasets; override with IdVar/TimeVar/ValueVar/ClassVar, each of
% which accepts a name or a column index.

vn = T.Properties.VariableNames;
idVar    = localPickVar(opt.IdVar,1,vn);
timeVar  = localPickVar(opt.TimeVar,2,vn);
valVar   = localPickVar(opt.ValueVar,3,vn);
classVar = localPickVar(opt.ClassVar,4,vn);

isLong = ~isempty(idVar) && height(T) > numel(unique(T.(idVar)));

if ~isLong
    fprintf('    table looks WIDE (one column per series)\n');
    num = varfun(@isnumeric,T,'OutputFormat','uniform');
    cols = vn(num);
    tsData = cell(numel(cols),1);
    for i = 1:numel(cols)
        v = double(T.(cols{i}));
        tsData{i} = v(isfinite(v));
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
            groups{i} = u(1);
        else
            u = unique(string(cv));
            groups{i} = char(u(1));
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
c1 = M(:,1);
if any(~isfinite(c1)), return; end
if any(abs(c1 - round(c1)) > 0), return; end          % ids are integers
u = unique(c1);
if numel(u) < 2 || numel(u) > size(M,1)/10, return; end % repeated many times
tf = true;
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
        if any(strcmp(vn,userVal)), name = userVal;
        else, error('HCTSA2_Analyse:badVar','No column named %s',userVal);
        end
    end
elseif defaultIdx >= 1 && defaultIdx <= numel(vn)
    name = vn{defaultIdx};
end
end

