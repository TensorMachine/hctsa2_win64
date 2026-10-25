function [F,info] = HCTSA2_ComputeFeatures(x,features,varargin)
% HCTSA2_ComputeFeatures   Compute chosen hctsa features of one time series.
%
%   F = HCTSA2_ComputeFeatures(x,{'mean','SB_MotifTwo_mean_hhh'})
%   F = HCTSA2_ComputeFeatures(x,[1011 4809])           % the same, by ID (3.0.x)
%   [F,info] = HCTSA2_ComputeFeatures(x,features,'FeatureSet','catch22')
%   F = HCTSA2_ComputeFeatures(x,features,'NormalizeLike','myRun\HCTSA_N.mat')
%
% For a feature-extraction pipeline: give it one time series and the features
% chosen from an HCTSA2_Analyse run (from FeatureAUC, say), and it returns a
% struct with one field per feature, named after the feature:
%
%   sel = R.FeatureAUC.Name(1:20);           % R = HCTSA2_Analyse(...)
%   F = HCTSA2_ComputeFeatures(x,sel);
%   F.mean, F.SB_MotifTwo_mean_hhh, ...
%   v = cell2mat(struct2cell(F))';           % a row vector, in the order of sel
%
% Each master operation the features need is computed once (by
% TS_CalculateFeatureVector, as TS_Compute does), and only the requested features
% are kept. The values are those HCTSA2_Analyse stores for the same series in its
% HCTSA.mat: identical on the same computer and MATLAB release, equal to rounding
% on another.
%
%---INPUTS
% x          one time series: a numeric (or logical) vector. As in HCTSA2_Analyse,
%            NaN/Inf padding at either end is trimmed. With NaN/Inf left inside it,
%            or nothing left, every feature is NaN and a warning is given (so a
%            pipeline keeps going). A series longer than 50,000 samples, which
%            HCTSA2_Analyse leaves out, is computed, with a warning.
% features   which features, either
%            - names (char, string array, or a cell array of them), as in
%              FeatureTable.Name or FeatureAUC.Name, or their field names below
%              (FeatureTable.Column);
%            - or IDs (numbers), as in FeatureTable.ID.
%            Prefer names in a pipeline you will keep: an ID is the feature's row
%            in the library's list (FeatureSets/INP_ops_hctsa.txt), and IDs move
%            when the library changes -- only 14% kept their feature across
%            2.5.0's upstream merge, and 2% across 3.0.0's -- with no error. A name
%            that no longer exists is an error (docs/HCTSA2_FeatureRenames_3.0.0.csv
%            maps 2.5.x names to 3.0.0's).
%
%---OPTIONS (name/value)
% 'FeatureSet'  the library the names and IDs belong to: 'hctsa' (default),
%               'catch22' or 'catch24' -- the 'FeatureSet' of the HCTSA2_Analyse run
% 'Parallel'    spread the master operations over a parallel pool (default false).
%               Needs the Parallel Computing Toolbox (without it the computation is
%               serial) and opens a pool if none is open; its workers need hctsa on
%               their path. Leave it false inside your own parfor.
% 'NormalizeLike'  return the values normalised exactly as a run normalised its
%               data (HCTSA2 3.0.0): the run's HCTSA_N.mat (or its output folder, or
%               the normalizationInfo struct it holds). Since 3.0.0, TS_Normalize
%               stores each feature's fitted normalisation there, and it is applied
%               to this series with the run's parameters (BF_ApplyTransform; nothing
%               is fitted to the series). A series of the run gets its HCTSA_N.mat
%               values; a new series the values it would have had, which for the
%               scaled sigmoids can lie a little outside [0,1]. Every requested
%               feature must be one the run kept (an error names any it filtered
%               out); a run made before 3.0.0, or with another feature library, is an
%               error. NaN stays NaN ('NaNToZero' is not applied). Default: raw values.
%
%---OUTPUTS
% F      a struct with one field per requested feature, in the requested order,
%        named after the feature. Where a name is not a valid field name (72 of
%        the 6,392 hctsa names: 66 contain '-', 6 are longer than 63 characters)
%        the field is the name HCTSA2_Analyse gives the feature's SeriesTable
%        column (FeatureTable.Column): other characters become '_', and a long
%        name keeps its first 24 and last 38 characters, e.g.
%        'firstMax_mi-gaussian' -> firstMax_mi_gaussian. (Like SeriesTable's,
%        the fields avoid the names Name, Group and Length.)
%        The values are raw by default, as FeatureAUC uses them, or normalised
%        like a run's ('NormalizeLike'). A feature that fails, or gives NaN, Inf, a
%        complex number or nothing, is NaN, as in HCTSA2_Analyse; if any feature
%        fails with an error, a warning says which, with the first error message
%        (identifier 'HCTSA2_ComputeFeatures:errors': turn it off with warning
%        if you expect it).
% info   a table with one row per requested feature: ID, Name, Field, MasterID
%        (features with the same MasterID come from one computation) and Quality
%        (0 a value; 1 an error; 2 NaN; 3 Inf; 4 -Inf; 5 complex; 6 an empty
%        output; 7 the output has no such field; NaN not computed, because the
%        series could not be used); with 'NormalizeLike', also RawValue.
%
% Run startup first. Features whose operations need compiled components (TISEAN,
% the KSG MEX) need them installed, as for HCTSA2_Analyse; without them they fail
% (Quality 1, and the warning above). A problem that stops the whole computation
% (hctsa not on the path, for example) is an error. The library is read on the
% first call and kept ('clear HCTSA2_ComputeFeatures' reads it again). Each call
% leaves the warning settings as it found them, and prints nothing.
%
% HCTSA2 2.5.1; 'NormalizeLike' and the first error message in the warning, 3.0.0.

%-------------------------------------------------------------------------------
%% Check the inputs
%-------------------------------------------------------------------------------
narginchk(2,Inf);
p = inputParser;
addParameter(p,'FeatureSet','hctsa',@(s) ischar(s) || (isstring(s) && isscalar(s)));
addParameter(p,'Parallel',false,@(v) (islogical(v) || isnumeric(v)) && isscalar(v));
addParameter(p,'NormalizeLike',[],@(v) isempty(v) || ischar(v) || (isstring(v) && isscalar(v)) ...
    || isstruct(v));
parse(p,varargin{:});
featureSet = lower(char(p.Results.FeatureSet));
if ~ismember(featureSet,{'hctsa','catch22','catch24'})
    error('HCTSA2_ComputeFeatures:featureSet', ...
        '''FeatureSet'' must be ''hctsa'', ''catch22'' or ''catch24''.');
end
doParallel = logical(p.Results.Parallel);
if ~((isnumeric(x) || islogical(x)) && isreal(x) && (isvector(x) || isempty(x)))
    error('HCTSA2_ComputeFeatures:data', ...
        'x must be one time series: a real numeric vector (one series per call).');
end

lib = localLibrary(featureSet);
idx = localResolve(features,lib.ID,lib.Name,lib.Field,featureSet);
Ops = lib.Operations(idx,:);
numFeatures = numel(idx);
% The run's normalisation of these features (checked before any computation)
normT = [];
if ~isempty(p.Results.NormalizeLike)
    normT = localTransform(p.Results.NormalizeLike,lib,idx,featureSet);
end

%-------------------------------------------------------------------------------
%% The series, prepared as HCTSA2_Analyse's ingestion prepares it
%-------------------------------------------------------------------------------
x = double(x(:));
first = find(isfinite(x),1);
last = find(isfinite(x),1,'last');
if isempty(first)
    x = zeros(0,1);
else
    x = x(first:last);
end

%-------------------------------------------------------------------------------
%% Compute
%-------------------------------------------------------------------------------
values = nan(numFeatures,1);
quality = nan(numFeatures,1);
if isempty(x)
    warning('HCTSA2_ComputeFeatures:empty', ...
        'The time series has no finite values: every feature is NaN.');
elseif ~all(isfinite(x))
    warning('HCTSA2_ComputeFeatures:nonFinite', ...
        ['The time series has %u NaN/Inf values inside it (not padding at an end): ' ...
         'every feature is NaN.'],sum(~isfinite(x)));
else
    if numel(x) > 50000
        warning('HCTSA2_ComputeFeatures:long', ['The time series has %u samples; ' ...
            'HCTSA2_Analyse leaves out series longer than 50,000, and some operations ' ...
            'take a long time on them.'],numel(x));
    end
    if all(x == x(1))
        warning('HCTSA2_ComputeFeatures:constant', ...
            'The time series is constant: most features will be NaN.');
    end
    % Operations switch warnings on and off; their messages are captured. The only
    % one left on screen would be TS_CalculateFeatureVector's own about a constant
    % series (given above, with an identifier). Restored on the way out.
    warningState = warning;
    restoreWarnings = onCleanup(@() warning(warningState));
    warning('off','all');
    % (3.0.0: the fourth output holds each failed master operation and its message)
    [fv,~,q,masterErrors] = TS_CalculateFeatureVector(x,doParallel,Ops,lib.MasterOperations, ...
                    false,'fast',false);
    clear('restoreWarnings');
    fv(q > 0) = NaN;          % error, NaN, Inf, -Inf, complex, empty: NaN, as in HCTSA2_Analyse
    values = real(fv(:));
    quality = q(:);
    failed = find(quality == 1);
    if ~isempty(failed)
        firstError = '';
        if ~isempty(masterErrors)
            firstError = sprintf(' (%s: %s)',masterErrors{1,1}, ...
                strrep(strtrim(masterErrors{1,2}),newline,' | '));
        end
        warning('HCTSA2_ComputeFeatures:errors','%u of the %u features failed with an error: %s%s%s', ...
            numel(failed),numFeatures,strjoin(lib.Name(idx(failed(1:min(5,end))))',', '), ...
            localIf(numel(failed) > 5,', ...',''),firstError);
    end
end

rawValues = values;
if ~isempty(normT)
    % The run's stored parameters, applied as TS_Normalize applied them (no imputation)
    values = BF_ApplyTransform(values',normT,lib.ID(idx),false)';
end

F = cell2struct(num2cell(values),lib.Field(idx),1);
if nargout > 1
    info = table(lib.ID(idx),string(lib.Name(idx)),string(lib.Field(idx)),Ops.MasterID, ...
        quality,'VariableNames',{'ID','Name','Field','MasterID','Quality'});
    if ~isempty(normT)
        info.RawValue = rawValues;
    end
end
end

%===============================================================================
function T = localTransform(src,lib,idx,featureSet)
% The normalisation a run fitted (TS_Normalize, upstream v3.0.0's normalizationInfo
% .transform), cut down to the requested features, in the requested order.
if isstruct(src)
    ni = src;
    where = 'the given struct';
else
    f = char(src);
    if isfolder(f), f = fullfile(f,'HCTSA_N.mat'); end   % an HCTSA2_Analyse output folder
    if ~isfile(f)
        error('HCTSA2_ComputeFeatures:normalizeLike','''NormalizeLike'': %s not found.',f);
    end
    vars = whos('-file',f);
    if ~ismember('normalizationInfo',{vars.name})
        error('HCTSA2_ComputeFeatures:normalizeLike', ...
            '''NormalizeLike'': %s holds no normalizationInfo (give a run''s HCTSA_N.mat).',f);
    end
    S = load(f,'normalizationInfo');
    ni = S.normalizationInfo;
    where = f;
end
if isfield(ni,'transform') && isstruct(ni.transform)
    T = ni.transform;
elseif isfield(ni,'featureParams') && isfield(ni,'spec')
    T = ni;                                              % a transform itself
else
    error('HCTSA2_ComputeFeatures:noTransform', ['''NormalizeLike'': %s does not store ' ...
        'the fitted normalisation (a run made before HCTSA2 3.0.0): run the analysis ' ...
        'again with 3.0.0.'],where);
end
fp = T.featureParams;
[inRun,row] = ismember(lib.ID(idx),fp.ID);
sameName = true(size(inRun));
sameName(inRun) = strcmp(cellstr(string(fp.Name(row(inRun)))),lib.Name(idx(inRun)));
if ~all(sameName)
    k = find(~sameName,1);
    error('HCTSA2_ComputeFeatures:otherLibrary', ['''NormalizeLike'': the run''s feature ' ...
        '%u is %s, but this %s library''s is %s: the run was made with another feature ' ...
        'library.'],lib.ID(idx(k)),char(string(fp.Name(row(k)))),featureSet,lib.Name{idx(k)});
end
if ~all(inRun)
    bad = lib.Name(idx(~inRun));
    error('HCTSA2_ComputeFeatures:notNormalized', ['''NormalizeLike'': %u requested ' ...
        'feature(s) were filtered out by the run''s normalisation: %s%s'],numel(bad), ...
        strjoin(bad(1:min(10,end))',', '),localIf(numel(bad) > 10,', ...',''));
end
T.featureParams = fp(row,:);
if isfield(T,'trainMedian') && numel(T.trainMedian) == height(fp)
    T.trainMedian = T.trainMedian(row);
end
end

%===============================================================================
function lib = localLibrary(featureSet)
% The feature set's operations and master operations, as TS_Init builds them
% (so the IDs are a run's), with each master's code compiled once, and every
% feature's field name. Read once per session.
persistent cache
if isempty(cache), cache = struct(); end
if ~isfield(cache,featureSet)
    [mopsFile,opsFile] = TS_GiveMeInputFiles(featureSet);
    % (quietly: TS_ReadInputFile always prints how many it read)
    evalc(['Operations = TS_ReadInputFile(''ops'',opsFile,false); ' ...
           'MasterOperations = TS_ReadInputFile(''mops'',mopsFile,false); ' ...
           '[Operations,MasterOperations] = TS_LinkOperationsWithMasters(Operations,MasterOperations);']);
    % TISEAN temporary files left by an earlier MATLAB that was stopped mid-run,
    % as TS_Compute sweeps them (only sets untouched for a day)
    evalc('BF_SweepStaleTempFiles();');
    MasterOperations.Fn = cell(height(MasterOperations),1);
    for i = 1:height(MasterOperations)
        MasterOperations.Fn{i} = str2func(['@(x,x_z) ',MasterOperations.Code{i}]);
    end
    L = struct();
    L.Operations = Operations;
    L.MasterOperations = MasterOperations;
    L.ID = Operations.ID;
    L.Name = cellstr(Operations.Name);
    L.Field = localFieldNames(L.Name,{'Name','Group','Length'});
    cache.(featureSet) = L;
end
lib = cache.(featureSet);
end

%===============================================================================
function idx = localResolve(req,ids,names,fields,featureSet)
% Rows of the library for the requested names, field names or IDs, in the
% requested order. (Plain arrays: the Octave check runs it.)
if isnumeric(req)
    req = req(:);
    if any(~isfinite(req) | req ~= round(req))
        error('HCTSA2_ComputeFeatures:features','Feature IDs must be whole numbers.');
    end
    [found,idx] = ismember(req,ids);
    shown = arrayfun(@(v) sprintf('%g',v),req,'UniformOutput',false);
elseif ischar(req) || isstring(req) || (iscell(req) && all(cellfun(@(c) ischar(c) || ...
        (isstring(c) && isscalar(c)),req(:))))
    if iscell(req)
        req = cellfun(@char,req,'UniformOutput',false);
    end
    req = cellstr(req);
    req = req(:);
    [found,idx] = ismember(req,names);
    [isField,idxField] = ismember(req,fields);
    useField = ~found & isField;
    idx(useField) = idxField(useField);
    found = found | isField;
    shown = req;
else
    error('HCTSA2_ComputeFeatures:features', ...
        'Give the features as names (text) or IDs (numbers).');
end
if isempty(req)
    error('HCTSA2_ComputeFeatures:features','No features were requested.');
end
if ~all(found)
    bad = shown(~found);
    error('HCTSA2_ComputeFeatures:unknown','%u requested feature(s) not in the %s library: %s%s', ...
        numel(bad),featureSet,strjoin(bad(1:min(10,end))',', '),localIf(numel(bad) > 10,', ...',''));
end
[~,first] = unique(idx,'stable');
if numel(first) < numel(idx)
    again = setdiff(1:numel(idx),first);
    error('HCTSA2_ComputeFeatures:duplicate','Feature(s) requested more than once: %s', ...
        strjoin(names(idx(again))',', '));
end
idx = idx(:);
end

%===============================================================================
function f = localFieldNames(names,reserved)
% The rule HCTSA2_Analyse uses for SeriesTable's column names (localColumnNames,
% with SeriesTable's other columns reserved), applied to the whole library so
% that a feature's field never depends on which others are requested: valid
% names of at most 63 characters; other characters become '_', and a long name
% keeps its first 24 and last 38 characters (the end tells the MF_GP_* variants
% apart). Unique over the whole library.
f = regexprep(cellstr(names(:)),'[^A-Za-z0-9_]','_');
for i = 1:numel(f)
    if isempty(regexp(f{i},'^[A-Za-z]','once')), f{i} = ['x' f{i}]; end
    if numel(f{i}) > 63
        f{i} = [f{i}(1:24) '_' f{i}(end-37:end)];
    end
end
f = matlab.lang.makeUniqueStrings(f,reserved,63);
end

%===============================================================================
function out = localIf(cond,a,b)
if cond, out = a; else, out = b; end
end
