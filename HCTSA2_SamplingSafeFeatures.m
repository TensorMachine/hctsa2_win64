function T = HCTSA2_SamplingSafeFeatures(dataFile,varargin)
% HCTSA2_SamplingSafeFeatures   Which features survive irregular sampling?
%
%   T = HCTSA2_SamplingSafeFeatures('mydata.mat')
%   T = HCTSA2_SamplingSafeFeatures('mydata.mat','NumSeries',8,'Save','safe.mat')
%
% hctsa never sees your timepoints -- every operation receives a bare value
% vector and assumes the samples are evenly spaced. On irregularly sampled data
% many features still return a number, but that number no longer means what the
% feature name says: "autocorrelation at lag 1" is not a fixed time interval,
% and every spectral feature maps to no real frequency axis.
%
% Rather than guess from operation names, this classifies each feature by
% experiment, on YOUR data:
%
%   ORDER-INDEPENDENT  f(x) == f(shuffle(x)) for every series tested.
%                      The feature depends only on the multiset of values, so
%                      no amount of gap or jitter can affect it. Safe outright.
%
%   GAP-ROBUST         f(x) ~= f(shuffle(x)) -- so it does use the ordering --
%                      but f(x) stays close to f(x with 30% of samples dropped):
%                      the change is small compared with how much the feature
%                      varies between your series (median over the series).
%                      Usable, with the caveat that results describe per-SAMPLE
%                      structure, not per-second structure.
%
%   SAMPLING-DEPENDENT changes materially under decimation. On unevenly spaced
%                      data these should not be trusted.
%
% Deleting random samples is a direct simulation of gaps: it is exactly what an
% irregular series looks like relative to a regular one.
%
%---OUTPUT
% A table, one row per feature: Name, OrderIndependent, GapRobust, the measured
% discrepancies, and a Recommendation. To compute only the usable features:
%   keep = T.Name(T.OrderIndependent | T.GapRobust);
%   R = HCTSA2_Analyse('mydata.mat','Features',keep);
% ('Features' is applied before extraction, so the rest is never computed.)
%
%---COST
% Three feature computations per sampled series. With the default 5 series and
% the full library that is roughly 5-10 minutes; 'Parallel',true spreads each
% computation over a parallel pool. Use 'FeatureSet','catch22' for a quick look
% at the method.

p = inputParser;
addParameter(p,'NumSeries',5,@isnumeric);
addParameter(p,'DropFraction',0.3,@isnumeric);
addParameter(p,'Tolerance',1e-8,@isnumeric);
addParameter(p,'GapTolerance',0.1,@isnumeric);
addParameter(p,'Parallel',false,@(x) islogical(x) || isnumeric(x));
addParameter(p,'FeatureSet','hctsa',@ischar);
addParameter(p,'Save','',@ischar);
addParameter(p,'Seed',1,@isnumeric);
% passed through to the ingester
addParameter(p,'DataVar','',@ischar);
addParameter(p,'IdVar',[]); addParameter(p,'TimeVar',[]);
addParameter(p,'ValueVar',[]); addParameter(p,'ClassVar',[]);
parse(p,varargin{:});
opt = p.Results;

fprintf('\n========= feature robustness to irregular sampling =========\n');

%-------------------------------------------------------------------------------
%% Get some representative series
%-------------------------------------------------------------------------------
ing = opt;
ing.Groups = []; ing.Labels = {};
tsData = HCTSA2_Ingest(dataFile,ing);
nAll = numel(tsData);
rng(opt.Seed);
pick = randperm(nAll,min(opt.NumSeries,nAll));
fprintf('Testing on %u of %u series.\n',numel(pick),nAll);

%-------------------------------------------------------------------------------
%% Build the operation list once
%-------------------------------------------------------------------------------
tmpDir = tempname; mkdir(tmpDir);
cleanupObj = onCleanup(@() rmdir(tmpDir,'s'));
inpFile = fullfile(tmpDir,'INP.mat');
timeSeriesData = tsData(pick); %#ok<NASGU>
labels = arrayfun(@(i) sprintf('probe_%u',i),1:numel(pick),'UniformOutput',false); %#ok<NASGU>
keywords = repmat({'probe'},1,numel(pick)); %#ok<NASGU>
save(inpFile,'timeSeriesData','labels','keywords');
probeFile = fullfile(tmpDir,'HCTSA_probe.mat');
TS_Init(inpFile,opt.FeatureSet,[false,false,false],probeFile);
L = load(probeFile,'Operations','MasterOperations');
Operations = L.Operations;
MasterOperations = L.MasterOperations;
nOps = height(Operations);
fprintf('Evaluating %u features.\n\n',nOps);

%-------------------------------------------------------------------------------
%% Compute each series three ways
%-------------------------------------------------------------------------------
nS = numel(pick);
F0 = nan(nS,nOps); Fp = nan(nS,nOps); Fd = nan(nS,nOps);
for i = 1:nS
    x = tsData{pick(i)};
    x = x(isfinite(x));
    rng(opt.Seed + i);
    xPerm = x(randperm(numel(x)));
    keep = sort(randperm(numel(x),round((1-opt.DropFraction)*numel(x))));
    xDrop = x(keep);

    fprintf('  series %u/%u (n=%u) ',i,nS,numel(x));
    F0(i,:) = localCompute(x,Operations,MasterOperations,opt.Parallel);      fprintf('.');
    Fp(i,:) = localCompute(xPerm,Operations,MasterOperations,opt.Parallel);   fprintf('.');
    Fd(i,:) = localCompute(xDrop,Operations,MasterOperations,opt.Parallel);   fprintf('. done\n');
end

%-------------------------------------------------------------------------------
%% Classify
%-------------------------------------------------------------------------------
scale = @(a,b) max(max(abs(a),abs(b)),1);
permDiff = abs(F0-Fp)./scale(F0,Fp);
% A dropped-sample change counts relative to how much the feature varies
% between the series tested (HCTSA2 2.2.1; was relative to its magnitude, which
% for values below 1 made the tolerance absolute). A feature equal on every
% series tested keeps the magnitude scale.
spread = std(F0,0,1,'omitnan');
dropScale = repmat(spread,nS,1);
noSpread = ~(dropScale > 0);
magScale = scale(F0,Fd);
dropScale(noSpread) = magScale(noSpread);
dropDiff = abs(F0-Fd)./dropScale;

% a feature counts only where it computed on every variant
valid = isfinite(F0) & isfinite(Fp) & isfinite(Fd);
nValid = sum(valid,1)';

permMax = nan(nOps,1); dropMed = nan(nOps,1);
for j = 1:nOps
    v = valid(:,j);
    if any(v)
        permMax(j) = max(permDiff(v,j));
        dropMed(j) = median(dropDiff(v,j));
    end
end

orderIndependent = (permMax <= opt.Tolerance) & (nValid > 0);
gapRobust = ~orderIndependent & (dropMed <= opt.GapTolerance) & (nValid > 0);

rec = repmat("sampling-dependent",nOps,1);
rec(gapRobust) = "gap-robust (per-sample interpretation)";
rec(orderIndependent) = "safe (order-independent)";
rec(nValid == 0) = "not evaluated";

T = table();
T.Name = string(Operations.Name);
try T.Keywords = string(Operations.Keywords); catch, end
T.OrderIndependent = orderIndependent;
T.GapRobust = gapRobust;
T.PermDiff = permMax;
T.DropDiff = dropMed;
T.nSeriesEvaluated = nValid;
T.Recommendation = rec;
T = sortrows(T,{'OrderIndependent','GapRobust','DropDiff'},{'descend','descend','ascend'});

%-------------------------------------------------------------------------------
%% Report
%-------------------------------------------------------------------------------
nOI = sum(orderIndependent); nGR = sum(gapRobust);
nSD = sum(strcmp(rec,"sampling-dependent")); nNE = sum(nValid==0);
fprintf('\n-----------------------------------------------------------\n');
fprintf('  safe (order-independent) : %5u  (%.0f%%)\n',nOI,100*nOI/nOps);
fprintf('  gap-robust               : %5u  (%.0f%%)\n',nGR,100*nGR/nOps);
fprintf('  sampling-dependent       : %5u  (%.0f%%)\n',nSD,100*nSD/nOps);
if nNE > 0
    fprintf('  not evaluated            : %5u  (failed on every variant)\n',nNE);
end
fprintf('-----------------------------------------------------------\n');
fprintf('Order-independent features cannot be affected by gaps at all.\n');
fprintf('Gap-robust ones use the ordering, so read them as per-SAMPLE\n');
fprintf('structure rather than per-second.\n');
if nOI + nGR < 0.1*nOps
    fprintf(2,'\nNote: very few features survive. Irregular sampling is a real\n');
    fprintf(2,'constraint on what hctsa can tell you about this data.\n');
end
fprintf('\n');

if ~isempty(opt.Save)
    save(opt.Save,'T');
    fprintf('Saved to %s\n\n',opt.Save);
end
end

%===============================================================================
function fv = localCompute(x,Operations,MasterOperations,doParallel)
% One feature vector, tolerating individual operation failures. Any output that
% is not a normal value (TS_Quality code > 0: error, NaN/Inf, complex, empty,
% link error) is NaN, so it cannot pass as a stable feature.
x = (x - mean(x))/std(x);
try
    [fv,~,q] = TS_CalculateFeatureVector(x,logical(doParallel),Operations, ...
                    MasterOperations,false,'fast');
    fv(q > 0) = NaN;
catch
    fv = nan(height(Operations),1);
end
fv = real(fv(:)');
end
