function [report,ok] = HCTSA2_ValidateKnownSignals(verbose)
% HCTSA2_ValidateKnownSignals   Do the features recover properties we know?
%
%   report = HCTSA2_ValidateKnownSignals
%   [report,ok] = HCTSA2_ValidateKnownSignals(false)
%
% The reference export predates the merge: 766 of its features no longer exist,
% 757 current features never appeared in it, and 141 operations were redefined.
% For those, comparing against the reference is not possible even in principle.
%
% This validates differently: it feeds the library signals whose dynamical
% properties are known analytically or from the literature, and checks that the
% features recover them. A feature that returns a plausible-looking number on
% EEG but reports a correlation dimension of 5 for a sine wave is broken, and no
% amount of reference comparison would reveal it.
%
%---SIGNALS AND THEIR KNOWN PROPERTIES
%
%   sine        period 50, amplitude 1
%               largest Lyapunov 0, autocorrelation first zero at a quarter
%               period (~12-13), perfectly predictable
%   logistic    x(n+1) = 4x(n)(1-x(n)), fully chaotic
%               largest Lyapunov = ln(2) = 0.6931 exactly
%   henon       a=1.4, b=0.3, the classic attractor
%               largest Lyapunov ~0.419
%   whitenoise  iid Gaussian
%               no deterministic structure: autocorrelation ~0 beyond lag 0
%   ar1         x(n) = 0.8x(n-1) + e(n)
%               autocorrelation at lag k is exactly 0.8^k
%
% NOT checked (HCTSA2 2.2.1: the header used to promise these): correlation
% dimensions -- 1 for the sine and the logistic map, ~1.22 for Henon, "not
% converging low" for white noise. NL_BoxCountEntropyRate (NL_BoxCorrDim before
% 2.3.0) reports statistics of the Renyi entropy curves, not one dimension
% estimate, so it is only run as a smoke test.
%
% Each check states the expected value and the tolerance it is judged against,
% so a failure says what was expected rather than only that something differed.
%
% A check is SKIPPED, not failed, only for a named gap in the environment: the
% TISEAN executables are not built (error hctsa:tisean:notInstalled), or MATLAB's
% histcounts2 is missing (Octave). Any other error is a FAIL that carries its
% message, and so is a missing operation or result field: a feature that now
% throws is the thing this suite exists to find. A skipped check is not evidence
% of anything either way; the summary counts them separately and ok is false when
% more than three were skipped or nothing was checked at all.

if nargin < 1, verbose = true; end
% (HCTSA2 2.2.1: rng(1) re-seeds the global stream; leave it as the caller had it.
% The checks are in a local function because this one holds the onCleanup, and
% Octave does not restore from an onCleanup in a function that has nested ones.)
rngState = rng;
restoreRng = onCleanup(@() rng(rngState)); %#ok<NASGU>
[report,ok] = runChecks(verbose);
end

%-------------------------------------------------------------------------------
function [report,ok] = runChecks(verbose)
report = struct('signal',{},'check',{},'expected',{},'got',{},'pass',{},'skipped',{});

if verbose
    fprintf('\n=========== known-signal validation ===========\n');
end

%-------------------------------------------------------------------------------
%% Build the signals
%-------------------------------------------------------------------------------
N = 4000;
rng(1);

t = (1:N)';
sig.sine = sin(2*pi*t/50);

x = zeros(N,1); x(1) = 0.3;
for i = 2:N, x(i) = 4*x(i-1)*(1-x(i-1)); end
sig.logistic = x;

hx = zeros(N,1); hy = zeros(N,1); hx(1) = 0.1; hy(1) = 0.1;
for i = 2:N
    hx(i) = 1 - 1.4*hx(i-1)^2 + hy(i-1);
    hy(i) = 0.3*hx(i-1);
end
sig.henon = hx;

sig.whitenoise = randn(N,1);

e = randn(N,1); a = zeros(N,1);
for i = 2:N, a(i) = 0.8*a(i-1) + e(i); end
sig.ar1 = a;

names = fieldnames(sig);
for i = 1:numel(names)
    z = sig.(names{i});
    sig.(names{i}) = (z - mean(z))/std(z);   % the library z-scores its input
end

%-------------------------------------------------------------------------------
%% Checks with known answers
%-------------------------------------------------------------------------------

% --- AR(1): autocorrelation is exactly 0.8^k -------------------------------
y = sig.ar1;
for k = [1 2 3]
    try
        got = CO_AutoCorr(y,k,'Fourier');
        add('ar1',sprintf('autocorrelation at lag %u',k),0.8^k,got,abs(got-0.8^k) < 0.08);
    catch err
        onError('ar1',sprintf('autocorrelation at lag %u',k),err);
    end
end

% --- sine: autocorrelation first zero at a quarter period ------------------
y = sig.sine;
try
    got = CO_FirstCrossing(y,'ac',0,'discrete');
    add('sine','first zero of autocorrelation (period/4 = 12.5)',12.5,got,abs(got-12.5) <= 1.5);
catch err
    onError('sine','first zero of autocorrelation',err);
end

% --- sine: a linear model should predict it almost perfectly ---------------
try
    r = MF_steps_ahead(y,'ar',5,2);
    got = r.stde_h1;          % was rmserr_1 before upstream's rename (7f498720)
    add('sine','1-step AR prediction error (near zero)',0,got,got < 0.05);
catch err
    onError('sine','1-step AR prediction error',err);
end

% --- white noise: autocorrelation vanishes beyond lag 0 --------------------
y = sig.whitenoise;
try
    got = max(abs(arrayfun(@(k) CO_AutoCorr(y,k,'Fourier'),1:10)));
    add('whitenoise','max |autocorrelation|, lags 1-10 (~0)',0,got,got < 0.1);
catch err
    onError('whitenoise','max |autocorrelation| lags 1-10',err);
end

% --- white noise vs sine: entropy ordering ---------------------------------
try
    eN = EN_ApEn(sig.whitenoise,1,0.2);
    eS = EN_ApEn(sig.sine,1,0.2);
    add('ordering','ApEn(white noise) > ApEn(sine)',1,eN > eS,eN > eS);
catch err
    onError('ordering','ApEn(white noise) > ApEn(sine)',err);
end

% --- Lyapunov exponent: the strongest check in this suite ------------------
% The logistic map at r = 4 has largest Lyapunov exponent exactly ln 2 =
% 0.693147, and a sine wave (a periodic orbit) has exactly 0.
%
% This began as an ordering check only, on the assumption that the estimator's
% absolute output would depend on its fitting range. It does not: measured on
% MATLAB with TISEAN's lyap_r, the logistic map returns **0.692** -- 0.17% from
% the analytic value -- and the sine returns 0.0007. That is a genuine
% quantitative validation of the whole path (TISEAN binary, embedding, scaling
% range selection, BF_fit regression), so the value is now checked directly.
%
% Tolerance 0.05 absolute: wide enough for estimator variance across embedding
% choices, narrow enough that a broken pipeline cannot pass.
%
% The call is the shipped operation NL_LargestLyap_n1_ac1e30_ac1e1_3_1_4 of
% FeatureSets/INP_mops_hctsa.txt: NL_LargestLyap(y, Nref, maxtstep, past, NNR,
% embedParams) with maxtstep {'ac1e',30} and past {'ac1e',1} (2.3.0; until 2.2.1
% {'ac',20} and {'ac',1}; the values here agree within 0.001). (HCTSA2 2.2.1: an
% earlier version passed the proportions 0.1 and 0.01 there, a parameter set the
% library does not ship; both give the analytic values -- logistic 0.692,
% sine 0.0007 against 0.6926, 0.0004 -- but this one is what users get.) An
% earlier version also passed 'cao' as `past`, which made a 1x3 array reach an
% && and the failure was reported as a missing TISEAN binary, which it was not.
if exist('NL_LargestLyap','file') ~= 2
    add('logistic','NL_LargestLyap exists',1,0,false);   % a missing operation is a FAIL
else
    try
        % Name the field explicitly. An earlier version used a "first finite
        % scalar field" helper, which returned `maxp` or a crossing count --
        % the value came back as exactly 0 and the ordering check failed for a
        % reason that had nothing to do with the library. The Lyapunov exponent
        % is the SLOPE of the fitted log-divergence curve, which this operation
        % reports as vse_gradient (a penalised linear regression over the
        % scaling range, with ve_gradient as the fixed-endpoint variant).
        lyapOf = @(z) namedField(NL_LargestLyap(z,-1,{'ac1e',30},{'ac1e',1},3,{1,4}), ...
                                 {'vse_gradient','ve_gradient'});
        Lchaos = lyapOf(sig.logistic);
        Lsine  = lyapOf(sig.sine);
        add('logistic','largest Lyapunov exponent (ln 2)',log(2),Lchaos, ...
            abs(Lchaos - log(2)) < 0.05);
        add('sine','largest Lyapunov exponent (periodic, 0)',0,Lsine, ...
            abs(Lsine) < 0.05);
        add('ordering','Lyapunov(logistic) > Lyapunov(sine)',1,Lchaos > Lsine,Lchaos > Lsine);
        % The Henon map's largest Lyapunov exponent is 0.419 (Henon 1976; Sprott).
        % Measured here: 0.424 (2.2.x), 0.423 (2.3.0, Octave). Same tolerance as the others.
        Lhenon = lyapOf(sig.henon);
        add('henon','largest Lyapunov exponent (0.419)',0.419,Lhenon, ...
            abs(Lhenon - 0.419) < 0.05);
    catch err
        onError('logistic','NL_LargestLyap',err);
    end
end

% --- Henon: correlation dimension around 1.22 ------------------------------
if exist('NL_BoxCountEntropyRate','file') ~= 2
    add('henon','NL_BoxCountEntropyRate exists',1,0,false);   % a missing operation is a FAIL
else
    try
        r = NL_BoxCountEntropyRate(sig.henon,50,{'ac1e',5});   % signature from INP_mops_hctsa.txt
        f = fieldnames(r);
        vals = [];
        for k = 1:numel(f)
            v = r.(f{k});
            if isnumeric(v) && isscalar(v) && isfinite(v), vals(end+1) = v; end %#ok<AGROW>
        end
        % A smoke test, not a value check. NL_BoxCountEntropyRate reports statistics of
        % the Renyi entropy surface (meand<i>, mediand<i>, stdmean, ...) rather
        % than a single correlation-dimension estimate, so there is no field to
        % compare against the Henon map's known ~1.22. Recorded as such rather
        % than dressed up as a stronger check than it is.
        add('henon','NL_BoxCountEntropyRate finite output [smoke test]',1,numel(vals),numel(vals) > 0);
    catch err
        onError('henon','NL_BoxCountEntropyRate',err);
    end
end

% --- invariance: z-scored features must not change under affine rescaling --
y = sig.ar1;
try
    o1 = CO_AutoCorr(y,1,'Fourier');
    o2 = CO_AutoCorr(3*y + 7,1,'Fourier');
    add('invariance','autocorrelation unchanged by 3y+7',o1,o2,abs(o1-o2) < 1e-10);
    o1 = EN_ApEn(y,1,0.2);
    o2 = EN_ApEn(3*y + 7,1,0.2);
    add('invariance','ApEn unchanged by 3y+7',o1,o2,abs(o1-o2) < 1e-10);
catch err
    onError('invariance','affine rescaling',err);
end

% --- determinism: the same input twice must give the same answer -----------
% Guarded: CO_HistogramAMI uses histcounts2, which exists in MATLAB but not in
% Octave. A check that cannot run is reported as skipped rather than failing,
% and must not abort the rest of the suite.
try
    a1 = CO_HistogramAMI(y,1,'even',10);
    a2 = CO_HistogramAMI(y,1,'even',10);
    add('determinism','CO_HistogramAMI repeatable',a1,a2,isequaln(a1,a2));
catch err
    onError('determinism','CO_HistogramAMI repeatable',err);
end

% --- bounded quantities stay in range --------------------------------------
try
    p = DN_WithinP(sig.whitenoise,1);
    add('range','DN_WithinP is a proportion in [0,1]',NaN,p,p >= 0 && p <= 1);
catch err
    onError('range','DN_WithinP in [0,1]',err);
end

try
    r = SP_SinusoidFit(sig.sine,'sin1');
    if isstruct(r) && isfield(r,'r2')
        add('range','SP_SinusoidFit sin1 r^2 <= 1 on a sine',1,r.r2,r.r2 <= 1+1e-9 && r.r2 > 0.9);
    else
        add('range','SP_SinusoidFit returns a struct with r2',1,0,false);   % renamed field = FAIL
    end
catch err
    onError('range','SP_SinusoidFit sin1 r^2',err);
end

%-------------------------------------------------------------------------------
%% Report
%-------------------------------------------------------------------------------
% ok: nothing failed, something passed, and no more than three checks (the
% TISEAN Lyapunov block, NL_BoxCountEntropyRate and histcounts2, each one skip) were skipped.
nSkipped = sum([report.skipped]);
ok = all([report.pass]) && (numel(report) - nSkipped) > 0 && nSkipped <= 3;
if verbose
    fprintf('\n%-12s %-46s %12s %12s %6s\n','signal','check','expected','got','');
    for i = 1:numel(report)
        e = report(i).expected; g = report(i).got;
        es = fmt(e); gs = fmt(g);
        fprintf('  %-10s %-46s %12s %12s %6s\n',report(i).signal, ...
            report(i).check(1:min(46,end)),es,gs,tf(report(i).pass));
    end
    nF = sum(~[report.pass]);
    nS = sum([report.skipped]);
    fprintf('\n----------------------------------------------\n');
    if nF > 0
        fprintf(2,'  %u of %u checks FAILED.\n',nF,numel(report));
    elseif ~ok
        fprintf(2,'  NOT VALIDATED: %u checks passed, %u skipped. Too much was skipped\n', ...
            numel(report)-nS,nS);
        fprintf(2,'  (or nothing ran); build the TISEAN executables and re-run.\n');
    else
        fprintf('  All %u known-signal checks passed',numel(report)-nS);
        if nS > 0, fprintf(' (%u skipped)',nS); end
        fprintf('.\n');
    end
    fprintf('==============================================\n\n');
end

    function add(s,c,e,g,p)
        report(end+1) = struct('signal',s,'check',c,'expected',e,'got',g, ...
                               'pass',logical(p),'skipped',false); %#ok<AGROW>
    end
    function v = namedField(r,candidates)
        % Take the first NAMED field that exists and is finite. Never "whatever
        % comes first" -- field order is not a contract.
        v = NaN;
        if ~isstruct(r), return; end
        for q = 1:numel(candidates)
            if isfield(r,candidates{q})
                x = r.(candidates{q});
                if isnumeric(x) && isscalar(x) && isfinite(x), v = x; return; end
            end
        end
    end
    function onError(s,c,err)
        % Skip only for a named gap in the environment; any other error is a FAIL
        % that carries its message (HCTSA2 2.2.1: every error used to be a skip).
        if strcmp(err.identifier,'hctsa:tisean:notInstalled') || ...
                ~isempty(strfind(err.message,'histcounts2'))
            addSkip(s,c,err.message);
        else
            msg = strtrim(regexprep(err.message,'\n.*','','once'));
            add(s,[c ' [ERROR]'],NaN,NaN,false);
            if verbose
                fprintf(2,'  (ERROR in %s: %s)\n',c,msg);
            end
        end
    end
    function addSkip(s,c,msg)
        report(end+1) = struct('signal',s,'check',[c ' [SKIPPED]'],'expected',NaN, ...
                               'got',NaN,'pass',true,'skipped',true); %#ok<AGROW>
        if verbose
            fprintf('  (skipped %s: %s)\n',c,strtrim(regexprep(msg,'\n.*','','once')));
        end
    end
end

function s = fmt(v)
if isempty(v) || (isnumeric(v) && ~isscalar(v)), s = '--';
elseif ~isnumeric(v), s = '--';
elseif isnan(v), s = '--';
else, s = sprintf('%.4g',v);
end
end

function s = tf(p)
if p, s = 'ok'; else, s = 'FAIL'; end
end
