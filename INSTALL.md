# hctsa2 3.0.0 — Installation (offline Windows 11)

hctsa2 is hctsa with the **Curve Fitting, Wavelet, Econometrics, System Identification and
Financial** toolbox dependencies removed. It needs only **Statistics and Machine Learning** and
**Signal Processing**. Version 3.0.0 is merged with upstream hctsa v3.0.0 (4581e712).
**Upgrading from an earlier version:** rebuild the TISEAN programs (`HCTSA2_BuildTisean`,
section 4a: in 3.0.0 four of their sources changed, and with a `false_nearest` built before
3.0.0 every `'fnn'` feature fails) and the MEX components (`HCTSA2_BuildBinaries`, section 4b:
two sources changed), run `startup` and `HCTSA2_CheckInstall`, then recompute: the feature set
changed (`docs/HCTSA2_FeatureRenames_3.0.0.csv` maps 2.5.x names to 3.0.0's; each release's notes,
`docs/hctsa2_<version>/README.md`, list what changed in it).

Written for a machine with no internet and no administrator rights.

---

## 1. Requirements

| | |
|---|---|
| MATLAB | R2024a or later (the recorded Windows validations ran on R2025b; R2024a is the first release whose MinGW-w64 add-on includes gfortran, which four TISEAN programs need) |
| Required toolboxes | Statistics and Machine Learning, Signal Processing |
| Disk | about 20 MB unpacked, plus the programs and MEX files it builds |
| Admin rights | not needed |
| Internet | not needed |

---

## 2. Install

1. Unzip anywhere you can write — e.g. `C:\Users\<you>\Documents\hctsa2`.
   Avoid `C:\Program Files` (needs admin) and avoid paths with spaces if you can.
2. Start MATLAB, `cd` to that folder, and run:

```matlab
cd C:\Users\<you>\Documents\hctsa2
startup
```

`startup` adds every path hctsa needs. It prints what it finds, including whether the TISEAN
executables are present. **Re-run `startup` in each new MATLAB session**, or add the folder to your
MATLAB startup file to make it permanent.

---

## 3. Verify

```matlab
HCTSA2_CheckInstall
```

This is the single most useful command in the package. It reports MATLAB version, required
toolboxes, every compiled component and TISEAN, runs a live smoke test, and then
tells you **how many of the 6,392 library features are actually available**. (3.0.0: it also
runs `false_nearest` once to tell whether it is a build from before 3.0.0, which makes every
`'fnn'` feature fail.)

Run it before trusting any results. If a feature returns `NaN`, this tells you whether that is a
real result or a missing binary.

A fresh install with no compiled components reports about **5,294 of 6,392 (83%)**. That is a
perfectly usable installation — see section 5.

This fork's own checks are in `Tests/` too (3.0.0; `startup` puts the folder on the path, so they
run by name from any folder): `HCTSA2_TestRun` (the analysis pipeline, about a minute; output in
`Tests/testrun_output/`), `HCTSA2_TestTisean` (the TISEAN executables), `HCTSA2_TestSolveChol`,
`HCTSA2_ValidateKnownSignals`, `HCTSA2_TiseanProbe` and `HCTSA2_AuditTiseanCleanup`.

Upstream's unit tests run from the `Tests` folder, after the binaries are built (several tests
call TISEAN, ripser and the catch22 MEX):

```matlab
cd Tests
results = runtests('OperationsUnitTests')
results = runtests('RobustnessTests')     % a few minutes: the whole library, several series
```

`BasicPipelineTests` has 17 tests that need upstream's `Tests/HCTSA_Bonn_EEG_OG.mat` (50 MB),
which this repository leaves out: copy it from an upstream hctsa checkout into `Tests` to run
them. Without it they are skipped and the rest of the class still runs. The class works in `Tests` whatever
the current folder is, writes its working files and logs there, and deletes every other `.mat`
file there when it finishes (the fixtures and this fork's own `hctsa2_reference_input.mat`
and `hctsa2_testdata.mat` are kept).

---

## 4. Compiled components

Nothing in this package is shipped as a prebuilt executable, MEX file or Java archive.
Everything compiled is built from the bundled source with MATLAB's own MinGW-w64 add-on.
(Until 2.3.1 the tree held one compiled file, `infodynamics.jar`, upstream's Java
information-dynamics toolkit, for the Kraskov mutual-information estimators. 2.4.0 replaces
it with the KSG MEX, section 4b.)

### 4a. TISEAN — build from source

    HCTSA2_ProbeToolchain     % confirms gcc AND gfortran are present
    HCTSA2_BuildTisean        % builds 11 TISEAN programs and ripser into Toolboxes/Tisean_bin
    startup

Eleven TISEAN programs are needed -- seven C, four Fortran -- plus `ripser`
(C++) for `NL_PersistentHomology`, so **gfortran is required** for a complete
build; MathWorks ships it, and g++, with the MinGW-w64 support package
(verified on R2025b, gcc/gfortran 8.1.0). Run `HCTSA2_ProbeToolchain` first --
it compiles a test file in each of the three languages (C, Fortran, C++) rather than
trusting a version string.

| language | programs |
|---|---|
| C | `d2`, `nstat_z`, `false_nearest`, `boxcount`, `lyap_r`, `lyap_spec`, `nrlazy` |
| Fortran | `c1`, `c2d`, `c2g`, `c2t` |
| C++ | `ripser` |

The executables are linked with `-static`, so they carry no MinGW DLL
dependencies and MinGW does not need to be on `PATH` at run time. Upstream's
`install` and `compile_tisean`/`compile_ripser` call the same builder on Windows.

Without gfortran the C and C++ programs still build. The four Fortran programs
feed `NL_c1` (`c1`, `c2d`) and `NL_d2` (`d2`'s output goes through `c2g` and
`c2t`), so a C-only build loses those two operations' 120 features and nothing
else.

### 4b. MEX components

    HCTSA2_BuildBinaries      % or HCTSA2_BuildBinaries('mex_win64') to stage them

| component | features | needs |
|---|---:|---|
| KSG (Kraskov mutual information; replaced the JIDT jar in 2.4.0) | 440 | C compiler; 104 of these features also need TISEAN |
| Physionet `sampen` | 91 | C compiler |
| Michael_Small | 47 | C compiler |
| Max_Little | 36 | C and C++ compiler |
| gpml | (30) | **optional** -- `solve_chol.m` is a working MATLAB fallback, so these features run either way; the MEX is a free speed-up |
| catch22 | -- | C compiler; used by the separate catch22/catch24 feature sets, not by the main hctsa set |

`HCTSA2_BuildBinaries` compiles each component independently, never aborting the
whole run because one failed, and reports what each failure costs in features.
Passing an output directory also stages every `.mexw64` produced with its paths
preserved, for carrying to a machine without a compiler. MEX files are not
portable across major MATLAB releases, so build on a machine matching the target.

> **If the compiler seems to fail with no error message**, MinGW is not on
> `PATH`. `gcc.exe` launches `cc1.exe`, `as.exe` and `ld.exe` from its own
> installation and finds them via `PATH`; without it the driver starts but its
> subprocesses do not, and it reports nothing. `HCTSA2_BuildTisean` and
> `HCTSA2_BuildBinaries` handle this themselves.

## 5. What you get if you skip the binaries

The library computes **6,392 features from 927 master operations**. Roughly
83% are pure MATLAB and run as soon as `startup` has been called; the remainder
need the compiled components in section 4.

The pure-MATLAB baseline includes **every feature the five removed toolboxes used
to gate** — that is the entire point of hctsa2. What remains outstanding
without compilation is third-party code that stock hctsa needs too, not a
limitation introduced by this work.

Run `HCTSA2_CheckInstall` for the exact count in your installation rather than
relying on a figure in this document.

---

## 5b. Runs never abort on a missing binary

hctsa wraps every master operation in `try`/`catch` and records the outcome per feature in
`TS_Quality`. A missing MEX file or TISEAN executable produces a *caught error* for the affected
features; the run continues and everything else is computed and stored normally. hctsa2 does not
change this — its only change in `Calculation/` is a sweep of abandoned TISEAN temporary files
(older than a day) when `TS_Compute` starts.

`TS_Quality` codes: `0` computed, `1` fatal error, `2` NaN, `3` +Inf, `4` -Inf, `5` complex,
`6` empty output, `7` link error (the feature's registered output field does not exist).
Codes 2-6 are usually legitimate results for a particular series, not installation problems.
Code 7 is different: it means the feature set and the code disagree (an out-of-date `HCTSA.mat`,
or a modified operation), whatever the series.

Because failures scroll past during a long run, use:

```matlab
HCTSA2_ReportSkipped('HCTSA.mat')
```

It groups the outcome by cause afterwards and, where the failure is a missing compiled component,
names it:

```
-- fatal errors, grouped by likely cause --
  TISEAN/ripser (executables)    20 operations,     538 feature values

  Missing compiled components are the likely cause above.
  Run HCTSA2_CheckInstall to confirm, then HCTSA2_BuildBinaries.
```

---

## 6. Troubleshooting

**Features return `NaN` and you are not sure why** — run `HCTSA2_CheckInstall`. If everything
reports OK, the `NaN` is a genuine result for that time series, not an installation problem.

**MEX files copied from another machine do not load** — the MATLAB release or architecture differs.
Rebuild on a machine matching the target.

**Path with spaces causes odd failures in TISEAN calls** — the `system()` calls do not always quote
paths. Install under a path without spaces.

---

## 7. What is in this package

| path | contents |
|---|---|
| `Operations/`, `PeripheryFunctions/` | the feature library, including the 41 `BF_*` files this fork adds: the 36 toolbox replacements, `BF_mafdr` (3.0.0, for the Bioinformatics Toolbox's `mafdr`), the TISEAN runner's three and the KSG estimator's `BF_MutualInfoCalculator` (`HCTSA2_CheckProtected` checks each) |
| `Toolboxes/` | bundled third-party code, including the vendored `catch22` |
| `startup.m` | path setup, platform-aware |
| `HCTSA2_CheckInstall.m` | installation self-check |
| `HCTSA2_BuildTisean.m` | TISEAN and ripser executables, built from source |
| `HCTSA2_BuildBinaries.m` | MEX build and staging |
| `HCTSA2_ReportSkipped.m` | post-run report of what did not compute, and why |
| `HCTSA2_Analyse.m` | end-to-end pipeline: data in, three tables and plots out (`HCTSA2_Ingest.m` reads the data for it) |
| `HCTSA2_TopN.m` | narrow SeriesTable to the N best-separating features |
| `HCTSA2_ComputeFeatures.m` | chosen features of one time series, as a struct, for your own pipeline |
| `HCTSA2_SamplingSafeFeatures.m` | which features survive irregular sampling |
| `HCTSA2_CheckProtected.m` | after a merge from upstream: has anything this build depends on been reverted? |
| `HCTSA2_ProbeToolchain.m`, `HCTSA2_WhyTiseanFails.m` | is a compiler present; at which layer TISEAN stops working |
| `Tests/` | upstream's unit tests, and (3.0.0) this fork's test and debug tools: `HCTSA2_TestRun` (with `hctsa2_testdata.mat`), `HCTSA2_TestTisean`, `HCTSA2_TestSolveChol`, `HCTSA2_ValidateKnownSignals`, `HCTSA2_TiseanProbe`, `HCTSA2_AuditTiseanCleanup` |
