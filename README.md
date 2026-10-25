# hctsa3_win64 (Fewer Toolbox Dependencies / Windows build)

A fork of [hctsa](https://github.com/benfulcher/hctsa) that removes five MATLAB
toolbox dependencies and makes the whole library build and run on Windows from
source alone. No pre-compiled binaries.

**Version 3.0.0** merges upstream hctsa's **v3.0.0 release** (`4581e712`, 2026-10-20)
Note: Many file names still retain the HCTSA2 naming to differentiate from upstream files

Upstream's own README is preserved as `README_upstream_hctsa.md`. Upstream notes there that
most new development of the feature library is moving to
[_pyhctsa_](https://github.com/DynamicsAndNeuralSystems/pyhctsa), a native Python
implementation; _hctsa_ itself remains maintained.

---

## A note on authorship

This fork was written collaboratively with Claude (Anthropic). The toolbox
replacements, the build tooling, the validation harnesses and this documentation
were produced in that collaboration, with every numerical claim checked against
real MATLAB output rather than assumed.

## Why this fork exists

Stock hctsa extracts thousands of time-series features (6,392 in v3.0.0), but a large fraction of them
are gated behind toolbox licences, and parts of it do not build on Windows. This
fork solves both problems:

**1. Five toolboxes removed.** Curve Fitting, Wavelet, Econometrics, System
Identification and Financial are no longer required. Every function hctsa used
from them was reimplemented in `BF_*` replacements (36 files in `PeripheryFunctions/`; the
fork adds 41 there in all, with the TISEAN runner's three, `BF_mafdr` for the one Bioinformatics
Toolbox call in upstream's analysis code, and the KSG estimator's `BF_MutualInfoCalculator`),
matching MathWorks' behaviour

Only **Statistics and Machine Learning** and **Signal Processing** are still
needed. (Parallel Computing is optional and used only for parallel extraction.)

**2. It builds on Windows.** All compiled components — the TISEAN executables and
every MEX file — are built from bundled source with the free compiler that ships
with MATLAB (MinGW-w64). Nothing is distributed as a prebuilt executable, DLL or MEX
file, which matters if your environment does not allow importing executables. The Java 
information-dynamics toolkit (`infodynamics.jar`), which upstream uses for its 
mutual-information estimators, is replaced by C code compiled like the rest.

The library is currently **927 master operations producing 6,392 features**.
Ask `HCTSA2_CheckInstall.m` for the number in your installation rather than
trusting any figure written down.

---

## Requirements

| | |
|---|---|
| MATLAB | R2024a or later |
| Required toolboxes | Statistics and Machine Learning, Signal Processing |
| Compiler | **MinGW-w64**, from the MATLAB Add-On Explorer |
| Disk | about 20 MB unpacked, plus the programs and MEX files it builds |
| Admin rights | not needed |
| Internet | not needed after the compiler is installed |
| Java | Java is required for various functionality

### The compiler, specifically

Install **"MATLAB Support for MinGW-w64 C/C++ Compiler"** from the Add-On
Explorer. Despite the name, MathWorks' package also includes **gfortran**, and
you need it: four of the eleven TISEAN programs are Fortran (and ripser is C++,
also included).

MSVC is not recommended. The bundled sources are patched for MinGW (see
*Third-party changes* below), and MSVC has not been tested.

**R2024a is the floor** because that is the release in which MathWorks added
gfortran to the supported compilers on Windows, and four of the eleven TISEAN
programs are Fortran. On an earlier release the C and C++ programs still build.

---

## Getting started — run these in order

```matlab
cd <hctsa2 root>

startup                   % 1. set paths (run once per MATLAB session)
HCTSA2_ProbeToolchain     % 2. is a usable compiler present?
HCTSA2_BuildTisean        % 3. build the TISEAN programs and ripser (.exe)
HCTSA2_BuildBinaries      % 4. build the MEX components
startup                   % 5. discover what was built
HCTSA2_CheckInstall       % 6. what is available
```

**Upgrading from an earlier version:** run steps 3 to 6 again. In 3.0.0 four TISEAN sources
changed (`false_nearest`, which now writes a fifth column that `NL_FNN` requires, `lyap_r`, `d2`
and `c1`'s `d1`) and two MEX sources (`sampen_mex.c`, `MS_complexitybs.c`). With an old
`false_nearest` every `'fnn'` feature fails; `HCTSA2_CheckInstall` and `HCTSA2_TestTisean` report
it. (Each release's notes, `docs/hctsa2_<version>/README.md`, list what changed in it.)

Order matters: `startup` reports whether the TISEAN executables exist, so build
them first. Re-run `startup` in every new MATLAB session, or add the folder to
your MATLAB startup file. (Upstream's `install`, `compile_tisean` and
`compile_ripser` also work: on Windows they call the same builders.)

**Step 2** compiles a test file in each of C, Fortran and C++ rather than trusting a
version string. If it reports a failure with an *empty* error message, MinGW is
not on `PATH` — `gcc.exe` launches `cc1.exe`, `as.exe` and `ld.exe` from its own
installation and finds them via `PATH`, so without it the driver starts and its
subprocesses do not. The build scripts handle this themselves; the probe reports
it.

**Step 6** is the one to trust. It reports how many features are actually
available and which components are missing. If a feature later returns `NaN`,
this tells you whether that is a real result or an incomplete install.

### Then, to analyse data

```matlab
results = HCTSA2_Analyse('mydata.mat','Groups',g,'OutputDir','run1');
results = HCTSA2_Analyse(X,'Groups',g);                    % a matrix, cell array or table
results = HCTSA2_Analyse([],'ExistingHCTSA','run1\HCTSA.mat','Groups',g2,'OutputDir','run2')
```

Runs the whole pipeline unattended (extraction, normalisation, ranking,
classification, figures) and returns three tables:

- `FeatureTable`: one row per feature of the run, kept or not. `Status` says why a
  feature was dropped. Then `PctComputed`, `Rank` and `Separability`; for two classes
  `AUC` and `HigherIn`; per group the normalised mean and the raw median and IQR;
  `RedundantWith` for top features that duplicate a better-ranked one, and `Distinct`
  (yes, no, or not checked: the check covers the top 200, `'RedundancyTopN'`); with
  `'ReduceFeatureSet'`, `Cluster`, `ClusterSize` and `Representative` (below); and
  `Column`, the feature's column in `SeriesTable`.
- `SeriesTable`: one row per time series: its group and every kept feature's
  normalised value. Long feature names are shortened to MATLAB's 63-character limit
  without colliding.
- `FeatureAUC` (also `FeatureAUC.csv`): one row per kept feature, in rank order: how well
  it separates the classes. `MeanAUC` and `MinAUC` over every pair of classes; per pair
  `AUC_A_vs_B`, the class with the higher values and a Mann-Whitney q-value; with three or
  more classes, each class against the rest; a Kruskal-Wallis p- and q-value across all
  classes; the % of series the feature computed on, and on how many of each class; and
  the redundancy columns. An AUC of 0.5 is no separation, 1 is perfect; with few series per class it sits
  above 0.5 by chance (about 0.68 at 4 per class, 0.57 at 20), so read it with the q-values
  (Benjamini-Hochberg over the features).
- `DroppedSeries`: the series left out, and why.

`'ReduceFeatureSet',true` also reduces the whole feature set: every feature kept by
normalisation is clustered by complete linkage on 1 − |Spearman ρ|, so that every pair within a
cluster has |ρ| ≥ `'RedundancyThreshold'` (0.8). Each cluster keeps its best-ranked member (with
fewer than two classes, its most central one). `HCTSA_N_reduced.mat` is `HCTSA_N.mat` with only
the kept features, readable by every `TS_*` function, plus a `reducedFeatureSet` table giving
every feature's cluster and the feature kept for it; `FeatureTable` and `FeatureAUC` gain
`Cluster`, `ClusterSize` and `Representative` ("" for the kept feature itself).

It accepts `.mat` and `.csv` files and data in memory: numeric matrices, cell arrays
of vectors, and wide or long-format tables (`id, time, value, class`). NaN padding at
either end of a series is trimmed. A series with a NaN inside it is dropped and
listed. It never prompts, and it never modifies its input: an `'ExistingHCTSA'` file
is copied into `OutputDir` first, so a computed dataset can be re-analysed under new
groups without recomputing. Everything goes to `OutputDir`, including a log of the
console output. Running again into the same `OutputDir` with the same data resumes
the extraction where it stopped.

The main options (the function's help lists them all):

| option | default | |
|---|---|---|
| `'Features'`, `'ExcludeFeatures'` | all | names, wildcards allowed (`'ExcludeFeatures','NL_DVV*'`), applied before extraction |
| `'Parallel'`, `'ParallelMode'` | false, `'series'` | `'series'` for many series, `'operations'` for a few very long ones |
| `'FilterThresholds'` | `[0.7 0.7]` | the fraction of good values a series and a feature need. Features that fail on up to 30% of series are kept, so they can show whether they separate the classes (`[0.7 1]` is hctsa's own default) |
| `'NaNToZero'` | false | set the NaNs left after normalisation to 0 in `HCTSA_N.mat`. 0 is the bottom of the sigmoid scale (0.5 is the median), so "not computed" reads as "lowest value" |
| `'ClassVarFilter'` | false | also drop features that are constant within a class |
| `'TestStatistic'` | `'ustatP'` (2 classes), `'meanAUC'` (3+) | ranking statistic: the Mann-Whitney p-value; for more classes the mean of the pairwise AUCs, ties broken by the lowest. `'classification'` (2.3.0's default for 3+) remains available |
| `'Classifier'` | `'tree'` | any classifier `TS_Classify` accepts. Most classifiers other than the tree leave out every series with a NaN feature; the tree uses them, and so does every classifier after `'NaNToZero'` |
| `'NumFolds'`, `'NumRepeats'`, `'NumNulls'` | hctsa's `HowManyFolds`, 1, 0 | cross-validation folds (2, 5 or 10 by the size of the smallest class), repeats, permutation nulls |
| `'DoPlots'`, `'ShowFigures'`, `'CloseFigures'`, `'FigureFormat'` | true, true, true, `'png'` | figures are saved in `OutputDir` |

`HCTSA2_TopN(results, 20)` narrows `SeriesTable` to the 20 best-separating features.

`HCTSA2_ComputeFeatures(x, features)` computes chosen features of one time series, for
use inside your own feature-extraction pipeline: give it the names (or IDs) picked from
`FeatureAUC`, and it returns a struct with one field per feature, named after it. Each master
operation is computed once and only the requested features kept; the values are the ones
`HCTSA2_Analyse` stores in `HCTSA.mat` for the same series (raw by default; NaN where a
feature fails). With `'NormalizeLike'` (3.0.0) they are normalised exactly as a run normalised
its data: since v3.0.0 `TS_Normalize` stores each feature's fitted normalisation in `HCTSA_N.mat`,
and it is applied to the new series with the run's parameters, so a series of the run gets its
`HCTSA_N.mat` values. Store names rather than IDs: an ID is a row of the library's list and moves
when the library changes (only 2% of IDs kept their feature across the 3.0.0 merge).

```matlab
sel = R.FeatureAUC.Name(1:20);              % chosen from an HCTSA2_Analyse run
F = HCTSA2_ComputeFeatures(x, sel);         % F.mean, F.SB_MotifTwo_mean_hhh, ...
v = cell2mat(struct2cell(F))';              % a row vector, in the order of sel
Fn = HCTSA2_ComputeFeatures(x, sel, 'NormalizeLike', 'run1\HCTSA_N.mat');   % normalised
```

The accuracy it reports is a rough indication of how well the features separate the
classes, not a tuned model. A decision tree needs at least 10 series in a node to
split, so on a very small dataset (a handful of series per class) it reports chance;
choose `'Classifier','svm-linear-lowdim'` there. `HCTSA2_TestRun` runs nine small
scenarios through `HCTSA2_Analyse` and checks the results, then `HCTSA2_ComputeFeatures` and
`TS_TopFeatures`' p- and q-values with permutation nulls (about a minute). It lives in `Tests/` since 3.0.0, with the other test and debug tools
(`HCTSA2_TestTisean`, `HCTSA2_TestSolveChol`, `HCTSA2_TiseanProbe`, `HCTSA2_AuditTiseanCleanup`,
`HCTSA2_ValidateKnownSignals`); `startup` puts `Tests/` on the path, so they run by name from any
folder, and `HCTSA2_TestRun` writes its output to `Tests/testrun_output/`.

---

## What gets built

### TISEAN — 11 executables, and ripser

Built from the bundled TISEAN 3.0.1 source (and `Toolboxes/ripser`) by
`HCTSA2_BuildTisean`, which replaces TISEAN's `./configure && make` — that needs
a POSIX shell Windows does not have. Output goes to `Toolboxes/Tisean_bin/`. (On
Linux and macOS, upstream's `compile_tisean`/`compile_ripser` build into
`Toolboxes/Tisean_3.0.1/bin` and `Toolboxes/ripser/bin`, and `startup` finds
them there or in `~/bin`.)

| language | programs |
|---|---|
| C | `d2`, `nstat_z`, `false_nearest`, `boxcount`, `lyap_r`, `lyap_spec`, `nrlazy` |
| Fortran | `c1`, `c2d`, `c2g`, `c2t` |
| C++ | `ripser` (persistent homology; not TISEAN, run the same way) |

They account for **538 feature values across 20 master operations**, plus 50
more that reach `false_nearest` indirectly through an `'fnn'` embedding dimension.

Every operation calls them through `BF_RunTisean` rather than `system()`, so
the timeout, the not-installed check and the temporary-file cleanup below apply
uniformly. The executables are linked with `-static`, so they carry no MinGW
DLL dependencies.

An earlier version of this project also built the TISEAN programs as MEX files, to
avoid shipping executables to a machine with a strict security posture. That
was removed. Running in-process means there is no child process to kill, and
therefore no way to stop a `d2` call that will not finish for ten hours — and
as the next section shows, that is not a hypothetical. A recoverable error beat
an unkillable one.

### The TISEAN timeout

Two TISEAN programs can run for hours on long series. Neither loops forever —
every loop in both is bounded — but both grow with roughly the *square* of the
series length, and `d2`'s early exit is switched off by the flags hctsa passes
it:

| program | why it grows | which operations |
|---|---|---|
| `d2` | `-N0` sets the pair limit to `ULONG_MAX` (`d2.c:121`), making the only early exit (`d2.c:579`) unreachable | `NL_Dimensions`, `NL_GPCorrSum` (via `Nref = -1`) |
| `false_nearest` | no pair limit exists; points farther than `varianz/rt` from any neighbour are never satisfied (`false_nearest.c:196`) and exact duplicates (since 3.0.0, anything within round-off) are rejected (`:185`), so flat stretches are rescanned at every epsilon step | `NL_FNN`, and anything using `BF_Embed(...,'fnn')` |

Measured at `-M1,10`: `d2 -N0` takes 2.3 s at N=4000, 9.5 s at 8000 and 38.4 s
at 16000 — four times the work for twice the data, which puts 500k points at
about ten hours. `false_nearest` takes 4.1 s on 16000 clean points and 41.8 s
when half the series is constant.

`BF_RunTisean` therefore gives each executable a wall-clock budget, 300 seconds
by default, and raises an error if it is exceeded. hctsa records that as a NaN
feature, which is the intended behaviour: an error is recoverable, a hung
session is not. Change it with

```matlab
setenv('HCTSA2_TISEAN_TIMEOUT','1800')      % seconds; 0 disables the limit
pctRunOnAll setenv('HCTSA2_TISEAN_TIMEOUT','1800')   % on a parallel pool
```

This is a trade-off, not a free win. A `d2 -N0` call on 100k points needs about
25 minutes to finish legitimately and the default kills it at five, so those
features become NaN. Raise the budget if you want them.

The timeout works because there is a child process to kill. This is the main
reason the project ships executables rather than MEX files.

### Temporary files

Every TISEAN call writes a temporary time-series file named `hc######` in the
system temporary directory, and the programs write their results next to it
(`.c2`, `.d2`, `.h2`, `.stat`, `.box`, `.ros`, `.c1`, `.fnn`). In the
original hctsa, seven of the ten TISEAN operations of the time never deleted that
input file, so a full run left one behind per operation per time series, and an
aborted call left five.

Names are drawn from a private random stream seeded per process and claimed
atomically (`java.io.File.createNewFile`), so parallel workers — or two MATLAB
sessions sharing one temporary directory — can never be handed the same name.

`BF_CleanTempFiles` removes a whole set given any member of it. Each operation
calls it on exit, and `BF_RunTisean` calls it on any path that raises — a
timeout, a missing backend — because that abort takes the operation down before
its own cleanup can run. Cleanup does **not** happen on a normal return from
`BF_RunTisean`: the caller usually still needs what TISEAN just wrote, since
`NL_d2` runs `d2` and then feeds the resulting `.c2` to `c2g` and `c2t`.

`BF_CleanTempFiles` only ever matches `hc` followed by **exactly six digits**,
optionally then `.` or `_` and a suffix. The obvious implementation, deleting
`tempdir/hc*`, would be unsafe: a parallel pool shares one temporary directory,
so a wildcard sweep would delete the input files of other workers mid-run, and
the resulting failures — a program reporting a missing file it was handed
moments ago, on a random subset of series — would be very hard to trace back.

For what in-process cleanup cannot catch — MATLAB killed outright, a crash, a
power cut, or files from an older hctsa — `TS_Compute` sweeps the temporary
directory once before it starts (`BF_SweepStaleTempFiles`): a file set is removed
only when its newest file is more than 24 hours old, so a run in progress, in
this MATLAB or any other, is never touched. Call `BF_SweepStaleTempFiles(hours)`
yourself to use another threshold.

Each operation guards its temporary file with a `try`/`catch` around everything
after `BF_WriteTempFile`, plus a cleanup call before every `return`. Both are
needed, and it took three attempts to get right — each one missed a whole
*class* of exit:

| attempt | mechanism | what it missed |
|---|---|---|
| 1 | `onCleanup` alone | does not fire in a file containing a **nested** function, because the parent workspace outlives the return — six of the ten operations of the time had one |
| 2 | explicit calls before `return` and at the end of the body | all **43** `error()` calls |
| 3 | `try`/`catch` around the body | `try`/`catch` does not intercept **`return`**, and nine were written `out = NaN; return` on one line, so a regex matching `return` alone on a line found none of them |

`try`/`catch` covers `error()` and anything a called function throws; the
per-return calls cover what `try`/`catch` cannot see. There is no third way out
of a MATLAB function, so together they are exhaustive.

Because the rule is easy to break and the breakage is invisible until a
temporary directory fills up, it is checked mechanically rather than by review.
`HCTSA2_AuditTiseanCleanup` reads the source of every operation that calls
`BF_WriteTempFile` — discovering them rather than working from a list, so one
added later is audited without anyone remembering — and proves the file is
removed on every path. `HCTSA2_TestTisean` runs it as one of its checks.

### Five crashes and three hangs, fixed in the TISEAN sources

The programs were swept over hctsa's own command shapes across series lengths,
delays, embedding dimensions and Theiler windows.
**145 of 756 runs terminated abnormally.** Every one was an unsigned subtraction
that wrapped when an embedding window exceeded the series length, or an epsilon
search with no exit when the requested neighbours could never be found:

| program | failure | fix |
|---|---|---|
| `lyap_r` | segfault (`put_in_boxes`, `length-del-steps` underflow) | signed bound plus a length check |
| `lyap_r` | infinite epsilon growth | stop at the data diameter |
| `nstat_z` | segfault (`clength` underflow); its own guard was written in unsigned arithmetic and passed exactly when it should have fired | length check, and the guard evaluated signed |
| `nstat_z` | infinite epsilon growth | stop at the data diameter |
| `boxcount` | segfault (`LENGTH-(maxembed-1)*DELAY` underflow) | length check |
| `poincare` | segfault (`length-(dim-comp)*delay-1` underflow) | length check (no longer built since 2.5.0) |
| `c1` | segfault and hang (`log` of a non-positive count when `-t` reaches the series length) | length check, and a bound on `d1`'s sweep (since 2.5.0 upstream's own `c1`/`d1` fixes, which skip what cannot be computed) |
| `c1` (2.5.0) | segfault on a constant series: zero spread made `d1`'s box size zero, and the neighbour search divided by it (stock upstream reads out of bounds there too) | the mass is skipped, as for too few points |
| `lyap_spec` (new in 3.0.0-alpha) | endless loop: when every neighbourhood so far had size 0 (exact duplicates), the search radius started at 0 and stayed there | start from the default radius in that case |

The `poincare` one is worth a note on method: the first sweep missed it because
it varied `-q` together with `-m`, the way hctsa calls it. A later randomised
fuzz that varied every option independently found it immediately. A grid built
from how the caller behaves tests the caller, not the program.

The current position: **0 abnormal terminations in 756 grid runs, in 3,000
randomised runs across lengths from 3 to 2,000 points, and in 45 malformed-input
runs** through the three Fortran post-processors; and, for the two programs the
3.0.0-alpha merge added, 0 in 1,500 randomised runs each of `lyap_spec` and
`nrlazy`. 2.2.0's `lyap_spec` Theiler option (`-t`) and the multi-column input
`NL_LyapSpec` now passes were fuzzed over 1,500 more runs: no crash, and one run
past the fuzzer's 60 s limit, a spike train with a starting radius (`-r 4.8e-5`)
hctsa never passes, which finishes in about 140 s with the same result on the
2.1.1 build. 3.0.0's four changed programs (`false_nearest`, `lyap_r`, `d2`, `c1`) were fuzzed
over 1,500 runs (0 abnormal), and, built with bounds checking (`-fsanitize=address,undefined`,
`-fcheck=all`), over 600 more: that build found the two out-of-bounds reads in `d2` listed under
*Third-party changes*, and after the guards none. These are fixes to the TISEAN
sources, so the executable build gets them too. Every guard fires only on an
input that previously crashed or hung — output on all working inputs is
byte-for-byte unchanged, which is what `HCTSA2_TestTisean` establishes.

### MEX components

| component | features | needs |
|---|---:|---|
| Physionet `sampen` | 91 | C |
| Michael_Small | 47 | C |
| KSG (Kraskov mutual information) | 440 | C; 104 of these features also need TISEAN |
| Max_Little | 36 | C and C++ (`rpde`, `steps_bumps`) |
| gpml | (30) | optional — a MATLAB fallback works, the MEX is a speed-up |
| catch22 | — | C; used by the separate catch22/catch24 feature sets, not the main one |

`HCTSA2_BuildBinaries` compiles each independently, never aborting the whole run
because one failed, and reports what each failure costs in features.

---

## Third-party libraries and the changes made to them

Everything below is bundled in `Toolboxes/`. The changes were needed either to
compile under MinGW on 64-bit Windows or to fix genuine defects found while
validating.

| library | what was changed | why |
|---|---|---|
| **TISEAN 3.0.1** | `false_nearest.c`: applied the authors' official 2009-03-08 patch (`vemb[i] = i*delay`), plus a stricter unsigned-underflow guard | the univariate branch omitted the delay from its embedding, so every call with tau > 1 embedded at the wrong lags. Upstream hctsa adopted the fix in `f54e4d1` (2.5.0) |
| | `lyap_r.c`, `nstat_z.c`, `boxcount.c`, `poincare.c`, `c1.f`, `d1.f`, `lyap_spec.c`: length guards and bounded epsilon searches (see above) | five segfaults and three hangs on series hctsa can hand them; `nstat_z`'s own guard was written in unsigned arithmetic and never fired |
| | `d2.c` (pair scramble) and `routines/rand.c`, `tsa.h`, with their callers `lyap_spec.c`, `makenoise.c`, `recurr.c`, `rand_arb_dist.c`: the random-number generators use `unsigned long long` | both chose their generator with `sizeof(long) == 8`, and `long` is 32-bit on Windows, so Windows builds ran a different generator: `NL_d2`'s correlation sums differed from Linux/macOS by up to 36%. Now every platform runs the 64-bit generator and the Windows programs compute the same numbers as Linux (on native Windows the C programs write exponents with three digits, `e-001`, and the Fortran post-processors `c2d`/`c2g`/`c2t` differ in the last single-precision digit); Linux and macOS output is unchanged |
| | `d2.c` (3.0.0): two out-of-bounds reads guarded — the 2-D box coordinate past the end of the series with an embedding of 1, and the `.stat` file's epsilon one past its array | found by fuzzing a bounds-checked build (stock TISEAN reads them too); neither reached a feature value |
| | `c1.f`, `c2d.f`, `c2g.f`, `c2t.f`, `help.f`: filename buffers widened from `character*72` to `character*512` | gfortran truncates a longer name **silently**, so the program writes to a path nobody asked for and the caller then looks for a file that was never created. A Windows `tempname` is 72 characters on its own. Only the filename variables were widened — `aline` and the option-list buffers are sliced with a hardcoded `:72` and keep their size |
| **gpml 4.2** | `solve_chol.c`: added MinGW to the LAPACK symbol test | MSVC exports `dpotrs`, MinGW needs `dpotrs_`; without this the build calls an undeclared symbol |
| | `make.m`: corrected the import-library path | it built `.../win64/gnu/` from the compiler manufacturer string, but MathWorks ships MinGW libraries under `.../win64/mingw64/` |
| **Max_Little** | removed `fastdfa/ML_fastdfa_core.dll` | `HCTSA2_BuildBinaries` compiles from their `.c` sources |
| **Michael_Small** | `MS_shannon.c`, `MS_complexitybs.c`: fixed the `qsort` comparator | it cast 8-byte doubles to 4-byte `int`, so the array was not sorted at all — the percentile thresholds built on it were meaningless, and the result varied by compiler |
| **ripser** (added upstream in 3.0.0-alpha) | `ripser.cpp`: `index_t(1)` instead of `1l` in one shift | `long` is 32-bit on Windows (LLP64), so the shift overflowed; the value is only read with `USE_COEFFICIENTS`, which is not set, so output is unchanged |
| **catch22** | vendored at commit `fc45e608` (the one upstream pins) instead of a git submodule, and not fetched by `startup`; upstream test fixtures and images removed | a submodule cannot be fetched on an offline machine |

---

## How the toolbox removal was validated

Reimplementing a toolbox function is easy to do *plausibly* and hard to do
*correctly*. The validation was built around that distinction.

**1. Against real MATLAB output.** The reference is an export of 50 time series ×
7,702 features computed on a machine with all five toolboxes licensed. Every
`BF_*` replacement was checked against it feature by feature, and every
disagreement had to be explained — not merely tolerated. That process resolved
conventions that are not documented anywhere, such as `armax`'s default initial
conditions and `n4sid`'s horizon selection.

For the 3.0.0-alpha merge, upstream's own code was run on licensed MATLAB R2025b
over 20 reference series, with probes of the toolbox functions it uses. The new
wavelet and Morlet features match on every series to within 2.5e-9; the bounded
curve fits match to their stopping tolerance; the GARCH fits start from exactly
MATLAB's starting points and agree in 59 of 72 cases, reaching a higher
likelihood in the rest, where MATLAB's optimiser stopped early

**2. Every difference attributed.** The final state is **zero unexplained
mismatches**. Differences that remain are each traced to a named cause:
documented deviations where our optimiser converges and MathWorks' stops early,
stochastic operations, chaotically sensitive operations, and operations upstream
redefined after the reference was made.

**3. Against analytically known answers.** Reference comparison cannot catch a
feature that returns a plausible number for the wrong reason, so the library is
also run on signals whose properties are known in closed form.
`HCTSA2_ValidateKnownSignals` checks, among others, that an AR(1) process with
φ = 0.8 gives autocorrelations of 0.8ᵏ, and that the logistic map's largest
Lyapunov exponent comes back as **ln 2 = 0.693** — measured 0.692, within 0.17%.
That single number exercises the TISEAN binary, the embedding, the scaling-range
selection and the replacement regression all at once.

Three genuine bugs in hctsa itself were found along the way and fixed: the
`false_nearest` delay embedding, a Pearson median skewness formula missing a pair
of parentheses, and an off-by-`d` index in `CO_TranslateShape`. Where a fix makes
a feature differ from the reference, the reference is the one that was wrong.

---
