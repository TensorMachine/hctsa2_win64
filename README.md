# hctsa2_win64 (Fewer Toolbox Dependencies / Windows build)

A fork of [hctsa](https://github.com/benfulcher/hctsa) that removes five MATLAB
toolbox dependencies and makes the whole library build and run on Windows from
source alone.

Upstream's own README is preserved as `README_upstream_hctsa.md`.

---

## A note on authorship

This fork was written collaboratively with Claude (Anthropic). The toolbox
replacements, the build tooling, the validation harnesses and this documentation
were produced in that collaboration, with every numerical claim checked against
real MATLAB output rather than assumed.

## Why this fork exists

Stock hctsa extracts ~7,700 time-series features, but a large fraction of them
are gated behind toolbox licences, and parts of it do not build on Windows. This
fork solves both problems:

**1. Five toolboxes removed.** Curve Fitting, Wavelet, Econometrics, System
Identification and Financial are no longer required. Every function hctsa used
from them was reimplemented as a `BF_*` replacement — 93 of them — matching
MathWorks' behaviour, including undocumented conventions, rather than merely
producing something plausible.

Only **Statistics and Machine Learning** and **Signal Processing** are still
needed. (Parallel Computing is optional and used only for parallel extraction.)

**2. It builds on Windows.** All compiled components — the TISEAN executables and
every MEX file — are built from bundled source with the free compiler that ships
with MATLAB (MinGW-w64). Nothing is distributed as a prebuilt binary, which matters 
if your environment does not allow importing executables.

**3. Merged with upstream v2.0.0**, so it also carries upstream's TSTOOL removal,
fixes and cleanups.

The library is currently **1,080 master operations producing 7,693 features**.
Ask `HCTSA2_CheckInstall.m` for the number in your installation rather than
trusting any figure written down.

---

## Requirements

| | |
|---|---|
| MATLAB | R2018b or later (developed and validated on R2024a and R2025b) |
| Required toolboxes | Statistics and Machine Learning, Signal Processing |
| Compiler | **MinGW-w64**, from the MATLAB Add-On Explorer |
| Disk | ~350 MB unpacked |
| Admin rights | not needed |
| Internet | not needed after the compiler is installed |

### The compiler, specifically

Install **"MATLAB Support for MinGW-w64 C/C++ Compiler"** from the Add-On
Explorer. Despite the name, MathWorks' package also includes **gfortran**, and
you need it: four of the ten TISEAN programs are Fortran. Verified on R2025b with
gcc/gfortran 8.1.0.

MSVC is not recommended. The bundled sources are patched for MinGW (see
*Third-party changes* below), and MSVC has not been tested.

---

## Getting started — run these in order

```matlab
cd <hctsa2 root>

startup                   % 1. set paths (run once per MATLAB session)
HCTSA2_ProbeToolchain     % 2. is a usable compiler present?
HCTSA2_BuildTisean        % 3. build 10 TISEAN executables from source
HCTSA2_BuildBinaries      % 4. build the MEX components
startup                   % 5. discover built executables
HCTSA2_CheckInstall       % 6. what is available
```

Order matters: `startup` reports whether the TISEAN executables exist, so build
them first. Re-run `startup` in every new MATLAB session, or add the folder to
your MATLAB startup file.

**Step 1** compiles a test file in both C and Fortran rather than trusting a
version string. If it reports a failure with an *empty* error message, MinGW is
not on `PATH` — `gcc.exe` launches `cc1.exe`, `as.exe` and `ld.exe` from its own
installation and finds them via `PATH`, so without it the driver starts and its
subprocesses do not. The build scripts handle this themselves; the probe reports
it.

**Step 5** is the one to trust. It reports how many features are actually
available and which components are missing. If a feature later returns `NaN`,
this tells you whether that is a real result or an incomplete install.

### Then, to analyse data

```matlab
results = HCTSA2_Analyse('mydata.mat');
```

Takes a `.mat` containing your time series and returns two tables:
`FeatureTable` (one row per feature: separability, rank, per-class means) and
`SeriesTable` (one row per series: group plus all feature values). It accepts
long-format tables (`id, time, value, class`), wide tables, numeric matrices and
cell arrays.

`HCTSA2_TopN(results, 20)` narrows `SeriesTable` to the 20 best-separating
features.

---

## What gets built

### TISEAN — 10 standalone executables

Built from the bundled TISEAN 3.0.1 source into `Toolboxes/Tisean_bin`. hctsa
calls these as command-line programs, not MEX.

| language | programs |
|---|---|
| C | `d2`, `nstat_z`, `false_nearest`, `boxcount`, `lyap_r`, `poincare` |
| Fortran | `c1`, `c2d`, `c2g`, `c2t` |

Linked with `-static`, so they carry no MinGW DLL dependencies and MinGW does not
need to be on `PATH` at run time.

`HCTSA2_BuildTisean` replaces TISEAN's `./configure && make`, which needs a POSIX
shell that Windows does not have.

### MEX components

| component | features | needs |
|---|---:|---|
| catch22 | 240 | C |
| Physionet `sampen` | 49 | C |
| Michael_Small | 33 | C |
| Max_Little | 10 | C and C++ |
| gpml | (72) | optional — a MATLAB fallback works, the MEX is a speed-up |

`HCTSA2_BuildBinaries` compiles each independently, never aborting the whole run
because one failed, and reports what each failure costs in features.

---

## Third-party libraries and the changes made to them

Everything below is bundled in `Toolboxes/`. The changes were needed either to
compile under MinGW on 64-bit Windows or to fix genuine defects found while
validating.

| library | what was changed | why |
|---|---|---|
| **TISEAN 3.0.1** | `false_nearest.c`: applied the authors' official 2009-03-08 patch (`vemb[i] = i*delay`), plus a stricter unsigned-underflow guard | the univariate branch omitted the delay from its embedding, so every call with tau > 1 embedded at the wrong lags |
| **gpml 4.2** | `solve_chol.c`: added MinGW to the LAPACK symbol test | MSVC exports `dpotrs`, MinGW needs `dpotrs_`; without this the build calls an undeclared symbol |
| | `make.m`: corrected the import-library path | it built `.../win64/gnu/` from the compiler manufacturer string, but MathWorks ships MinGW libraries under `.../win64/mingw64/` |
| **Michael_Small** | `MS_shannon.c`, `MS_complexitybs.c`: fixed the `qsort` comparator | it cast 8-byte doubles to 4-byte `int`, so the array was not sorted at all — the percentile thresholds built on it were meaningless, and the result varied by compiler |
| **catch22** | vendored at commit `9ff9da73` instead of a git submodule; upstream test fixtures and images removed | a submodule cannot be fetched on an offline machine |

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

**4. Stochastic features tested statistically.** Operations that draw random
numbers cannot match exactly, so they are run repeatedly and the reference is
scored against the resulting distribution. 95% fall within 2σ, 98% within 3σ.

Three genuine bugs in hctsa itself were found along the way and fixed: the
`false_nearest` delay embedding, a Pearson median skewness formula missing a pair
of parentheses, and an off-by-`d` index in `CO_TranslateShape`. Where a fix makes
a feature differ from the reference, the reference is the one that was wrong.

---