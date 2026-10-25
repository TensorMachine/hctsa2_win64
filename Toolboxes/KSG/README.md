# KSG — Kraskov mutual information (HCTSA2 2.4.0)

`KSG_mi.c` estimates mutual information with the Kraskov–Stögbauer–Grassberger method:
A. Kraskov, H. Stögbauer, P. Grassberger, *Estimating mutual information*, Phys. Rev. E
**69**, 066138 (2004), algorithms 1 and 2 (their equations 8 and 9).

It replaces the Java Information Dynamics Toolkit (JIDT 1.3, `infodynamics.jar`). Upstream hctsa
uses JIDT for its Kraskov estimators, and hctsa2 shipped JIDT as a precompiled jar until 2.3.1.
The new code is written from the paper, it is compiled locally like every other MEX file in the
tree, and it has no Java dependency.

## Use

    I = KSG_mi(X, Y, k, algorithm)      % nats

| input | meaning |
|---|---|
| `X` | n × dx source (one observation per row) |
| `Y` | n × 1 target |
| `k` | nearest neighbours, 1 ≤ k < n (hctsa uses 4) |
| `algorithm` | 1 or 2 |

The operations do not call it directly. `IN_Initialize_MI` returns a `BF_MutualInfoCalculator`,
a MATLAB class with JIDT's method names, so the calling code is the same as it was with JIDT.

```matlab
miCalc = IN_Initialize_MI('kraskov1','4');
miCalc.setObservations(y1, y2);
mi = miCalc.computeAverageLocalOfObservations();
```

## Conventions

These are the conventions of JIDT 1.3 as hctsa configured it: k = 4, no added noise, one thread.

- **Normalisation.** Each column is normalised to mean 0 and standard deviation 1 (n − 1). A
  constant column is only centred.
- **Distances.** All distances use the max norm: over every column in the joint space, and over
  the source columns in the source space.
- **Algorithm 1.**
  - ε is the distance to the k-th nearest neighbour.
  - n_x and n_y count the points strictly closer than ε in each marginal space.
  - I = ψ(k) − ⟨ψ(n_x+1) + ψ(n_y+1)⟩ + ψ(n).
- **Algorithm 2.**
  - ε_x and ε_y are the largest marginal distances among the k nearest neighbours.
  - n_x and n_y count the points within or on them.
  - I = ψ(k) − 1/k − ⟨ψ(n_x) + ψ(n_y)⟩ + ψ(n).
- **Errors.** Non-finite data, or n ≤ k, stops it with an error, as it stopped JIDT.
- **No noise.** Nothing random happens inside. For data with many repeated values, the callers
  add a tiny, reproducible jitter first (`BF_TieBreakNoise`).

## Agreement with JIDT

Algorithm 1 gives JIDT's values to about 1e-14 (relative): that is the order in which the
terms are summed.

Algorithm 2 agrees to the same precision, except where two neighbours sit at exactly the k-th
distance. Ties can occur only when the data repeat values. Algorithm 2 depends on which of those
neighbours are taken:

- here, the earlier point in time is taken first;
- JIDT's choice follows the layout of its k-d tree.

On series with 1–4% repeated values (too few for `BF_TieBreakNoise` to jitter), 79 of 300
algorithm-2 calls differed from JIDT, by at most 1.5e-3. The registered features that use
algorithm 2 are `CO_FirstMin`/`CO_FirstMax` `mi-kraskov2` and `IN_MultivariateAMI` `kraskov2`.
This is deviation D14 in `docs/HCTSA2_Engineering_Notes.md`.

## Method

The points are sorted by the first source column. For each point, the search runs outwards along
that order and stops on each side once that column alone is farther away than the current k-th
distance. The marginal counts are binary searches in sorted columns. A source space with more
than one column checks the other columns inside the window that the first column allows.

A grid (box-assisted) search was tried. It was faster on evenly spread data, but 1.5–5× slower
on heavy-tailed or strongly autocorrelated series, so it was not adopted.

## Build

`HCTSA2_BuildBinaries` compiles it with the other MEX files (`mex KSG_mi.c`). In Octave, use
`mkoctfile --mex KSG_mi.c`.

## Licence

GNU GPL v3, like the operations that use it (`Operations/IN_*`).
