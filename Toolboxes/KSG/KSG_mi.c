/* KSG_mi.c -- Kraskov-Stoegbauer-Grassberger estimate of mutual information (HCTSA2)
 *
 *   I = KSG_mi(X, Y, k, algorithm)
 *
 *   X          n x dx source (dx >= 1)
 *   Y          n x 1 target
 *   k          number of nearest neighbours, 1 <= k < n
 *   algorithm  1 or 2: equations 8 and 9 of Kraskov, Stoegbauer & Grassberger,
 *              "Estimating mutual information", Phys. Rev. E 69, 066138 (2004)
 *   I          mutual information in nats
 *
 * The conventions are those of JIDT 1.3 (MutualInfoCalculatorMultiVariateKraskov1/2),
 * which computed these features before HCTSA2 2.4.0, as hctsa configured it: k = 4,
 * no added noise, one thread. Written from the paper and tested against JIDT's
 * output; no JIDT code.
 *   - Each column is normalised to mean 0 and standard deviation 1 (n-1); a constant
 *     column is only centred.
 *   - Distances are max norms: over all columns in the joint space, over the source
 *     columns in the source space.
 *   - Algorithm 1: eps = distance to the k-th nearest neighbour; n_x, n_y = points
 *     strictly closer than eps in each marginal space;
 *     I = psi(k) - mean(psi(n_x + 1) + psi(n_y + 1)) + psi(n).
 *   - Algorithm 2: eps_x, eps_y = largest marginal distances among the k nearest
 *     neighbours; n_x, n_y = points within or on them;
 *     I = psi(k) - 1/k - mean(psi(n_x) + psi(n_y)) + psi(n).
 *   - Non-finite data, or n <= k, is an error (JIDT stops too).
 *   - Neighbours at exactly the k-th distance (only with repeated values) are taken
 *     in time order, the earlier first. Algorithm 2 depends on which are taken;
 *     JIDT's choice follows its search tree, so there the two can differ.
 *
 * Neighbour search: the points are sorted by the first source column and searched
 * outwards from each point, stopping on each side once that column alone is farther
 * than the current k-th distance. Marginal counts: binary search in sorted columns
 * (the source space, when it has more than one column, checks the others inside the
 * window the first column allows).
 *
 * Part of HCTSA2; GNU GPL v3, like the operations that use it.
 */
#include "mex.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct { double v; mwSize i; } Pair;

static int comparePairs(const void *a, const void *b)
{
    const Pair *p = (const Pair *)a, *q = (const Pair *)b;
    if (p->v < q->v) return -1;
    if (p->v > q->v) return 1;
    return (p->i > q->i) - (p->i < q->i);           /* ties: original order */
}

static int compareDoubles(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

/* The range [*first, *last) of sorted v[0..n) with |v - c| < r (strict) or <= r. */
static void window(const double *v, mwSize n, double c, double r, int strict,
                   mwSize *first, mwSize *last)
{
    mwSize lo = 0, hi = n, m;
    while (lo < hi) {                                /* first value close enough */
        m = lo + (hi - lo) / 2;
        double d = c - v[m];
        if (d > r || (strict && d == r)) lo = m + 1; else hi = m;
    }
    *first = lo;
    hi = n;
    while (lo < hi) {                                /* first value beyond it */
        m = lo + (hi - lo) / 2;
        double d = v[m] - c;
        if (d < r || (!strict && d == r)) lo = m + 1; else hi = m;
    }
    *last = lo;
}

/* Insert a neighbour into the k best so far, sorted by joint distance and, at equal
 * distance, by position in time (the earlier point first). */
static void keep(double dist, mwSize at, double dx, double dy, mwSize k, mwSize *nb,
                 double *bd, mwSize *bt, double *bx, double *by)
{
    mwSize t;
    if (*nb < k) t = (*nb)++;
    else if (dist < bd[k - 1] || (dist == bd[k - 1] && at < bt[k - 1])) t = k - 1;
    else return;
    while (t > 0 && (bd[t - 1] > dist || (bd[t - 1] == dist && bt[t - 1] > at))) {
        bd[t] = bd[t - 1]; bt[t] = bt[t - 1]; bx[t] = bx[t - 1]; by[t] = by[t - 1];
        t--;
    }
    bd[t] = dist; bt[t] = at; bx[t] = dx; by[t] = dy;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSize n, dx, D, k, p, q, d, i, nb, a, b;
    int alg;
    double *X, *Y, kd, ad;

    (void)nlhs;
    if (nrhs != 4)
        mexErrMsgIdAndTxt("KSG_mi:nargin", "Usage: I = KSG_mi(X, Y, k, algorithm)");
    for (i = 0; i < 4; i++)
        if (!mxIsDouble(prhs[i]) || mxIsComplex(prhs[i]) || mxIsSparse(prhs[i]))
            mexErrMsgIdAndTxt("KSG_mi:type", "All inputs must be real, full double arrays.");
    n = mxGetM(prhs[0]);
    dx = mxGetN(prhs[0]);
    if (n < 1 || dx < 1 || (mwSize)mxGetNumberOfElements(prhs[1]) != n)
        mexErrMsgIdAndTxt("KSG_mi:size", "X must be n x dx and Y must have n elements.");
    if (mxGetNumberOfElements(prhs[2]) != 1 || mxGetNumberOfElements(prhs[3]) != 1)
        mexErrMsgIdAndTxt("KSG_mi:size", "k and algorithm must be scalars.");
    kd = mxGetScalar(prhs[2]);
    ad = mxGetScalar(prhs[3]);
    if (!(kd >= 1) || kd != floor(kd))              /* (NaN fails too) */
        mexErrMsgIdAndTxt("KSG_mi:k", "k must be a positive integer.");
    if (ad != 1 && ad != 2)
        mexErrMsgIdAndTxt("KSG_mi:algorithm", "algorithm must be 1 or 2.");
    if (kd >= (double)n)                            /* before the cast: k may be huge */
        mexErrMsgIdAndTxt("KSG_mi:tooFew", "%.0f observations: need more than k = %g.",
                          (double)n, kd);
    k = (mwSize)kd;
    alg = (int)ad;
    X = mxGetPr(prhs[0]);
    Y = mxGetPr(prhs[1]);
    D = dx + 1;

    /* normalised columns, in the order of the first source column */
    double **z = (double **)mxMalloc(D * sizeof(double *));
    double *col = (double *)mxMalloc(n * sizeof(double));
    Pair *order = (Pair *)mxMalloc(n * sizeof(Pair));
    for (d = 0; d < D; d++) {
        const double *src = d < dx ? X + d * n : Y;
        double m = 0, s = 0;
        for (i = 0; i < n; i++) {
            if (!mxIsFinite(src[i]))
                mexErrMsgIdAndTxt("KSG_mi:nonFinite", "The data contain NaN or Inf.");
            m += src[i];
        }
        m /= n;
        for (i = 0; i < n; i++) s += (src[i] - m) * (src[i] - m);
        s = n > 1 ? sqrt(s / (n - 1)) : 0;
        z[d] = (double *)mxMalloc(n * sizeof(double));
        for (i = 0; i < n; i++) z[d][i] = s > 0 ? (src[i] - m) / s : src[i] - m;
    }
    for (i = 0; i < n; i++) { order[i].v = z[0][i]; order[i].i = i; }
    qsort(order, n, sizeof(Pair), comparePairs);
    for (d = 0; d < D; d++) {
        for (i = 0; i < n; i++) col[i] = z[d][order[i].i];
        memcpy(z[d], col, n * sizeof(double));
    }
    double *ys = (double *)mxMalloc(n * sizeof(double));   /* target, sorted */
    memcpy(ys, z[dx], n * sizeof(double));
    qsort(ys, n, sizeof(double), compareDoubles);

    double *psi = (double *)mxMalloc((n + 1) * sizeof(double));   /* psi(1..n) */
    psi[1] = -0.57721566490153286061;
    for (i = 1; i < n; i++) psi[i + 1] = psi[i] + 1.0 / (double)i;

    double *bd = (double *)mxMalloc(k * sizeof(double));
    mwSize *bt = (mwSize *)mxMalloc(k * sizeof(mwSize));
    double *bx = (double *)mxMalloc(k * sizeof(double));
    double *by = (double *)mxMalloc(k * sizeof(double));
    const double *z0 = z[0], *zy = z[dx];
    double sum = 0;

    for (p = 0; p < n; p++) {
        /* k nearest neighbours in the joint space */
        nb = 0;
        for (q = p + 1; q < n; q++) {
            double d0 = z0[q] - z0[p];
            if (nb == k && d0 > bd[k - 1]) break;
            double ex = d0;
            for (d = 1; d < dx; d++) { double t = fabs(z[d][q] - z[d][p]); if (t > ex) ex = t; }
            double ey = fabs(zy[q] - zy[p]);
            keep(ex > ey ? ex : ey, order[q].i, ex, ey, k, &nb, bd, bt, bx, by);
        }
        for (q = p; q-- > 0; ) {
            double d0 = z0[p] - z0[q];
            if (nb == k && d0 > bd[k - 1]) break;
            double ex = d0;
            for (d = 1; d < dx; d++) { double t = fabs(z[d][q] - z[d][p]); if (t > ex) ex = t; }
            double ey = fabs(zy[q] - zy[p]);
            keep(ex > ey ? ex : ey, order[q].i, ex, ey, k, &nb, bd, bt, bx, by);
        }

        /* the radii, and the points within them in each marginal space */
        double rx, ry;
        int strict = alg == 1;
        if (alg == 1) rx = ry = bd[k - 1];
        else {
            rx = ry = 0;
            for (i = 0; i < k; i++) { if (bx[i] > rx) rx = bx[i]; if (by[i] > ry) ry = by[i]; }
        }
        mwSize nx, ny;
        window(z0, n, z0[p], rx, strict, &a, &b);
        if (dx == 1) nx = b - a;
        else {
            nx = 0;
            for (q = a; q < b; q++) {
                int in = 1;
                for (d = 1; d < dx && in; d++) {
                    double t = fabs(z[d][q] - z[d][p]);
                    in = strict ? t < rx : t <= rx;
                }
                nx += in;
            }
        }
        window(ys, n, zy[p], ry, strict, &a, &b);
        ny = b - a;
        /* not the point itself (at distance 0: counted unless strictly within 0) */
        if (!strict || rx > 0) nx--;
        if (!strict || ry > 0) ny--;
        sum += alg == 1 ? psi[nx + 1] + psi[ny + 1] : psi[nx] + psi[ny];
    }

    double I = alg == 1 ? psi[k] - sum / n + psi[n]
                        : psi[k] - 1.0 / (double)k - sum / n + psi[n];
    plhs[0] = mxCreateDoubleScalar(I);

    for (d = 0; d < D; d++) mxFree(z[d]);
    mxFree(z); mxFree(col); mxFree(order); mxFree(ys); mxFree(psi);
    mxFree(bd); mxFree(bt); mxFree(bx); mxFree(by);
}
