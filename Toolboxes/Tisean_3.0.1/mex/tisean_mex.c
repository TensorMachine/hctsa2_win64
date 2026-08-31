/* tisean_mex.c -- shared gateway runtime for the TISEAN MEX files.
 * See tisean_mex.h for why the exit() and allocator shims exist.
 */
#include "mex.h"
#include "tisean_mex.h"

/* This file implements the shims, so it must see the real functions. */
#undef exit
#undef malloc
#undef calloc
#undef realloc
#undef free

#include <string.h>
#include <stdio.h>

#ifdef _WIN32
#  include <io.h>
#  define TS_DUP    _dup
#  define TS_DUP2   _dup2
#  define TS_CLOSE  _close
#  define TS_FILENO _fileno
#else
#  include <unistd.h>
#  define TS_DUP    dup
#  define TS_DUP2   dup2
#  define TS_CLOSE  close
#  define TS_FILENO fileno
#endif

TISEAN_TLS jmp_buf tisean_jmp;
TISEAN_TLS int     tisean_status;

void tisean_exit(int code)
{
    tisean_status = code;
    longjmp(tisean_jmp, 1);
}

/* Fortran entry point for the same thing. The generated MEX copies of the
   Fortran sources call this where the originals say `stop`. Verified that a
   longjmp out of a gfortran frame unwinds cleanly: these routines hold no
   allocatables and no finalizable state, and tisean_reset closes the i/o
   units afterwards. */
void tsfatal_(int *code) { tisean_exit(*code); }
void tsfatal(int *code)  { tisean_exit(*code); }   /* -fno-underscoring builds */

/*---------------------------------------------------------------------------
 * Allocation arena
 *
 * Every block carries a 16-byte header holding its registry slot, so free()
 * is O(1) and a block freed by TISEAN itself is never swept a second time.
 * The sweep at the end of a run releases whatever the program did not.
 *-------------------------------------------------------------------------*/
#define TS_HEAD 16

static void   **ts_reg;        /* slot -> block base, NULL if already freed */
static size_t  *ts_hole;       /* stack of reusable slots */
static size_t   ts_reg_n, ts_reg_cap, ts_hole_n, ts_hole_cap;

static int ts_grow(void)
{
    if (ts_reg_n == ts_reg_cap) {
        size_t cap = ts_reg_cap ? ts_reg_cap * 2 : 64;
        void **p = (void**)realloc(ts_reg, cap * sizeof(void*));
        if (p == NULL) return 0;
        ts_reg = p; ts_reg_cap = cap;
    }
    return 1;
}

static size_t ts_slot(void)
{
    if (ts_hole_n > 0) return ts_hole[--ts_hole_n];
    if (!ts_grow()) return (size_t)-1;
    return ts_reg_n++;
}

static void ts_release(size_t slot)
{
    ts_reg[slot] = NULL;
    if (ts_hole_n == ts_hole_cap) {
        size_t cap = ts_hole_cap ? ts_hole_cap * 2 : 64;
        size_t *p = (size_t*)realloc(ts_hole, cap * sizeof(size_t));
        if (p == NULL) return;          /* slot is simply not reused */
        ts_hole = p; ts_hole_cap = cap;
    }
    ts_hole[ts_hole_n++] = slot;
}

void *tisean_malloc(size_t n)
{
    size_t slot;
    char *base;
    if (n == 0) n = 1;
    slot = ts_slot();
    if (slot == (size_t)-1) return NULL;
    base = (char*)malloc(n + TS_HEAD);
    if (base == NULL) { ts_release(slot); return NULL; }
    *(size_t*)base = slot;
    ts_reg[slot] = base;
    return base + TS_HEAD;
}

void *tisean_calloc(size_t n, size_t s)
{
    size_t bytes = n * s;
    void *p;
    if (s != 0 && bytes / s != n) return NULL;      /* overflow */
    p = tisean_malloc(bytes);
    if (p != NULL) memset(p, 0, bytes ? bytes : 1);
    return p;
}

void *tisean_realloc(void *p, size_t n)
{
    char *base, *nbase;
    size_t slot;
    if (p == NULL) return tisean_malloc(n);
    base = (char*)p - TS_HEAD;
    slot = *(size_t*)base;
    nbase = (char*)realloc(base, n + TS_HEAD);
    if (nbase == NULL) return NULL;                 /* old block still tracked */
    ts_reg[slot] = nbase;
    return nbase + TS_HEAD;
}

void tisean_free(void *p)
{
    char *base;
    if (p == NULL) return;
    base = (char*)p - TS_HEAD;
    ts_release(*(size_t*)base);
    free(base);
}

static void ts_sweep(void)
{
    size_t i;
    for (i = 0; i < ts_reg_n; i++)
        if (ts_reg[i] != NULL) { free(ts_reg[i]); ts_reg[i] = NULL; }
    ts_reg_n = 0;
    ts_hole_n = 0;
}

/*---------------------------------------------------------------------------
 * stderr capture
 *
 * TISEAN writes its banner, its progress notes and every diagnostic to stderr,
 * and hctsa's callers test the captured text for emptiness as a did-it-run
 * check. stdout is deliberately left alone: the data-producing programs are
 * given -o by BF_RunTisean instead, so no MEX file ever has to touch fd 1.
 *-------------------------------------------------------------------------*/
typedef struct { int saved; FILE *fp; } ts_cap;

static void cap_begin(ts_cap *c, const char *path)
{
    c->saved = -1; c->fp = NULL;
    if (path == NULL || path[0] == '\0') return;
    c->fp = fopen(path, "w+");
    if (c->fp == NULL) return;
    fflush(stderr);
    c->saved = TS_DUP(TS_FILENO(stderr));
    if (c->saved < 0 || TS_DUP2(TS_FILENO(c->fp), TS_FILENO(stderr)) < 0) {
        if (c->saved >= 0) { TS_CLOSE(c->saved); c->saved = -1; }
        fclose(c->fp); c->fp = NULL;                /* capture unavailable */
    }
}

/* Returns an mxArray char row vector; never NULL. */
static mxArray *cap_end(ts_cap *c, const char *path)
{
    char *buf = NULL;
    long n = 0;
    mxArray *out;

    if (c->fp == NULL) return mxCreateString("");
    fflush(stderr);
    if (c->saved >= 0) { TS_DUP2(c->saved, TS_FILENO(stderr)); TS_CLOSE(c->saved); }
    fflush(c->fp);
    if (fseek(c->fp, 0L, SEEK_END) == 0) {
        n = ftell(c->fp);
        if (n > 1 << 20) n = 1 << 20;               /* a runaway banner is not data */
        if (n > 0 && fseek(c->fp, 0L, SEEK_SET) == 0) {
            buf = (char*)malloc((size_t)n + 1);
            if (buf != NULL) { n = (long)fread(buf, 1, (size_t)n, c->fp); buf[n] = '\0'; }
        }
    }
    fclose(c->fp); c->fp = NULL;
    if (path != NULL) remove(path);
    out = mxCreateString(buf != NULL ? buf : "");
    free(buf);
    return out;
}

/*---------------------------------------------------------------------------
 * Gateway
 *-------------------------------------------------------------------------*/
static char **argv_build(const mxArray *a, const char *prog, int *pargc, size_t *pslots)
{
    mwSize n, i;
    char **argv;

    n = (mwSize)mxGetNumberOfElements(a);
    *pslots = (size_t)n + 2;
    *pargc  = -1;
    argv = (char**)calloc(*pslots, sizeof(char*));
    if (argv == NULL) return NULL;

    /* argv[0] is what TISEAN prints in its banner and usage text. */
    argv[0] = (char*)malloc(strlen(prog) + 1);
    if (argv[0] == NULL) return argv;
    strcpy(argv[0], prog);

    for (i = 0; i < n; i++) {
        const mxArray *e = mxGetCell(a, i);
        size_t len;
        char *s;
        if (e == NULL || !mxIsChar(e)) return argv;
        len = (size_t)mxGetNumberOfElements(e) + 1;
        s = (char*)malloc(len);
        if (s == NULL) return argv;
        argv[i + 1] = s;
        if (mxGetString(e, s, (mwSize)len) != 0) return argv;
    }
    *pargc = (int)n + 1;
    return argv;
}

static void argv_free(char **argv, size_t slots)
{
    size_t i;
    if (argv == NULL) return;
    for (i = 0; i < slots; i++) free(argv[i]);
    free(argv);
}

/* One call is in flight per process. Process-based parallel pools give each
   worker its own copy of this file, so this only ever fires if a thread-based
   pool calls in concurrently -- which TISEAN's globals could not survive. */
static int ts_busy;

void tisean_mexmain(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[],
                    const char *prog, void (*reset)(void), int (*entry)(int, char**))
{
    char **argv = NULL;
    char *capPath = NULL;
    int argc = -1, status = 0, bad = 0;
    size_t slots = 0;
    ts_cap cap;

    (void)nlhs;
    cap.saved = -1; cap.fp = NULL;

    if (nrhs < 1 || nrhs > 2 || !mxIsCell(prhs[0]))
        mexErrMsgIdAndTxt("hctsa:tisean:usage",
            "Usage: [status,stderrText] = tisean_%s(argsCell, capturePath)", prog);
    if (nrhs == 2 && !mxIsChar(prhs[1]) && !mxIsEmpty(prhs[1]))
        mexErrMsgIdAndTxt("hctsa:tisean:usage", "capturePath must be a string.");
    if (ts_busy)
        mexErrMsgIdAndTxt("hctsa:tisean:reentrant",
            "tisean_%s was called while already running. The TISEAN programs "
            "keep global state and cannot be run concurrently in one process; "
            "use a process-based parallel pool, not a thread-based one.", prog);

    if (nrhs == 2 && mxIsChar(prhs[1]) && !mxIsEmpty(prhs[1]))
        capPath = mxArrayToString(prhs[1]);

    ts_busy = 1;
    argv = argv_build(prhs[0], prog, &argc, &slots);
    bad = (argv == NULL || argc < 0);

    if (!bad) {
        cap_begin(&cap, capPath);
        reset();                       /* known state before the run */
        tisean_status = 0;
        if (setjmp(tisean_jmp) == 0) {
            entry(argc, argv);         /* normal return == success */
            status = 0;
        } else {
            status = tisean_status;    /* an exit() we intercepted */
        }
        reset();                       /* drop dangling globals */
        ts_sweep();                    /* free whatever the program did not */
    }
    argv_free(argv, slots);
    ts_busy = 0;

    if (bad) {
        mxArray *t = cap_end(&cap, capPath);
        mxDestroyArray(t);
        if (capPath) mxFree(capPath);
        mexErrMsgIdAndTxt("hctsa:tisean:args",
            "tisean_%s: argsCell must be a cell array of character vectors.", prog);
    }

    plhs[0] = mxCreateDoubleScalar((double)status);
    if (nlhs > 1) {
        plhs[1] = cap_end(&cap, capPath);
    } else {
        mxArray *t = cap_end(&cap, capPath);
        mxDestroyArray(t);
    }
    if (capPath) mxFree(capPath);
}
