/* tisean_mex.h -- make a TISEAN command-line program safe to run inside MATLAB.
 *
 * Force-included (-include) into every TISEAN translation unit of a MEX build.
 * It rewrites three things the programs do that are fine in a process of their
 * own and fatal inside a host process:
 *
 *   exit()    -> longjmp back to the gateway. A TISEAN exit() inside a MEX file
 *                terminates MATLAB. There are 46 of them across the six C
 *                programs and the support routines they link, and the path is
 *                reached in normal use -- false_nearest's "not enough points"
 *                guard is an exit().
 *
 *   malloc()  -> a tracking arena. TISEAN allocates on every path and frees on
 *   calloc()     some of them: d2, nstat_z and false_nearest free their globals
 *   realloc()    at the end of main, boxcount and lyap_r do not, and NO path
 *   free()       frees anything when the program exits early. A per-program
 *                "free the globals" reset would therefore double-free on the
 *                normal path and leak on the early one. Tracking every block
 *                and sweeping what is still outstanding at the end of the run
 *                is correct on both paths and needs no per-program audit.
 *
 * The reset of the file-scope globals themselves is per-program and lives in
 * the individual gateways -- see tisean_d2.c for why it is not optional.
 *
 * Nothing here is used by the standalone executable build, which still compiles
 * the same sources unmodified.
 */
#ifndef TISEAN_MEX_H
#define TISEAN_MEX_H

#include <stdlib.h>
#include <stddef.h>
#include <setjmp.h>

#if defined(_MSC_VER)
#  define TISEAN_TLS __declspec(thread)
#  define TISEAN_NORETURN __declspec(noreturn)
#  define TISEAN_EXPORT __declspec(dllexport)
#elif defined(_WIN32)
#  define TISEAN_TLS __thread
#  define TISEAN_NORETURN __attribute__((noreturn))
#  define TISEAN_EXPORT __declspec(dllexport)
#else
#  define TISEAN_TLS __thread
#  define TISEAN_NORETURN __attribute__((noreturn))
#  define TISEAN_EXPORT __attribute__((visibility("default")))
#endif

/* Why the ELF builds compile with -fvisibility=hidden.
 *
 * It has to be the compiler flag, not `#pragma GCC visibility push(hidden)` in
 * this header. The pragma applies to every declaration that follows it,
 * including the ones in system headers the TISEAN sources include afterwards,
 * so `extern FILE *stderr` became a hidden reference to a symbol defined in
 * libc and the link failed outright. -fvisibility=hidden hides what this object
 * DEFINES and leaves what it merely references alone, which is what is wanted.
 *
 * d2.c, nstat_z.c, boxcount.c, lyap_r.c, false_nearest.c and poincare.c all
 * declare globals with the same names -- series, verbosity, length, list, box,
 * DIM, DELAY -- at file scope. On Windows each MEX file is a DLL and keeps its
 * own copy, but on Linux and macOS the dynamic loader interposes across shared
 * objects: whichever TISEAN MEX MATLAB happened to load first wins, and every
 * other one silently reads and writes its neighbour's state. Measured: with
 * default visibility, eleven interleaved calls succeeded and the twelfth
 * segfaulted inside lyap_r. Hiding everything but mexFunction removes the
 * possibility entirely. Harmless on Windows, where it is already the rule. */

/* Thread-local so that a stray thread-based parallel pool corrupts nothing
   silently; the gateway's re-entrancy guard turns that case into a clean
   MATLAB error rather than a longjmp across stacks. */
extern TISEAN_TLS jmp_buf tisean_jmp;
extern TISEAN_TLS int     tisean_status;

TISEAN_NORETURN void tisean_exit(int code);

void *tisean_malloc(size_t n);
void *tisean_calloc(size_t n, size_t s);
void *tisean_realloc(void *p, size_t n);
void  tisean_free(void *p);

#define exit(c)      tisean_exit(c)
#define malloc(n)    tisean_malloc((size_t)(n))
#define calloc(n,s)  tisean_calloc((size_t)(n),(size_t)(s))
#define realloc(p,n) tisean_realloc((void*)(p),(size_t)(n))
#define free(p)      tisean_free((void*)(p))

#endif /* TISEAN_MEX_H */
