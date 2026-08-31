/* MEX gateway for TISEAN lyap_r. See tisean_reset.h for why the reset matters. */
/* mex.h before tisean_mex.h: that header pushes hidden visibility, and a
   mexFunction declared under it stays hidden however this definition is
   marked -- the MEX file would load with no entry point. */
#include "mex.h"
#include "tisean_mex.h"
#undef exit
#undef malloc
#undef calloc
#undef realloc
#undef free
#include "tisean_reset.h"

extern char *outfile;
extern char *infile;
extern char epsset;
extern double *series, *lyap;
extern long *list;
extern unsigned int dim, delay, steps, mindist;
extern unsigned int column;
extern unsigned int verbosity;
extern unsigned long length, exclude;
extern long *found;
extern double eps0, eps, epsinv;

/* box[NMAX][NMAX] is set to -1 throughout by put_in_boxes on every run.
   nmax is const. */
#define LYAP_R_GLOBALS(X)                                                 \
    X(char*, outfile)         X(char*, infile)                            \
    X(char, epsset)                                                       \
    X(double*, series)        X(double*, lyap)                            \
    X(long*, list)            X(long*, found)                             \
    X(unsigned int, dim)      X(unsigned int, delay)                      \
    X(unsigned int, steps)    X(unsigned int, mindist)                    \
    X(unsigned int, column)   X(unsigned int, verbosity)                  \
    X(unsigned long, length)  X(unsigned long, exclude)                   \
    X(double, eps0)  X(double, eps)  X(double, epsinv)

TISEAN_DEFINE_RESET(LYAP_R_GLOBALS)

int tisean_lyap_r_main(int argc, char **argv);
void tisean_mexmain(int, mxArray**, int, const mxArray**, const char*,
                    void (*)(void), int (*)(int, char**));

TISEAN_EXPORT void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    tisean_mexmain(nlhs, plhs, nrhs, prhs, "lyap_r", reset, tisean_lyap_r_main);
}
