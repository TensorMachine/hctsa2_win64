/* MEX gateway for TISEAN boxcount. See tisean_reset.h for why the reset matters. */
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

extern unsigned long LENGTH, exclude;
extern unsigned int maxembed, dimension, DELAY, EPSCOUNT;
extern unsigned int verbosity;
extern double Q, EPSMIN, EPSMAX;
extern char dimset, epsminset, epsmaxset;
extern char *outfile;
extern char *column;
extern int epsi;
extern unsigned long length;
extern double EPSFAKTOR;
extern unsigned int **which_dims;
extern double *histo;
extern double **series;

#define BOXCOUNT_GLOBALS(X)                                               \
    X(unsigned long, LENGTH)  X(unsigned long, exclude)                   \
    X(unsigned int, maxembed) X(unsigned int, dimension)                  \
    X(unsigned int, DELAY)    X(unsigned int, EPSCOUNT)                   \
    X(unsigned int, verbosity)                                            \
    X(double, Q)  X(double, EPSMIN)  X(double, EPSMAX)                    \
    X(char, dimset)  X(char, epsminset)  X(char, epsmaxset)               \
    X(char*, outfile)         X(char*, column)                            \
    X(int, epsi)              X(unsigned long, length)                    \
    X(double, EPSFAKTOR)                                                  \
    X(unsigned int**, which_dims)                                         \
    X(double*, histo)         X(double**, series)

TISEAN_DEFINE_RESET(BOXCOUNT_GLOBALS)

int tisean_boxcount_main(int argc, char **argv);
void tisean_mexmain(int, mxArray**, int, const mxArray**, const char*,
                    void (*)(void), int (*)(int, char**));

TISEAN_EXPORT void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    tisean_mexmain(nlhs, plhs, nrhs, prhs, "boxcount", reset, tisean_boxcount_main);
}
