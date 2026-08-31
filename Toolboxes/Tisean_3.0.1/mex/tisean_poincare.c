/* MEX gateway for TISEAN poincare. See tisean_reset.h for why the reset matters. */
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

extern char *outfile, dimset, compset, whereset, stdo;
extern char *infile;
extern unsigned long length, count, exclude;
extern int dim, comp, delay, dir;
extern unsigned int column;
extern unsigned int verbosity;
extern double *series, min, max, average, where;

#define POINCARE_GLOBALS(X)                                               \
    X(char*, outfile)                                                     \
    X(char, dimset)  X(char, compset)  X(char, whereset)  X(char, stdo)   \
    X(char*, infile)                                                      \
    X(unsigned long, length)  X(unsigned long, count)                     \
    X(unsigned long, exclude)                                             \
    X(int, dim)  X(int, comp)  X(int, delay)  X(int, dir)                 \
    X(unsigned int, column)   X(unsigned int, verbosity)                  \
    X(double*, series)                                                    \
    X(double, min)  X(double, max)  X(double, average)  X(double, where)

TISEAN_DEFINE_RESET(POINCARE_GLOBALS)

int tisean_poincare_main(int argc, char **argv);
void tisean_mexmain(int, mxArray**, int, const mxArray**, const char*,
                    void (*)(void), int (*)(int, char**));

TISEAN_EXPORT void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    tisean_mexmain(nlhs, plhs, nrhs, prhs, "poincare", reset, tisean_poincare_main);
}
