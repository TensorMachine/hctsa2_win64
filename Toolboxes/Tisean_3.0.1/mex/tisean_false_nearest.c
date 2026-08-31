/* MEX gateway for TISEAN false_nearest. See tisean_reset.h for why the reset matters. */
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
extern char stdo, dimset;
extern char *column;
extern unsigned long length, exclude, theiler;
extern unsigned int delay, maxdim, minemb;
extern unsigned int comp, maxemb;
extern unsigned int verbosity;
extern double rt;
extern double eps0;
extern double **series;
extern double aveps, vareps;
extern double varianz;
extern int ibox;
extern long **box, *list;
extern unsigned int *vcomp, *vemb;
extern unsigned long toolarge;

#define FALSE_NEAREST_GLOBALS(X)                                          \
    X(char*, outfile)         X(char*, infile)                            \
    X(char, stdo)             X(char, dimset)                             \
    X(char*, column)                                                      \
    X(unsigned long, length)  X(unsigned long, exclude)                   \
    X(unsigned long, theiler)                                             \
    X(unsigned int, delay)    X(unsigned int, maxdim) X(unsigned int, minemb) \
    X(unsigned int, comp)     X(unsigned int, maxemb)                     \
    X(unsigned int, verbosity)                                            \
    X(double, rt)             X(double, eps0)                             \
    X(double**, series)                                                   \
    X(double, aveps)          X(double, vareps)     X(double, varianz)    \
    X(int, ibox)                                                          \
    X(long**, box)            X(long*, list)                              \
    X(unsigned int*, vcomp)   X(unsigned int*, vemb)                      \
    X(unsigned long, toolarge)

TISEAN_DEFINE_RESET(FALSE_NEAREST_GLOBALS)

int tisean_false_nearest_main(int argc, char **argv);
void tisean_mexmain(int, mxArray**, int, const mxArray**, const char*,
                    void (*)(void), int (*)(int, char**));

TISEAN_EXPORT void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    tisean_mexmain(nlhs, plhs, nrhs, prhs, "false_nearest", reset, tisean_false_nearest_main);
}
