/* MEX gateway for TISEAN nstat_z. See tisean_reset.h for why the reset matters. */
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

extern unsigned int nmax;
extern long **box, *list;
extern unsigned long *found;
extern double *series, *series1, *series2;
extern double interval, min, epsilon;
extern char epsset, causalset;
extern char *infile;
extern char *outfile, stdo, centerset;
extern char *firstwindow, *secondwindow, **window;
extern unsigned int COLUMN, pieces;
extern unsigned int verbosity;
extern int DIM, DELAY, MINN, STEP;
extern int firstoffset, secondoffset;
extern double EPS0, EPSF;
extern unsigned long LENGTH, exclude, center, causal;

#define NSTAT_Z_GLOBALS(X)                                                \
    X(unsigned int, nmax)                                                 \
    X(long**, box)            X(long*, list)                              \
    X(unsigned long*, found)                                              \
    X(double*, series)        X(double*, series1)  X(double*, series2)    \
    X(double, interval)       X(double, min)       X(double, epsilon)     \
    X(char, epsset)           X(char, causalset)                          \
    X(char*, infile)          X(char*, outfile)                           \
    X(char, stdo)             X(char, centerset)                          \
    X(char*, firstwindow)     X(char*, secondwindow) X(char**, window)    \
    X(unsigned int, COLUMN)   X(unsigned int, pieces)                     \
    X(unsigned int, verbosity)                                            \
    X(int, DIM)  X(int, DELAY)  X(int, MINN)  X(int, STEP)                \
    X(int, firstoffset)       X(int, secondoffset)                        \
    X(double, EPS0)           X(double, EPSF)                             \
    X(unsigned long, LENGTH)  X(unsigned long, exclude)                   \
    X(unsigned long, center)  X(unsigned long, causal)

TISEAN_DEFINE_RESET(NSTAT_Z_GLOBALS)

int tisean_nstat_z_main(int argc, char **argv);
void tisean_mexmain(int, mxArray**, int, const mxArray**, const char*,
                    void (*)(void), int (*)(int, char**));

TISEAN_EXPORT void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    tisean_mexmain(nlhs, plhs, nrhs, prhs, "nstat_z", reset, tisean_nstat_z_main);
}
