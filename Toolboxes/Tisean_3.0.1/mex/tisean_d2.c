/* MEX gateway for TISEAN d2. See tisean_reset.h for why the reset matters. */
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

extern double **series;
extern long *scr;
extern char dimset, rescale_set, eps_min_set, eps_max_set;
extern char *FOUT;
extern double epsfactor, epsinv, lneps, lnfac;
extern double EPSMAX, EPSMIN;
extern double min, interval;
extern int imax, howoften1, imin;
extern long *list, *listc1;
extern unsigned long nmax;
extern double **found, *norm;
extern unsigned long MINDIST, MAXFOUND;
extern unsigned long length, exclude;
extern unsigned int DIM, EMBED, HOWOFTEN, DELAY;
extern unsigned int verbosity;
extern char *column;
extern char *infile;

/* box[NMAX][NMAX] and boxc1[NMAX] are refilled by make_c2 on every run. */
#define D2_GLOBALS(X)                                                     \
    X(double**, series)      X(long*, scr)                                \
    X(char, dimset)          X(char, rescale_set)                         \
    X(char, eps_min_set)     X(char, eps_max_set)                         \
    X(char*, FOUT)                                                        \
    X(double, epsfactor)     X(double, epsinv)                            \
    X(double, lneps)         X(double, lnfac)                             \
    X(double, EPSMAX)        X(double, EPSMIN)                            \
    X(double, min)           X(double, interval)                          \
    X(int, imax)             X(int, howoften1)      X(int, imin)          \
    X(long*, list)           X(long*, listc1)                             \
    X(unsigned long, nmax)                                                \
    X(double**, found)       X(double*, norm)                             \
    X(unsigned long, MINDIST) X(unsigned long, MAXFOUND)                  \
    X(unsigned long, length)  X(unsigned long, exclude)                   \
    X(unsigned int, DIM)      X(unsigned int, EMBED)                      \
    X(unsigned int, HOWOFTEN) X(unsigned int, DELAY)                      \
    X(unsigned int, verbosity)                                            \
    X(char*, column)          X(char*, infile)

TISEAN_DEFINE_RESET(D2_GLOBALS)

int tisean_d2_main(int argc, char **argv);
void tisean_mexmain(int, mxArray**, int, const mxArray**, const char*,
                    void (*)(void), int (*)(int, char**));

TISEAN_EXPORT void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    tisean_mexmain(nlhs, plhs, nrhs, prhs, "d2", reset, tisean_d2_main);
}
