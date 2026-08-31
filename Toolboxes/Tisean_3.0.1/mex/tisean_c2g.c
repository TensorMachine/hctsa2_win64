/* MEX gateway for TISEAN c2g (Fortran).
 *
 * HCTSA2_BuildTiseanMex generates a subroutine form of c2g.f at build time --
 * gfortran compiles an unnamed main program to a LOCAL MAIN__ plus a global
 * main, so the program body is not callable from outside until it is turned
 * into a subroutine. The generated copy also routes every `stop` back here.
 *
 * _gfortran_set_args is what makes iargc()/getarg() -- and therefore every
 * TISEAN option -- read from the argv we build rather than from MATLAB's own
 * command line.
 */
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

void _gfortran_set_args(int argc, char **argv);
void c2g_main_(void);
void tsreset_(void);

/* Kept so libgfortran never holds a pointer to an argv the gateway has freed. */
static char  ts_idle_arg[] = "c2g";
static char *ts_idle[2] = { ts_idle_arg, 0 };

static void reset(void)
{
    tsreset_();
    _gfortran_set_args(1, ts_idle);
}

static int entry(int argc, char **argv)
{
    _gfortran_set_args(argc, argv);
    c2g_main_();
    return 0;
}

void tisean_mexmain(int, mxArray**, int, const mxArray**, const char*,
                    void (*)(void), int (*)(int, char**));

TISEAN_EXPORT void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    tisean_mexmain(nlhs, plhs, nrhs, prhs, "c2g", reset, entry);
}
