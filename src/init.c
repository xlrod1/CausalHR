#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: I think it might be fixed
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls  --- what was in the skeleton from tools::package_native_routine_registration_skeleton
extern void F77_NAME(knnhad)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(newhad)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
*/ 

/* #   my guess at what is needed here 
     Not sure whether the asterisks are correct
     Basing the int and double designations on the R calling routine
*/
/* .Fortran calls */
extern void F77_NAME(knnhad)(int *n, double *x, int *delta, int *ks, int *bwchoi, int *gridz, double *z, int *m, double *zz, double *bpilot, double *endl, double *endr, double *bsmo, int *kflag, double *fzz, int *kmin, int *kmax, double *bopt, double *bopt1, double *kimse);
extern void F77_NAME(newhad)(int *n, double *x, int *delta, int *ks, int *local, double *z, int *gridz,
double *zz, int *m, double *bpilot, double *bw, int *gridb, double *endl, double *endr, double *bsmo,  int *kflag, double *fzz, double *bopt, double *bopt1, double *msemin, double *biasmn, double *varmin, double * imsemn, double *globlb, double *glmse,double *w, double *naes, double *nw);

/* Not sure if I should be doing anything with "NULL"s  */
static const R_FortranMethodDef FortranEntries[] = {
    {"knnhad", (DL_FUNC) &F77_NAME(knnhad), 20},
    {"newhad", (DL_FUNC) &F77_NAME(newhad), 28},
    {NULL, NULL, 0}
};

void R_init_muhaz(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

/*  And now to figure out exactly where this should go ??? 
Writing R Extensions seems to say it should be in a file named init.c in the src/ directory
So that will be my first try
*/
