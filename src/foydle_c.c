#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>

void F77_NAME(rval)(double *, double *, double *, int *, int *, int *, double *);

SEXP compute_rvalue(SEXP xmat, SEXP ymat, SEXP zmat,
    SEXP xcol_, SEXP ycol_, SEXP zcol_, SEXP n_, SEXP rvalue_)
{
    int xcol = asInteger(xcol_);
    int ycol = asInteger(ycol_);
    int zcol = asInteger(zcol_);
    int n = asInteger(n_);
    double *rvalue = REAL(rvalue_);

    double *r = malloc(xcol * ycol * sizeof(*r));
    for (int i = 0; i < zcol; i++) {
        double *z = REAL(zmat) + i * n;
        F77_CALL(rval)(REAL(xmat), REAL(ymat), z, &xcol, &ycol, &n, r);
        for (int j = 0; j < xcol * ycol; j++)
            rvalue[i * xcol * ycol + j] = r[j];
    }
    free(r);

    return rvalue_;
}
