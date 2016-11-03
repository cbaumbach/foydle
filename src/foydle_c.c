#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <stdio.h>

static const char **extract_names(SEXP names, int n);

void F77_NAME(rval)(double *, double *, double *, int *, int *, int *, double *);

SEXP compute_and_save_rvalues(
    SEXP xmat, SEXP ymat, SEXP zmat,
    SEXP xcol_, SEXP ycol_, SEXP zcol_,
    SEXP xnames_, SEXP ynames_, SEXP znames_,
    SEXP n_, SEXP output_file)
{
    int xcol = asInteger(xcol_);
    int ycol = asInteger(ycol_);
    int zcol = asInteger(zcol_);
    const char **xnames = extract_names(xnames_, xcol);
    const char **ynames = extract_names(ynames_, ycol);
    const char **znames = extract_names(znames_, zcol);
    int n = asInteger(n_);

    FILE *fp = fopen(CHAR(asChar(output_file)), "wb");
    fprintf(fp, "x\ty\tz\tr\n");
    double *r = malloc(xcol * ycol * sizeof(*r));
    for (int i = 0; i < zcol; i++) {
        double *z = REAL(zmat) + i * n;
        F77_CALL(rval)(REAL(xmat), REAL(ymat), z, &xcol, &ycol, &n, r);
        for (int j = 0; j < xcol * ycol; j++)
            fprintf(fp, "%s\t%s\t%s\t%.9f\n",
                xnames[j % xcol], ynames[j / xcol], znames[i], r[j]);
    }
    free(r);
    fclose(fp);

    return R_NilValue;
}

static const char **extract_names(SEXP names, int n) {
    const char **result = malloc(n * sizeof(*result));
    for (int i = 0; i < n; i++)
        result[i] = CHAR(STRING_ELT(names, i));
    return result;
}
