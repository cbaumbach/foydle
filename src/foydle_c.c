#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <stdio.h>

static const char **extract_colnames(SEXP mat);

void F77_NAME(rval)(double *, double *, double *, int *, int *, int *, double *);

SEXP compute_and_save_rvalues(SEXP xmat, SEXP ymat, SEXP zmat, SEXP n_, SEXP output_file)
{
    int n = asInteger(n_);
    int xcol = length(xmat) / n;
    int ycol = length(ymat) / n;
    int zcol = length(zmat) / n;
    const char **xnames = extract_colnames(xmat);
    const char **ynames = extract_colnames(ymat);
    const char **znames = extract_colnames(zmat);

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
    free(xnames); free(ynames); free(znames);

    return R_NilValue;
}

static const char **extract_colnames(SEXP mat) {
    SEXP colnames = VECTOR_ELT(getAttrib(mat, R_DimNamesSymbol), 1);
    int n = length(colnames);
    const char **result = malloc(n * sizeof(*result));
    for (int i = 0; i < n; i++)
        result[i] = CHAR(STRING_ELT(colnames, i));
    return result;
}
