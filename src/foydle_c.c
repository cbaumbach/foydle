#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

static const char **extract_colnames(SEXP mat);
static void print_header(FILE *fp, SEXP colnames);
static void print_rvalues(FILE *fp, double *rvalue, int n, int xcol, int ycol,
    double threshold, const char **xnames, const char **ynames, const char *zname);

void F77_NAME(rval)(double *, double *, double *, int *, int *, int *, double *);

SEXP compute_and_save_rvalues(SEXP xmat, SEXP ymat, SEXP zmat,
    SEXP n_, SEXP output_file, SEXP colnames, SEXP rvalue_threshold)
{
    int n = asInteger(n_);
    int xcol = length(xmat) / n;
    int ycol = length(ymat) / n;
    int zcol = length(zmat) / n;
    const char **xnames = extract_colnames(xmat);
    const char **ynames = extract_colnames(ymat);
    const char **znames = extract_colnames(zmat);
    double *rvalue = malloc(xcol * ycol * sizeof(*rvalue));

    FILE *fp = fopen(CHAR(asChar(output_file)), "wb");
    print_header(fp, colnames);
    for (int i = 0; i < zcol; i++) {
        F77_CALL(rval)(REAL(xmat), REAL(ymat), REAL(zmat) + i * n,
            &xcol, &ycol, &n, rvalue);
        print_rvalues(fp, rvalue, xcol * ycol, xcol, ycol,
            asReal(rvalue_threshold), xnames, ynames, znames[i]);
    }
    fclose(fp);

    free(rvalue);
    free(xnames);
    free(ynames);
    free(znames);

    return R_NilValue;
}

static void print_header(FILE *fp, SEXP colnames) {
    fprintf(fp, "%s\t%s\t%s\tr\n",
        CHAR(STRING_ELT(colnames, 0)),
        CHAR(STRING_ELT(colnames, 1)),
        CHAR(STRING_ELT(colnames, 2)));
}

static void print_rvalues(FILE *fp, double *rvalue, int n, int xcol, int ycol,
    double threshold, const char **xnames, const char **ynames, const char *zname)
{
    for (int i = 0; i < n; i++) {
        if (!R_FINITE(threshold) || fabs(rvalue[i]) >= threshold)
            fprintf(fp, "%s\t%s\t%s\t%.9f\n",
                xnames[i % xcol], ynames[i / xcol], zname, rvalue[i]);
    }
}

static const char **extract_colnames(SEXP mat) {
    SEXP colnames = VECTOR_ELT(getAttrib(mat, R_DimNamesSymbol), 1);
    int n = length(colnames);
    const char **result = malloc(n * sizeof(*result));
    for (int i = 0; i < n; i++)
        result[i] = CHAR(STRING_ELT(colnames, i));
    return result;
}
