#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

static const char **extract_colnames(SEXP mat);
static void compute_and_save(double *xmat, double *ymat, double *zmat,
    int xcol, int ycol, int zcol, int n, const char **xnames,
    const char **ynames, const char **znames, const char **colnames,
    const char *filename, double rvalue_threshold, int cores,
    int swap_y_and_z);
static void print_header(FILE *fp, const char **colnames);
static void print_rvalues(FILE *fp, double *rvalue, int n, int xcol, int ycol,
    double threshold, const char **xnames, const char **ynames, const char *zname,
    int swap_y_and_z);

void F77_NAME(rval)(double *, double *, double *, int *, int *, int *, double *, int *);
void F77_NAME(center)(double *, int *, int *);

SEXP compute_and_save_rvalues(SEXP xmat_, SEXP ymat_, SEXP zmat_, SEXP n,
    SEXP output_file, SEXP colnames_, SEXP rvalue_threshold, SEXP cores)
{
    int swap_y_and_z = ncols(zmat_) > ncols(ymat_);

    SEXP xmat = PROTECT(duplicate(xmat_));
    SEXP ymat = PROTECT(duplicate(swap_y_and_z ? zmat_ : ymat_));
    SEXP zmat = PROTECT(duplicate(swap_y_and_z ? ymat_ : zmat_));

    const char **xnames = extract_colnames(xmat_);
    const char **ynames = extract_colnames(swap_y_and_z ? zmat_ : ymat_);
    const char **znames = extract_colnames(swap_y_and_z ? ymat_ : zmat_);

    const char *colnames[] = {
        CHAR(STRING_ELT(colnames_, 0)),
        CHAR(STRING_ELT(colnames_, 1)),
        CHAR(STRING_ELT(colnames_, 2))
    };

    compute_and_save(REAL(xmat), REAL(ymat), REAL(zmat),
        ncols(xmat), ncols(ymat), ncols(zmat), asInteger(n),
        xnames, ynames, znames, colnames, CHAR(asChar(output_file)),
        asReal(rvalue_threshold), asInteger(cores), swap_y_and_z);

    UNPROTECT(3);
    free(xnames);
    free(ynames);
    free(znames);

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

static void compute_and_save(double *xmat, double *ymat, double *zmat,
    int xcol, int ycol, int zcol, int n, const char **xnames,
    const char **ynames, const char **znames, const char **colnames,
    const char *filename, double rvalue_threshold, int cores,
    int swap_y_and_z)
{
    F77_CALL(center)(xmat, &n, &xcol);
    F77_CALL(center)(ymat, &n, &ycol);
    F77_CALL(center)(zmat, &n, &zcol);

    double *rvalue = malloc(xcol * ycol * sizeof(*rvalue));
    FILE *fp = fopen(filename, "wb");
    print_header(fp, colnames);
    for (int i = 0; i < zcol; i++) {
        F77_CALL(rval)(xmat, ymat, zmat + i * n, &xcol, &ycol, &n, rvalue, &cores);
        print_rvalues(fp, rvalue, xcol * ycol, xcol, ycol,
            rvalue_threshold, xnames, ynames, znames[i], swap_y_and_z);
    }
    fclose(fp);
    free(rvalue);
}

static void print_header(FILE *fp, const char **colnames) {
    fprintf(fp, "%s\t%s\t%s\tr\n", colnames[0], colnames[1], colnames[2]);
}

static void print_rvalues(FILE *fp, double *rvalue, int n, int xcol, int ycol,
    double threshold, const char **xnames, const char **ynames, const char *zname,
    int swap_y_and_z)
{
    int no_threshold = !R_FINITE(threshold);
    for (int i = 0; i < n; i++) {
        if (no_threshold || fabs(rvalue[i]) >= threshold) {
            const char *x = xnames[i % xcol];
            const char *y = swap_y_and_z ? zname : ynames[i / xcol];
            const char *z = swap_y_and_z ? ynames[i / xcol] : zname;
            fprintf(fp, "%s\t%s\t%s\t%.9f\n", x, y, z, rvalue[i]);
        }
    }
}
