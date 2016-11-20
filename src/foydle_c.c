#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

static const char **extract_colnames(SEXP mat);
static SEXP compute_and_save(double *xmat, double *ymat, double *zmat,
    int xcol, int ycol, int zcol, int n, const char **xnames,
    const char **ynames, const char **znames, SEXP xnames2, SEXP ynames2,
    SEXP znames2, const char **colnames, const char *filename,
    double rvalue_threshold, int cores, int swap_y_and_z);
static void print_header(FILE *fp, const char **colnames);
static void print_rvalues(FILE *fp, double *rvalue, int n, int xcol, int ycol,
    double threshold, const char **xnames, const char **ynames, const char *zname,
    int swap_y_and_z);
static int annotate_rvalues(double *rvalue, int n, double threshold, int *xindex, int *yindex, int *zindex, double *rvalues, int xcol, int zi, int swap_y_and_z, int offset);

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

    SEXP xnames2 = VECTOR_ELT(getAttrib(xmat_, R_DimNamesSymbol), 1);
    SEXP ynames2 = VECTOR_ELT(getAttrib(ymat_, R_DimNamesSymbol), 1);
    SEXP znames2 = VECTOR_ELT(getAttrib(zmat_, R_DimNamesSymbol), 1);

    const char *colnames[] = {
        CHAR(STRING_ELT(colnames_, 0)),
        CHAR(STRING_ELT(colnames_, 1)),
        CHAR(STRING_ELT(colnames_, 2))
    };

    SEXP result = compute_and_save(REAL(xmat), REAL(ymat), REAL(zmat),
        ncols(xmat), ncols(ymat), ncols(zmat), asInteger(n),
        xnames, ynames, znames, xnames2, ynames2, znames2, colnames,
        (output_file == R_NilValue) ? NULL : CHAR(asChar(output_file)),
        asReal(rvalue_threshold), asInteger(cores), swap_y_and_z);

    UNPROTECT(3);
    free(xnames);
    free(ynames);
    free(znames);

    return result;
}

static const char **extract_colnames(SEXP mat) {
    SEXP colnames = VECTOR_ELT(getAttrib(mat, R_DimNamesSymbol), 1);
    int n = length(colnames);
    const char **result = malloc(n * sizeof(*result));
    for (int i = 0; i < n; i++)
        result[i] = CHAR(STRING_ELT(colnames, i));
    return result;
}

static SEXP compute_and_save(double *xmat, double *ymat, double *zmat,
    int xcol, int ycol, int zcol, int n, const char **xnames,
    const char **ynames, const char **znames, SEXP xnames2, SEXP ynames2,
    SEXP znames2, const char **colnames, const char *filename,
    double rvalue_threshold, int cores, int swap_y_and_z)
{
    int *xindex = NULL, *yindex = NULL, *zindex = NULL;
    double *rvalues = NULL;

    if (!filename) {
        xindex = malloc(xcol * ycol * zcol * sizeof(*xindex));
        yindex = malloc(xcol * ycol * zcol * sizeof(*yindex));
        zindex = malloc(xcol * ycol * zcol * sizeof(*zindex));
        rvalues = malloc(xcol * ycol * zcol * sizeof(*rvalues));
    }

    F77_CALL(center)(xmat, &n, &xcol);
    F77_CALL(center)(ymat, &n, &ycol);
    F77_CALL(center)(zmat, &n, &zcol);

    FILE *fp = NULL;
    if (filename) {
        fp = fopen(filename, "wb");
        print_header(fp, colnames);
    }
    double *rvalue = malloc(xcol * ycol * sizeof(*rvalue));
    int offset = 0;
    for (int i = 0; i < zcol; i++) {
        F77_CALL(rval)(xmat, ymat, zmat + i * n, &xcol, &ycol, &n, rvalue, &cores);
        if (filename)
            print_rvalues(fp, rvalue, xcol * ycol, xcol, ycol, rvalue_threshold, xnames, ynames, znames[i], swap_y_and_z);
        if (!filename)
            offset += annotate_rvalues(rvalue, xcol * ycol, rvalue_threshold, xindex, yindex, zindex, rvalues, xcol, i, swap_y_and_z, offset);
    }
    if (filename)
        fclose(fp);
    free(rvalue);

    if (filename)
        return R_NilValue;

    int number_of_results = offset;
    SEXP x = PROTECT(allocVector(STRSXP, number_of_results));
    SEXP y = PROTECT(allocVector(STRSXP, number_of_results));
    SEXP z = PROTECT(allocVector(STRSXP, number_of_results));
    SEXP r = PROTECT(allocVector(REALSXP, number_of_results));
    for (int i = 0; i < number_of_results; i++) {
        SET_STRING_ELT(x, i, STRING_ELT(xnames2, xindex[i]));
        SET_STRING_ELT(y, i, STRING_ELT(ynames2, yindex[i]));
        SET_STRING_ELT(z, i, STRING_ELT(znames2, zindex[i]));
        REAL(r)[i] = rvalues[i];
    }
    free(xindex);
    free(yindex);
    free(zindex);
    free(rvalues);
    SEXP result = PROTECT(allocVector(VECSXP, 4));
    SET_VECTOR_ELT(result, 0, x);
    SET_VECTOR_ELT(result, 1, y);
    SET_VECTOR_ELT(result, 2, z);
    SET_VECTOR_ELT(result, 3, r);
    SEXP class = PROTECT(allocVector(STRSXP, 1));
    SET_STRING_ELT(class, 0, mkChar("data.frame"));
    setAttrib(result, R_ClassSymbol, class);
    SEXP names = PROTECT(allocVector(STRSXP, 4));
    SET_STRING_ELT(names, 0, mkChar(colnames[0]));
    SET_STRING_ELT(names, 1, mkChar(colnames[1]));
    SET_STRING_ELT(names, 2, mkChar(colnames[2]));
    SET_STRING_ELT(names, 3, mkChar("r"));
    setAttrib(result, R_NamesSymbol, names);
    SEXP rownames = PROTECT(allocVector(INTSXP, number_of_results));
    for (int i = 0; i < number_of_results; i++)
        INTEGER(rownames)[i] = i + 1;
    setAttrib(result, R_RowNamesSymbol, rownames);
    UNPROTECT(8);

    return result;
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

static int annotate_rvalues(double *rvalue, int n, double threshold, int *xindex, int *yindex, int *zindex, double *rvalues, int xcol, int zi, int swap_y_and_z, int offset)
{
    int no_threshold = !R_FINITE(threshold);
    int result = 0;
    for (int i = 0; i < n; i++) {
        if (no_threshold || fabs(rvalue[i]) >= threshold) {
            xindex[i + offset] = i % xcol;
            yindex[i + offset] = swap_y_and_z ? zi : i / xcol;
            zindex[i + offset] = swap_y_and_z ? i / xcol : zi;
            rvalues[i + offset] = rvalue[i];
            ++result;
        }
    }
    return result;
}
