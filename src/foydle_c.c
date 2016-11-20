#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

static const char **extract_colnames(SEXP mat);
static SEXP compute_and_save(double *xmat, double *ymat, double *zmat,
    int xcol, int ycol, int zcol, int n, const char **xnames,
    const char **ynames, const char **znames, SEXP colnames,
    const char *filename, double rvalue_threshold, int cores,
    int with_return, int swap_y_and_z);
static void print_header(FILE *fp, SEXP colnames);
static void print_rvalues(FILE *fp, int *xindex, int *yindex, int *zindex,
    const char **xnames, const char **ynames, const char **znames,
    double *rvalues, int offset, int nsignif);
static int annotate_rvalues(double *rvalue, int n, double threshold,
    int *xindex, int *yindex, int *zindex, double *rvalues, int xcol,
    int zi, int swap_y_and_z, int offset);
static SEXP create_data_frame(int *xindex, int *yindex, int *zindex,
    const char **xnames, const char **ynames, const char **znames,
    double *rvalues, int offset, SEXP colnames);

void F77_NAME(rval)(double *, double *, double *, int *, int *, int *, double *, int *);
void F77_NAME(center)(double *, int *, int *);

SEXP compute_and_save_rvalues(SEXP xmat_, SEXP ymat_, SEXP zmat_, SEXP n,
    SEXP output_file, SEXP colnames, SEXP rvalue_threshold, SEXP cores,
    SEXP with_return)
{
    int swap_y_and_z = ncols(zmat_) > ncols(ymat_);
    SEXP xmat = PROTECT(duplicate(xmat_));
    SEXP ymat = PROTECT(duplicate(swap_y_and_z ? zmat_ : ymat_));
    SEXP zmat = PROTECT(duplicate(swap_y_and_z ? ymat_ : zmat_));

    const char **xnames = extract_colnames(xmat_);
    const char **ynames = extract_colnames(ymat_);
    const char **znames = extract_colnames(zmat_);


    SEXP result = compute_and_save(REAL(xmat), REAL(ymat), REAL(zmat),
        ncols(xmat), ncols(ymat), ncols(zmat), asInteger(n),
        xnames, ynames, znames, colnames,
        (output_file == R_NilValue) ? NULL : CHAR(asChar(output_file)),
        asReal(rvalue_threshold), asInteger(cores), asInteger(with_return),
        swap_y_and_z);

    free(xnames);
    free(ynames);
    free(znames);

    UNPROTECT(3);
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
    const char **ynames, const char **znames, SEXP colnames,
    const char *filename, double rvalue_threshold, int cores,
    int with_return, int swap_y_and_z)
{
    int *xindex = NULL, *yindex = NULL, *zindex = NULL;
    double *rvalues = NULL;

    int nelem = with_return ? xcol * ycol * zcol : xcol * ycol;
    xindex = malloc(nelem * sizeof(*xindex));
    yindex = malloc(nelem * sizeof(*yindex));
    zindex = malloc(nelem * sizeof(*zindex));
    rvalues = malloc(nelem * sizeof(*rvalues));

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
        int nsignif = annotate_rvalues(rvalue, xcol * ycol, rvalue_threshold,
            xindex, yindex, zindex, rvalues, xcol, i, swap_y_and_z,
            with_return ? offset : 0);
        if (filename)
            print_rvalues(fp, xindex, yindex, zindex, xnames, ynames, znames,
                rvalues, with_return ? offset : 0, nsignif);
        offset += nsignif;
    }
    if (filename)
        fclose(fp);
    free(rvalue);

    if (!with_return) {
        free(xindex);
        free(yindex);
        free(zindex);
        free(rvalues);
        return R_NilValue;
    }

    SEXP result = create_data_frame(xindex, yindex, zindex,
        xnames, ynames, znames, rvalues, offset, colnames);

    free(xindex);
    free(yindex);
    free(zindex);
    free(rvalues);

    return result;
}

static void print_header(FILE *fp, SEXP colnames) {
    fprintf(fp, "%s\t%s\t%s\tr\n",
        CHAR(STRING_ELT(colnames, 0)),
        CHAR(STRING_ELT(colnames, 1)),
        CHAR(STRING_ELT(colnames, 2)));
}

static void print_rvalues(FILE *fp, int *xindex, int *yindex, int *zindex,
    const char **xnames, const char **ynames, const char **znames,
    double *rvalues, int offset, int nsignif)
{
    for (int i = offset; i < offset + nsignif; i++) {
        fprintf(fp, "%s\t%s\t%s\t%.9f\n",
            xnames[xindex[i]],
            ynames[yindex[i]],
            znames[zindex[i]],
            rvalues[i]);
    }
}

static int annotate_rvalues(double *rvalue, int n, double threshold,
    int *xindex, int *yindex, int *zindex, double *rvalues, int xcol,
    int zi, int swap_y_and_z, int offset)
{
    int no_threshold = !R_FINITE(threshold);
    int initial_offset = offset;
    for (int i = 0; i < n; i++) {
        if (no_threshold || fabs(rvalue[i]) >= threshold) {
            xindex[offset] = i % xcol;
            yindex[offset] = swap_y_and_z ? zi : i / xcol;
            zindex[offset] = swap_y_and_z ? i / xcol : zi;
            rvalues[offset] = rvalue[i];
            ++offset;
        }
    }
    return offset - initial_offset;
}

static SEXP create_data_frame(int *xindex, int *yindex, int *zindex,
    const char **xnames, const char **ynames, const char **znames,
    double *rvalues, int nrow, SEXP colnames)
{
    SEXP x = PROTECT(allocVector(STRSXP, nrow));
    SEXP y = PROTECT(allocVector(STRSXP, nrow));
    SEXP z = PROTECT(allocVector(STRSXP, nrow));
    SEXP r = PROTECT(allocVector(REALSXP, nrow));
    for (int i = 0; i < nrow; i++) {
        SET_STRING_ELT(x, i, mkChar(xnames[xindex[i]]));
        SET_STRING_ELT(y, i, mkChar(ynames[yindex[i]]));
        SET_STRING_ELT(z, i, mkChar(znames[zindex[i]]));
        REAL(r)[i] = rvalues[i];
    }

    SEXP result = PROTECT(allocVector(VECSXP, 4));
    SET_VECTOR_ELT(result, 0, x);
    SET_VECTOR_ELT(result, 1, y);
    SET_VECTOR_ELT(result, 2, z);
    SET_VECTOR_ELT(result, 3, r);

    SEXP class = PROTECT(allocVector(STRSXP, 1));
    SET_STRING_ELT(class, 0, mkChar("data.frame"));
    setAttrib(result, R_ClassSymbol, class);

    SEXP names = PROTECT(allocVector(STRSXP, 4));
    SET_STRING_ELT(names, 0, STRING_ELT(colnames, 0));
    SET_STRING_ELT(names, 1, STRING_ELT(colnames, 1));
    SET_STRING_ELT(names, 2, STRING_ELT(colnames, 2));
    SET_STRING_ELT(names, 3, mkChar("r"));
    setAttrib(result, R_NamesSymbol, names);

    SEXP rownames = PROTECT(allocVector(INTSXP, nrow));
    for (int i = 0; i < nrow; i++)
        INTEGER(rownames)[i] = i + 1;
    setAttrib(result, R_RowNamesSymbol, rownames);

    UNPROTECT(8);
    return result;
}
