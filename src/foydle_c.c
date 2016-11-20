#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

static int filter_and_annotate_rvalues(double *rvalue, int n, double threshold,
    const char **xnames, const char **ynames, const char **znames,
    const char **x, const char **y, const char **z, double *r,
    int xcol, int zi, int swap_y_and_z, int initial_offset);
static SEXP compute_and_save(double *xmat, double *ymat, double *zmat,
    int xcol, int ycol, int zcol, int nrow, const char **xnames,
    const char **ynames, const char **znames, SEXP names,
    const char *filename, double rvalue_threshold, int cores,
    int with_return, int swap_y_and_z);
static SEXP create_data_frame(const char **x, const char **y,
    const char **z, double *r, int offset, SEXP names);
static const char **extract_colnames(SEXP mat);
static void print_header(FILE *fp, SEXP names);
static void print_rvalues(FILE *fp, const char **x, const char **y,
    const char **z, double *r, int offset, int nsignif);

void F77_NAME(rval)(double *xmat, double *ymat, double *z, int *xcol, int *ycol, int *nrow, double *rvalue, int *cores);
void F77_NAME(center)(double *matrix, int *nrow, int *ncol);

SEXP compute_and_save_rvalues(SEXP xmat_, SEXP ymat_, SEXP zmat_,
    SEXP output_file_, SEXP names, SEXP rvalue_threshold, SEXP cores,
    SEXP with_return)
{
    int swap_y_and_z = ncols(zmat_) > ncols(ymat_);
    SEXP xmat = PROTECT(duplicate(xmat_));
    SEXP ymat = PROTECT(duplicate(swap_y_and_z ? zmat_ : ymat_));
    SEXP zmat = PROTECT(duplicate(swap_y_and_z ? ymat_ : zmat_));

    const char **xnames = extract_colnames(xmat_);
    const char **ynames = extract_colnames(ymat_);
    const char **znames = extract_colnames(zmat_);

    const char *output_file = (output_file_ == R_NilValue) ? NULL : CHAR(asChar(output_file_));

    SEXP result = compute_and_save(REAL(xmat), REAL(ymat), REAL(zmat),
        ncols(xmat), ncols(ymat), ncols(zmat), nrows(xmat), xnames,
        ynames, znames, names, output_file, asReal(rvalue_threshold),
        asInteger(cores), asInteger(with_return), swap_y_and_z);

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
    int xcol, int ycol, int zcol, int nrow, const char **xnames,
    const char **ynames, const char **znames, SEXP names,
    const char *filename, double rvalue_threshold, int cores,
    int with_return, int swap_y_and_z)
{
    F77_CALL(center)(xmat, &nrow, &xcol);
    F77_CALL(center)(ymat, &nrow, &ycol);
    F77_CALL(center)(zmat, &nrow, &zcol);

    int nelem = with_return ? xcol * ycol * zcol : xcol * ycol;
    const char **x = malloc(nelem * sizeof(*x));
    const char **y = malloc(nelem * sizeof(*y));
    const char **z = malloc(nelem * sizeof(*z));
    double      *r = malloc(nelem * sizeof(*r));
    int offset = 0;             // index of next element in x, y, z, r

    double *rvalue = malloc(xcol * ycol * sizeof(*rvalue));

    FILE *fp = NULL;
    if (filename) {
        fp = fopen(filename, "wb");
        print_header(fp, names);
    }
    for (int i = 0; i < zcol; i++) {
        F77_CALL(rval)(xmat, ymat, zmat + i * nrow, &xcol, &ycol, &nrow, rvalue, &cores);
        int nsignif = filter_and_annotate_rvalues(rvalue, xcol * ycol, rvalue_threshold,
            xnames, ynames, znames, x, y, z, r, xcol, i, swap_y_and_z, offset);
        if (filename)
            print_rvalues(fp, x, y, z, r, offset, nsignif);
        if (with_return)
            offset += nsignif;
        R_CheckUserInterrupt();
    }
    if (filename)
        fclose(fp);

    SEXP result = with_return ? create_data_frame(x, y, z, r, offset, names) : R_NilValue;

    free(x); free(y); free(z); free(r); free(rvalue);

    return result;
}

static void print_header(FILE *fp, SEXP names) {
    fprintf(fp, "%s\t%s\t%s\t%s\n",
        CHAR(STRING_ELT(names, 0)),
        CHAR(STRING_ELT(names, 1)),
        CHAR(STRING_ELT(names, 2)),
        CHAR(STRING_ELT(names, 3)));
}

static void print_rvalues(FILE *fp, const char **x, const char **y,
    const char **z, double *r, int offset, int nsignif)
{
    for (int i = offset; i < offset + nsignif; i++)
        fprintf(fp, "%s\t%s\t%s\t%.9f\n", x[i], y[i], z[i], r[i]);
}

static int filter_and_annotate_rvalues(double *rvalue, int n, double threshold,
    const char **xnames, const char **ynames, const char **znames,
    const char **x, const char **y, const char **z, double *r,
    int xcol, int zi, int swap_y_and_z, int initial_offset)
{
    int no_threshold = !R_FINITE(threshold);
    int offset = initial_offset;
    for (int i = 0; i < n; i++) {
        if (no_threshold || fabs(rvalue[i]) >= threshold) {
            x[offset] = xnames[i % xcol];
            y[offset] = ynames[swap_y_and_z ? zi : i / xcol];
            z[offset] = znames[swap_y_and_z ? i / xcol : zi];
            r[offset] = rvalue[i];
            ++offset;
        }
    }
    return offset - initial_offset;
}

static SEXP create_data_frame(const char **x_, const char **y_,
    const char **z_, double *r_, int nrow, SEXP names_)
{
    SEXP x = PROTECT(allocVector(STRSXP, nrow));
    SEXP y = PROTECT(allocVector(STRSXP, nrow));
    SEXP z = PROTECT(allocVector(STRSXP, nrow));
    SEXP r = PROTECT(allocVector(REALSXP, nrow));
    for (int i = 0; i < nrow; i++) {
        SET_STRING_ELT(x, i, mkChar(x_[i]));
        SET_STRING_ELT(y, i, mkChar(y_[i]));
        SET_STRING_ELT(z, i, mkChar(z_[i]));
        REAL(r)[i] = r_[i];
    }

    SEXP result = PROTECT(allocVector(VECSXP, 4));
    SET_VECTOR_ELT(result, 0, x);
    SET_VECTOR_ELT(result, 1, y);
    SET_VECTOR_ELT(result, 2, z);
    SET_VECTOR_ELT(result, 3, r);

    SEXP class = PROTECT(allocVector(STRSXP, 1));
    SET_STRING_ELT(class, 0, mkChar("data.frame"));
    setAttrib(result, R_ClassSymbol, class);

    SEXP names = PROTECT(duplicate(names_));
    setAttrib(result, R_NamesSymbol, names);

    SEXP rownames = PROTECT(allocVector(INTSXP, nrow));
    for (int i = 0; i < nrow; i++)
        INTEGER(rownames)[i] = i + 1;
    setAttrib(result, R_RowNamesSymbol, rownames);

    UNPROTECT(8);
    return result;
}
