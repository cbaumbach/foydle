#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <Rmath.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

typedef struct StorageStruct {
    const char **xnames, **ynames, **znames;
    const char **x, **y, **z;
    double *pvalue, *rvalue;
} *Storage;

static SEXP compute_and_save(double *xmat, double *ymat, double *zmat,
    int xcol, int ycol, int zcol, int nrow, const char **xnames,
    const char **ynames, const char **znames, SEXP names,
    const char *filename, double rvalue_threshold, int cores,
    int with_return, int swap_y_and_z, Storage storage);
static SEXP create_data_frame(const char **x, const char **y,
    const char **z, double *pvalue, int offset, SEXP names);
static const char **extract_colnames(SEXP mat);
static int filter_and_convert_rvalues(double *rvalue, int n, double threshold,
    const char **xnames, const char **ynames, const char **znames,
    const char **x, const char **y, const char **z, double *pvalue,
    int xcol, int zi, int swap_y_and_z, int initial_offset, int nrow);
static void print_header(FILE *fp, SEXP names);
static void print_rvalues(FILE *fp, const char **x, const char **y,
    const char **z, double *pvalue, int offset, int nsignif);
static double r2p(double r, int df);

void F77_NAME(rval)(double *xmat, double *ymat, double *z, int *xcol, int *ycol, int *nrow, double *rvalue, int *cores);
void F77_NAME(center)(double *matrix, int *nrow, int *ncol);

SEXP compute_and_save_rvalues(SEXP xmat_, SEXP ymat_, SEXP zmat_,
    SEXP output_file_, SEXP names, SEXP rvalue_threshold, SEXP cores,
    SEXP with_return, SEXP storage_)
{
    int swap_y_and_z = ncols(zmat_) > ncols(ymat_);
    SEXP xmat = PROTECT(duplicate(xmat_));
    SEXP ymat = PROTECT(duplicate(swap_y_and_z ? zmat_ : ymat_));
    SEXP zmat = PROTECT(duplicate(swap_y_and_z ? ymat_ : zmat_));

    Storage storage = R_ExternalPtrAddr(storage_);
    const char **xnames = storage->xnames = extract_colnames(xmat_);
    const char **ynames = storage->ynames = extract_colnames(ymat_);
    const char **znames = storage->znames = extract_colnames(zmat_);

    const char *output_file = (output_file_ == R_NilValue) ? NULL : CHAR(asChar(output_file_));

    SEXP result = compute_and_save(REAL(xmat), REAL(ymat), REAL(zmat),
        ncols(xmat), ncols(ymat), ncols(zmat), nrows(xmat), xnames,
        ynames, znames, names, output_file, asReal(rvalue_threshold),
        asInteger(cores), asInteger(with_return), swap_y_and_z, storage);

    UNPROTECT(3);
    return result;
}

static const char **extract_colnames(SEXP mat) {
    SEXP colnames = VECTOR_ELT(getAttrib(mat, R_DimNamesSymbol), 1);
    int n = length(colnames);
    const char **result = Calloc(n, const char *);
    for (int i = 0; i < n; i++)
        result[i] = CHAR(STRING_ELT(colnames, i));
    return result;
}

static SEXP compute_and_save(double *xmat, double *ymat, double *zmat,
    int xcol, int ycol, int zcol, int nrow, const char **xnames,
    const char **ynames, const char **znames, SEXP names,
    const char *filename, double rvalue_threshold, int cores,
    int with_return, int swap_y_and_z, Storage storage)
{
    F77_CALL(center)(xmat, &nrow, &xcol);
    F77_CALL(center)(ymat, &nrow, &ycol);
    F77_CALL(center)(zmat, &nrow, &zcol);

    int minimum_capacity = xcol * ycol;
    int full_capacity = xcol * ycol * zcol;
    int capacity;
    int no_threshold = !R_FINITE(rvalue_threshold);
    if (with_return && no_threshold)
        capacity = full_capacity;
    else if (with_return)
        capacity = 2 * minimum_capacity;
    else
        capacity = minimum_capacity;

    const char **x = storage->x = Calloc(capacity, const char *);
    const char **y = storage->y = Calloc(capacity, const char *);
    const char **z = storage->z = Calloc(capacity, const char *);
    double *pvalue = storage->pvalue = Calloc(capacity, double);
    int offset = 0;             // index of next element in x, y, z, r
    double *rvalue = storage->rvalue = Calloc(minimum_capacity, double);

    FILE *fp = NULL;
    if (filename) {
        fp = fopen(filename, "wb");
        print_header(fp, names);
    }
    for (int i = 0; i < zcol; i++) {
        int required_capacity = offset + minimum_capacity;
        while (capacity < required_capacity) {
            x = storage->x = Realloc(x, capacity + minimum_capacity, const char *);
            y = storage->y = Realloc(y, capacity + minimum_capacity, const char *);
            z = storage->z = Realloc(z, capacity + minimum_capacity, const char *);
            pvalue = storage->pvalue = Realloc(pvalue, capacity + minimum_capacity, double);
            capacity += minimum_capacity;
        }
        F77_CALL(rval)(xmat, ymat, zmat + i * nrow, &xcol, &ycol, &nrow, rvalue, &cores);
        int nsignif = filter_and_convert_rvalues(rvalue, xcol * ycol, rvalue_threshold,
            xnames, ynames, znames, x, y, z, pvalue, xcol, i, swap_y_and_z, offset, nrow);
        if (filename)
            print_rvalues(fp, x, y, z, pvalue, offset, nsignif);
        if (with_return)
            offset += nsignif;
        R_CheckUserInterrupt();
    }
    if (filename)
        fclose(fp);

    return with_return ? create_data_frame(x, y, z, pvalue, offset, names) : R_NilValue;
}

static void print_header(FILE *fp, SEXP names) {
    fprintf(fp, "%s\t%s\t%s\t%s\n",
        CHAR(STRING_ELT(names, 0)),
        CHAR(STRING_ELT(names, 1)),
        CHAR(STRING_ELT(names, 2)),
        CHAR(STRING_ELT(names, 3)));
}

static void print_rvalues(FILE *fp, const char **x, const char **y,
    const char **z, double *pvalue, int offset, int nsignif)
{
    for (int i = offset; i < offset + nsignif; i++)
        fprintf(fp, "%s\t%s\t%s\t%.9f\n", x[i], y[i], z[i], pvalue[i]);
}

static int filter_and_convert_rvalues(double *rvalue, int n, double threshold,
    const char **xnames, const char **ynames, const char **znames,
    const char **x, const char **y, const char **z, double *pvalue,
    int xcol, int zi, int swap_y_and_z, int initial_offset, int nrow)
{
    int df = nrow - 4;
    int no_threshold = !R_FINITE(threshold);
    int offset = initial_offset;
    for (int i = 0; i < n; i++) {
        if (no_threshold || fabs(rvalue[i]) >= threshold) {
            x[offset] = xnames[i % xcol];
            y[offset] = ynames[swap_y_and_z ? zi : i / xcol];
            z[offset] = znames[swap_y_and_z ? i / xcol : zi];
            pvalue[offset] = r2p(rvalue[i], df);
            ++offset;
        }
    }
    return offset - initial_offset;
}

static double r2p(double r, int df) {
    double t = fabs(r * sqrt(df / (1 - r*r)));
    return 2 * pt(t, df, 0, 0);
}

static SEXP create_data_frame(const char **x_, const char **y_,
    const char **z_, double *pvalue_, int nrow, SEXP names_)
{
    SEXP x = PROTECT(allocVector(STRSXP, nrow));
    SEXP y = PROTECT(allocVector(STRSXP, nrow));
    SEXP z = PROTECT(allocVector(STRSXP, nrow));
    SEXP pvalue = PROTECT(allocVector(REALSXP, nrow));
    for (int i = 0; i < nrow; i++) {
        SET_STRING_ELT(x, i, mkChar(x_[i]));
        SET_STRING_ELT(y, i, mkChar(y_[i]));
        SET_STRING_ELT(z, i, mkChar(z_[i]));
        REAL(pvalue)[i] = pvalue_[i];
    }

    SEXP result = PROTECT(allocVector(VECSXP, 4));
    SET_VECTOR_ELT(result, 0, x);
    SET_VECTOR_ELT(result, 1, y);
    SET_VECTOR_ELT(result, 2, z);
    SET_VECTOR_ELT(result, 3, pvalue);

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

SEXP create_storage(void) {
    Storage storage = Calloc(1, struct StorageStruct);
    storage->xnames = NULL;
    storage->ynames = NULL;
    storage->znames = NULL;
    storage->x = NULL;
    storage->y = NULL;
    storage->z = NULL;
    storage->pvalue = NULL;
    storage->rvalue = NULL;
    return R_MakeExternalPtr(storage, R_NilValue, R_NilValue);
}

SEXP free_storage(SEXP storage_) {
    Storage storage = R_ExternalPtrAddr(storage_);
    if (!storage)
        return R_NilValue;
    if (storage->xnames)
        Free(storage->xnames);
    if (storage->ynames)
        Free(storage->ynames);
    if (storage->znames)
        Free(storage->znames);
    if (storage->x)
        Free(storage->x);
    if (storage->y)
        Free(storage->y);
    if (storage->z)
        Free(storage->z);
    if (storage->pvalue)
        Free(storage->pvalue);
    if (storage->rvalue)
        Free(storage->rvalue);
    R_ClearExternalPtr(storage_);
    return R_NilValue;
}
