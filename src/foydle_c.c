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
    int capacity, minimum_capacity;
} *Storage;

typedef struct ArgumentsStruct {
    double *xmat, *ymat, *zmat;
    const char *filename;
    SEXP names;
    double rvalue_threshold;
    int xcol, ycol, zcol, nrow;
    int swap_y_and_z, with_return, cores;
} *Arguments;

static SEXP compute_and_save(Arguments args, Storage storage);
static SEXP create_data_frame(Storage storage, int offset, SEXP names);
static const char **extract_colnames(SEXP mat);
static int filter_and_convert_rvalues(Storage storage, double threshold,
    int xcol, int zi, int swap_y_and_z, int initial_offset, int nrow);
static void initialize_storage(Storage storage, int xcol, int ycol, int zcol,
    int no_threshold, int with_return);
static void maybe_grow_storage(Storage storage, int size);
static void print_header(FILE *fp, SEXP names);
static void print_pvalues(FILE *fp, Storage storage, int offset, int nsignif);
static double r2p(double r, int df);

void F77_NAME(rval)(double *xmat, double *ymat, double *z, int *xcol, int *ycol, int *nrow, double *rvalue, int *cores);
void F77_NAME(center)(double *matrix, int *nrow, int *ncol);

SEXP compute_and_save_rvalues(SEXP xmat_, SEXP ymat_, SEXP zmat_,
    SEXP output_file_, SEXP names, SEXP rvalue_threshold, SEXP cores,
    SEXP with_return, SEXP storage_)
{
    struct ArgumentsStruct args;
    args.swap_y_and_z = ncols(zmat_) > ncols(ymat_);
    args.xmat = REAL(PROTECT(duplicate(xmat_)));
    args.ymat = REAL(PROTECT(duplicate(args.swap_y_and_z ? zmat_ : ymat_)));
    args.zmat = REAL(PROTECT(duplicate(args.swap_y_and_z ? ymat_ : zmat_)));
    args.xcol = ncols(xmat_);
    args.ycol = ncols(args.swap_y_and_z ? zmat_ : ymat_);
    args.zcol = ncols(args.swap_y_and_z ? ymat_ : zmat_);
    args.nrow = nrows(xmat_);
    args.filename = (output_file_ == R_NilValue) ? NULL : CHAR(asChar(output_file_));
    args.rvalue_threshold = asReal(rvalue_threshold);
    args.with_return = asInteger(with_return);
    args.cores = asInteger(cores);
    args.names = names;

    Storage storage = R_ExternalPtrAddr(storage_);
    storage->xnames = extract_colnames(xmat_);
    storage->ynames = extract_colnames(ymat_);
    storage->znames = extract_colnames(zmat_);

    SEXP result = compute_and_save(&args, storage);

    UNPROTECT(3);
    return result;
}

static SEXP compute_and_save(Arguments args, Storage storage)
{
    F77_CALL(center)(args->xmat, &args->nrow, &args->xcol);
    F77_CALL(center)(args->ymat, &args->nrow, &args->ycol);
    F77_CALL(center)(args->zmat, &args->nrow, &args->zcol);

    int no_threshold = !R_FINITE(args->rvalue_threshold);
    initialize_storage(storage, args->xcol, args->ycol, args->zcol, no_threshold, args->with_return);
    int offset = 0;             // index of next element in x, y, z, r

    FILE *fp = NULL;
    if (args->filename) {
        fp = fopen(args->filename, "wb");
        print_header(fp, args->names);
    }
    for (int i = 0; i < args->zcol; i++) {
        maybe_grow_storage(storage, offset);
        F77_CALL(rval)(args->xmat, args->ymat, args->zmat + i * args->nrow, &args->xcol, &args->ycol, &args->nrow, storage->rvalue, &args->cores);
        int nsignif = filter_and_convert_rvalues(storage, args->rvalue_threshold, args->xcol, i, args->swap_y_and_z, offset, args->nrow);
        if (args->filename && nsignif > 0)
            print_pvalues(fp, storage, offset, nsignif);
        if (args->with_return)
            offset += nsignif;
        R_CheckUserInterrupt();
    }
    if (args->filename)
        fclose(fp);

    return args->with_return ? create_data_frame(storage, offset, args->names) : R_NilValue;
}

static void initialize_storage(Storage storage, int xcol, int ycol, int zcol, int no_threshold, int with_return) {
    storage->minimum_capacity = xcol * ycol;
    int full_capacity = xcol * ycol * zcol;
    if (with_return && no_threshold)
        storage->capacity = full_capacity;
    else if (with_return)
        storage->capacity = 2 * storage->minimum_capacity;
    else
        storage->capacity = storage->minimum_capacity;
    storage->x = Calloc(storage->capacity, const char *);
    storage->y = Calloc(storage->capacity, const char *);
    storage->z = Calloc(storage->capacity, const char *);
    storage->pvalue = Calloc(storage->capacity, double);
    storage->rvalue = Calloc(storage->minimum_capacity, double);
}

static void maybe_grow_storage(Storage storage, int size) {
    int required_capacity = size + storage->minimum_capacity;
    while (storage->capacity < required_capacity) {
        storage->x = Realloc(storage->x, storage->capacity + storage->minimum_capacity, const char *);
        storage->y = Realloc(storage->y, storage->capacity + storage->minimum_capacity, const char *);
        storage->z = Realloc(storage->z, storage->capacity + storage->minimum_capacity, const char *);
        storage->pvalue = Realloc(storage->pvalue, storage->capacity + storage->minimum_capacity, double);
        storage->capacity += storage->minimum_capacity;
    }
}

static void print_header(FILE *fp, SEXP names) {
    fprintf(fp, "%s\t%s\t%s\t%s\n",
        CHAR(STRING_ELT(names, 0)),
        CHAR(STRING_ELT(names, 1)),
        CHAR(STRING_ELT(names, 2)),
        CHAR(STRING_ELT(names, 3)));
}

static void print_pvalues(FILE *fp, Storage storage, int offset, int nsignif)
{
    for (int i = offset; i < offset + nsignif; i++)
        fprintf(fp, "%s\t%s\t%s\t%.9f\n", storage->x[i], storage->y[i], storage->z[i], storage->pvalue[i]);
}

static const char **extract_colnames(SEXP mat) {
    SEXP colnames = VECTOR_ELT(getAttrib(mat, R_DimNamesSymbol), 1);
    int n = length(colnames);
    const char **result = Calloc(n, const char *);
    for (int i = 0; i < n; i++)
        result[i] = CHAR(STRING_ELT(colnames, i));
    return result;
}

static int filter_and_convert_rvalues(Storage storage, double threshold,
    int xcol, int zi, int swap_y_and_z, int initial_offset, int nrow)
{
    int df = nrow - 4;
    int no_threshold = !R_FINITE(threshold);
    int offset = initial_offset;
    for (int i = 0; i < storage->minimum_capacity; i++) {
        if (no_threshold || fabs(storage->rvalue[i]) >= threshold) {
            storage->x[offset] = storage->xnames[i % xcol];
            storage->y[offset] = storage->ynames[swap_y_and_z ? zi : i / xcol];
            storage->z[offset] = storage->znames[swap_y_and_z ? i / xcol : zi];
            storage->pvalue[offset] = r2p(storage->rvalue[i], df);
            ++offset;
        }
    }
    return offset - initial_offset;
}

static double r2p(double r, int df) {
    double t = fabs(r * sqrt(df / (1 - r*r)));
    return 2 * pt(t, df, 0, 0);
}

static SEXP create_data_frame(Storage storage, int nrow, SEXP names_)
{
    SEXP x = PROTECT(allocVector(STRSXP, nrow));
    SEXP y = PROTECT(allocVector(STRSXP, nrow));
    SEXP z = PROTECT(allocVector(STRSXP, nrow));
    SEXP pvalue = PROTECT(allocVector(REALSXP, nrow));
    for (int i = 0; i < nrow; i++) {
        SET_STRING_ELT(x, i, mkChar(storage->x[i]));
        SET_STRING_ELT(y, i, mkChar(storage->y[i]));
        SET_STRING_ELT(z, i, mkChar(storage->z[i]));
        REAL(pvalue)[i] = storage->pvalue[i];
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
