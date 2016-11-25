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
    int stored, previously_stored;
} *Storage;

typedef struct ArgumentsStruct {
    double *xmat, *ymat, *zmat;
    const char *filename;
    SEXP xnames, ynames, znames, names;
    double rvalue_threshold;
    int xcol, ycol, zcol, nrow;
    int swap_y_and_z, with_return, cores;
} *Arguments;

static void allocate_memory_for_saved_results(Storage storage, Arguments args);
static SEXP compute_and_save(Arguments args, Storage storage);
static void compute_results(Storage storage, Arguments args, int zcolumn);
static SEXP create_data_frame(Storage storage, SEXP names);
static const char **convert_column_names(SEXP colnames);
static void convert_column_names_to_string_arrays(Storage storage, Arguments args);
static void store_results(Storage storage, Arguments args, int zindex);
static void initialize_storage(Storage storage, Arguments args);
static void maybe_grow_storage(Storage storage);
static void print_header(FILE *fp, SEXP names);
static void print_results(FILE *fp, Storage storage);
static double r2p(double r, int df);
static void reset_storage(Storage storage);

void F77_NAME(rval)(double *xmat, double *ymat, double *z, int *xcol, int *ycol, int *nrow, double *rvalue, int *cores);
void F77_NAME(center)(double *matrix, int *nrow, int *ncol);

SEXP compute_and_save_rvalues(SEXP xmat, SEXP ymat, SEXP zmat,
    SEXP output_file, SEXP names, SEXP rvalue_threshold, SEXP cores,
    SEXP with_return, SEXP storage_)
{
    struct ArgumentsStruct args;
    args.swap_y_and_z = ncols(zmat) > ncols(ymat);
    args.xmat = REAL(PROTECT(duplicate(xmat)));
    args.ymat = REAL(PROTECT(duplicate(args.swap_y_and_z ? zmat : ymat)));
    args.zmat = REAL(PROTECT(duplicate(args.swap_y_and_z ? ymat : zmat)));
    args.xnames = VECTOR_ELT(getAttrib(xmat, R_DimNamesSymbol), 1);
    args.ynames = VECTOR_ELT(getAttrib(ymat, R_DimNamesSymbol), 1);
    args.znames = VECTOR_ELT(getAttrib(zmat, R_DimNamesSymbol), 1);
    args.xcol = ncols(xmat);
    args.ycol = ncols(args.swap_y_and_z ? zmat : ymat);
    args.zcol = ncols(args.swap_y_and_z ? ymat : zmat);
    args.nrow = nrows(xmat);
    args.filename = (output_file == R_NilValue) ? NULL : CHAR(asChar(output_file));
    args.rvalue_threshold = asReal(rvalue_threshold);
    args.with_return = asInteger(with_return);
    args.cores = asInteger(cores);
    args.names = names;

    Storage storage = R_ExternalPtrAddr(storage_);
    initialize_storage(storage, &args);

    SEXP result = compute_and_save(&args, storage);

    UNPROTECT(3);
    return result;
}

static SEXP compute_and_save(Arguments args, Storage storage)
{
    FILE *fp = NULL;
    if (args->filename) {
        fp = fopen(args->filename, "wb");
        print_header(fp, args->names);
    }

    F77_CALL(center)(args->xmat, &args->nrow, &args->xcol);
    F77_CALL(center)(args->ymat, &args->nrow, &args->ycol);
    F77_CALL(center)(args->zmat, &args->nrow, &args->zcol);

    for (int i = 0; i < args->zcol; i++) {
        maybe_grow_storage(storage);
        compute_results(storage, args, i);
        store_results(storage, args, i);
        if (args->filename && storage->stored - storage->previously_stored > 0)
            print_results(fp, storage);
        if (!args->with_return)
            reset_storage(storage);
        R_CheckUserInterrupt();
    }

    if (args->filename)
        fclose(fp);

    return args->with_return ? create_data_frame(storage, args->names) : R_NilValue;
}

static void compute_results(Storage storage, Arguments args, int zcolumn) {
    F77_CALL(rval)(args->xmat, args->ymat, args->zmat + zcolumn * args->nrow,
        &args->xcol, &args->ycol, &args->nrow, storage->rvalue, &args->cores);
}

static void initialize_storage(Storage storage, Arguments args) {
    convert_column_names_to_string_arrays(storage, args);
    allocate_memory_for_saved_results(storage, args);
}

static void convert_column_names_to_string_arrays(Storage storage, Arguments args) {
    storage->xnames = convert_column_names(args->xnames);
    storage->ynames = convert_column_names(args->ynames);
    storage->znames = convert_column_names(args->znames);
}

static const char **convert_column_names(SEXP colnames) {
    int n = length(colnames);
    const char **result = Calloc(n, const char *);
    for (int i = 0; i < n; i++)
        result[i] = CHAR(STRING_ELT(colnames, i));
    return result;
}

static void allocate_memory_for_saved_results(Storage storage, Arguments args) {
    int full_capacity = args->xcol * args->ycol * args->zcol;
    int minimum_capacity = args->xcol * args->ycol;
    int no_threshold = !R_FINITE(args->rvalue_threshold);
    if (args->with_return && no_threshold)
        storage->capacity = full_capacity;
    else if (args->with_return)
        storage->capacity = 2 * minimum_capacity;
    else
        storage->capacity = minimum_capacity;
    storage->minimum_capacity = minimum_capacity;
    storage->x = Calloc(storage->capacity, const char *);
    storage->y = Calloc(storage->capacity, const char *);
    storage->z = Calloc(storage->capacity, const char *);
    storage->pvalue = Calloc(storage->capacity, double);
    storage->rvalue = Calloc(minimum_capacity, double);
}

static void maybe_grow_storage(Storage storage) {
    int required_capacity = storage->stored + storage->minimum_capacity;
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

static void print_results(FILE *fp, Storage storage) {
    for (int i = storage->previously_stored; i < storage->stored; i++)
        fprintf(fp, "%s\t%s\t%s\t%.9f\n", storage->x[i], storage->y[i], storage->z[i], storage->pvalue[i]);
}

static void store_results(Storage storage, Arguments args, int zindex) {
    int df = args->nrow - 4;
    int no_threshold = !R_FINITE(args->rvalue_threshold);
    storage->previously_stored = storage->stored;
    for (int i = 0; i < storage->minimum_capacity; i++) {
        if (no_threshold || fabs(storage->rvalue[i]) >= args->rvalue_threshold) {
            storage->x[storage->stored] = storage->xnames[i % args->xcol];
            storage->y[storage->stored] = storage->ynames[args->swap_y_and_z ? zindex : i / args->xcol];
            storage->z[storage->stored] = storage->znames[args->swap_y_and_z ? i / args->xcol : zindex];
            storage->pvalue[storage->stored] = r2p(storage->rvalue[i], df);
            ++storage->stored;
        }
    }
}

static void reset_storage(Storage storage) {
    storage->previously_stored = 0;
    storage->stored = 0;
}

static double r2p(double r, int df) {
    double t = fabs(r * sqrt(df / (1 - r*r)));
    return 2 * pt(t, df, 0, 0);
}

static SEXP create_data_frame(Storage storage, SEXP names_) {
    SEXP x = PROTECT(allocVector(STRSXP, storage->stored));
    SEXP y = PROTECT(allocVector(STRSXP, storage->stored));
    SEXP z = PROTECT(allocVector(STRSXP, storage->stored));
    SEXP pvalue = PROTECT(allocVector(REALSXP, storage->stored));
    for (int i = 0; i < storage->stored; i++) {
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

    SEXP rownames = PROTECT(allocVector(INTSXP, storage->stored));
    for (int i = 0; i < storage->stored; i++)
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
    storage->capacity = 0;
    storage->minimum_capacity = 0;
    storage->stored = 0;
    storage->previously_stored = 0;
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
