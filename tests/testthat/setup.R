unlink("data", recursive = TRUE, force = TRUE)
dir.create("data", showWarnings = FALSE)

expect_same_contents <- function(actual, expected) {
    pvalue_column <- names(expected)[ncol(expected)]
    if (nrow(actual) != nrow(expected)) {
        fail(sprintf("nrow(actual) != nrow(expected): %d vs %d",
            nrow(actual), nrow(expected)))
        return()
    }
    if (length(missing_columns <- setdiff(names(expected), names(actual))) > 0) {
        fail(sprintf("missing columns in actual: %s",
            paste(missing_columns, collapse = ", ")))
        return()
    }
    if (length(extra_columns <- setdiff(names(actual), names(expected))) > 0) {
        fail(sprintf("extra columns in actual: %s",
            paste(extra_columns, collapse = ", ")))
        return()
    }
    if (length(duplicated_columns <- unique(names(actual)[duplicated(names(actual))])) > 0) {
        fail(sprintf("duplicate column names in actual: %s",
            paste(duplicated_columns, collapse = ", ")))
        return()
    }
    merge_columns <- setdiff(names(expected), pvalue_column)
    xyz <- setNames(as.list(names(actual)[1:3]), c("x", "y", "z"))
    combine <- function(columns) {
        do.call(function(...) paste(..., sep = "\r"), columns)
    }
    if (length(strange_combo <- setdiff(combine(actual[merge_columns]), combine(expected[merge_columns]))) > 0) {
        value <- setNames(as.list(strsplit(strange_combo[1], "\r", fixed = TRUE)[[1]]), names(xyz))
        fail(sprintf("unexpected (%s, %s, %s) combination(s) in actual, first: (%s, %s, %s)",
            xyz$x, xyz$y, xyz$z, value$x, value$y, value$z))
        return()
    }
    merged <- merge(actual, expected, by = merge_columns, suffixes = c("_actual", "_expected"))
    is_close <- function(x, y, tolerance = 1e-5) {
        relative_differences <- abs((x - y) / ((x + y) / 2))
        vapply(relative_differences, `<`, logical(1), tolerance)
    }
    p_actual <- paste0(pvalue_column, "_actual")
    p_expected <- paste0(pvalue_column, "_expected")
    if (any(not_close <- ! is_close(merged[[p_actual]], merged[[p_expected]]))) {
        first <- which.max(not_close)
        fail(sprintf("actual and expected p-values differ, first at (%s, %s, %s) = (%s, %s, %s): %s vs %s",
            xyz$x, xyz$y, xyz$z, merged[[xyz$x]][first], merged[[xyz$y]][first],
            merged[[xyz$z]][first], merged[[p_actual]][first], merged[[p_expected]][first]))
        return()
    }
    succeed()
}
