#' Fast P-Values For Linear Regression Interaction Term
#'
#' \code{foydle} computes p-values for the interaction term in the
#' linear regression
#'     \deqn{x = \beta_0 + \beta_1 y + \beta_2 z + \beta_3 y*z}{x = b0 + b1 y + b2 z + b3 y*z}
#' where \eqn{x}, \eqn{y}, and \eqn{z} are columns from \code{xmat},
#' \code{ymat}, and \code{zmat}, respectively.
#'
#' @param xmat,ymat,zmat Numeric matrices
#' @param output_file Path to output file
#' @param with_return If FALSE, return NULL (default: TRUE)
#' @param names Names for columns in output table
#' @param threshold Save only results with p-value below
#'     threshold (default: save all results)
#' @param cores Number of CPU cores to use (default 1)
#'
#' @return If \code{with_return} is \code{TRUE}, the function returns
#'     a data frame with columns "x", "y", "z", and "pvalue".  Columns
#'     "x", "y", and "z" contain the names of the columns of
#'     \code{xmat}, \code{ymat}, and \code{zmat} that were used in
#'     computing the p-value in the "pvalue" column.  The column names
#'     of the resulting table can be changed via the \code{names}
#'     argument.  The same data frame that is returned is also written
#'     to \code{output_file}.
#'
#' @export
foydle <- function(xmat, ymat, zmat, output_file = NULL, with_return = TRUE,
    names = c("x", "y", "z", "pvalue"), threshold = NULL, cores = 1)
{
    if (is.null(threshold))
        rvalue_threshold <- Inf
    else
        rvalue_threshold <- p2r(threshold, df = nrow(xmat) - 4)
    storage.mode(xmat) <- storage.mode(ymat) <- storage.mode(zmat) <- "double"
    storage <- create_storage()
    on.exit(free_storage(storage), add = TRUE)
    .Call("compute_and_save_rvalues", xmat, ymat, zmat, output_file, names,
        rvalue_threshold, as.integer(cores), as.logical(with_return), storage,
        PACKAGE = "foydle")
}

create_storage <- function() {
    .Call("create_storage", PACKAGE = "foydle")
}

free_storage <- function(storage) {
    .Call("free_storage", storage, PACKAGE = "foydle")
}

foydle_lm <- function(xmat, ymat, zmat) {
    find_pvalue <- function(xname, yname, zname) {
        x <- xmat[, xname]; y <- ymat[, yname]; z <- zmat[, zname]
        stats::coef(summary(stats::lm(x ~ y * z)))["y:z", "Pr(>|t|)"]
    }
    result <- expand.grid(x = colnames(xmat), y = colnames(ymat), z = colnames(zmat),
        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    result$pvalue <- as.double(Map(find_pvalue, result$x, result$y, result$z))
    result
}

t2r <- function(t, df) {
    t / sqrt(t^2 + df)
}

p2r <- function(p, df) {
    t2r(p2t(p, df), df)
}

p2t <- function(p, df) {
    # Find (1 - p/2)th quantile.
    stats::qt(log(p) - log(2), df, lower.tail = FALSE, log.p = TRUE)
}
