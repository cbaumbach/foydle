#' Fast Pearson Correlation For Linear Regression Interaction Term
#'
#' \code{foydle} computes Pearson correlation coefficients for the
#' interaction term in the linear regression
#'     \deqn{x = \beta_0 + \beta_1 y + \beta_2 z + \beta_3 y*z}{x = b0 + b1 y + b2 z + b3 y*z}
#' where \eqn{x}, \eqn{y}, and \eqn{z} are columns from \code{xmat},
#' \code{ymat}, and \code{zmat}, respectively.
#'
#' @param xmat,ymat,zmat Numeric matrices
#' @param output_file Path to output file
#' @param colnames Names for columns 1-3 in output table
#'
#' @return Returns \code{output_file} where a data frame with columns
#'     "x", "y", "z", and "r" was written.  Columns "x", "y", and "z"
#'     contain the names of the columns of \code{xmat}, \code{ymat},
#'     and \code{zmat} that were used in computing the corresponding
#'     Pearson correlation coefficient in the "r" column.  Column
#'     names other than "r" can be changed via the \code{colnames}
#'     argument.
#'
#' @export
foydle <- function(xmat, ymat, zmat, output_file, colnames = c("x", "y", "z")) {
    storage.mode(xmat) <- storage.mode(ymat) <- storage.mode(zmat) <- "double"
    .Call("compute_and_save_rvalues", xmat, ymat, zmat,
        nrow(xmat), output_file, colnames, PACKAGE = "foydle")
    output_file
}

foydle_lm <- function(xmat, ymat, zmat) {
    find_rvalue <- function(xname, yname, zname) {
        x <- xmat[, xname]
        y <- ymat[, yname]
        z <- zmat[, zname]
        model <- stats::lm(x ~ y * z)
        tvalue <- stats::coef(summary(model))["y:z", "t value"]
        t2r(tvalue, df = nrow(xmat) - 4L)
    }
    result <- expand.grid(
        x = colnames(xmat),
        y = colnames(ymat),
        z = colnames(zmat),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE)
    result$r <- as.double(Map(find_rvalue, result$x, result$y, result$z))
    result
}

t2r <- function(t, df) {
    t / sqrt(t^2 + df)
}

p2t <- function(p, df) {
    # Find (1 - p/2)th quantile.
    qt(log(p) - log(2), df, lower.tail = FALSE, log.p = TRUE)
}
