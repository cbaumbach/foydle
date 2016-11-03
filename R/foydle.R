#' Fast Pearson correlation for linear regression interaction term
#'
#' \code{foydle} computes Pearson correlation coefficients for the
#' interaction term in the linear regression
#'
#'     \deqn{x = 1 + y + z + y*z}
#'
#' where \eqn{x}, \eqn{y}, and \eqn{z} are columns from \code{xmat},
#' \code{ymat}, and \code{zmat}, respectively.
#'
#' @param xmat,ymat,zmat Numeric matrices
#' @param output_file Path to output file
#'
#' @return None.  Instead a data frame with columns "x", "y", "z", and
#'     "r" is written to \code{output_file}.  Columns "x", "y", and
#'     "z" contain the names of the columns of \code{xmat},
#'     \code{ymat}, and \code{zmat} that were used in computing the
#'     corresponding Pearson correlation coefficient in the "r"
#'     column.
#'
#' @export
foydle <- function(xmat, ymat, zmat, output_file) {
    .Call("compute_and_save_rvalues",
        xmat = as.double(xmat),
        ymat = as.double(ymat),
        zmat = as.double(zmat),
        xcol = ncol(xmat),
        ycol = ncol(ymat),
        zcol = ncol(zmat),
        xnames = colnames(xmat),
        ynames = colnames(ymat),
        znames = colnames(zmat),
        n = nrow(xmat),
        output_file = output_file,
        PACKAGE = "foydle")
}

foydle_lm <- function(xmat, ymat, zmat) {
    find_rvalue <- function(xname, yname, zname) {
        x <- xmat[, xname]
        y <- ymat[, yname]
        z <- zmat[, zname]
        tvalue <- stats::coef(summary(stats::lm(x ~ y * z)))["y:z", "t value"]
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
