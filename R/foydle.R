#' @export
foydle <- function(xmat, ymat, zmat) {
}

foydle_lm <- function(xmat, ymat, zmat) {
    find_rvalue <- function(xname, yname, zname) {
        x <- xmat[, xname]
        y <- ymat[, yname]
        z <- zmat[, zname]
        tvalue <- coef(summary(lm(x ~ y * z)))["y:z", "t value"]
        t2r(tvalue, df = nrow(xmat) - 4L)
    }
    result <- expand.grid(
        x = colnames(xmat),
        y = colnames(ymat),
        z = colnames(zmat),
        stringsAsFactors = FALSE)
    result$r <- as.double(Map(find_rvalue, result$x, result$y, result$z))
    result
}

t2r <- function(t, df) {
    t / sqrt(t^2 + df)
}
