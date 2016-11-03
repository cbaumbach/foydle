context("foydle")

test_that("foydle", {
    N <- 20
    xmat <- matrix(rnorm(N * 2), ncol = 2, dimnames = list(NULL, c("x1", "x2")))
    ymat <- matrix(rnorm(N * 3), ncol = 3, dimnames = list(NULL, c("y1", "y2", "y3")))
    zmat <- matrix(rnorm(N * 4), ncol = 4, dimnames = list(NULL, c("z1", "z2", "z3", "z4")))
    expect_equal(foydle(xmat, ymat, zmat), foydle_lm(xmat, ymat, zmat))
})
