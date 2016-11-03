unlink("data", recursive = TRUE, force = TRUE)
dir.create("data", showWarnings = FALSE)

context("foydle")

test_that("foydle", {
    N <- 20
    xmat <- matrix(rnorm(N * 2), ncol = 2, dimnames = list(NULL, c("x1", "x2")))
    ymat <- matrix(rnorm(N * 3), ncol = 3, dimnames = list(NULL, c("y1", "y2", "y3")))
    zmat <- matrix(rnorm(N * 4), ncol = 4, dimnames = list(NULL, c("z1", "z2", "z3", "z4")))

    result <- read.delim(foydle(xmat, ymat, zmat, output_file = "data/out.txt"), as.is = TRUE)

    expect_equal(result, foydle_lm(xmat, ymat, zmat))
})
