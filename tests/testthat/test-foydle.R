source("setup.R")

context("foydle")

test_that("happy path", {
    mat <- create_data()
    foydle(mat$x, mat$y, mat$z, "data/out.txt")
    actual <- read.delim("data/out.txt", as.is = TRUE)
    expect_equal(actual, foydle_lm(mat$x, mat$y, mat$z))
})

test_that("we can specify the column names for the result table", {
    mat <- create_data()
    foydle(mat$x, mat$y, mat$z, "data/out.txt", colnames = c("a", "b", "c"))
    actual <- colnames(read.delim("data/out.txt", as.is = TRUE))
    expect_equal(actual, c("a", "b", "c", "r"))
})

test_that("we can select results below a p-value threshold", {
    mat <- create_data()

    foydle(mat$x, mat$y, mat$z, "data/out.txt", pvalue_threshold = .7)
    actual <- read.delim("data/out.txt", as.is = TRUE)

    d <- foydle_lm(mat$x, mat$y, mat$z)
    rvalue_threshold <- p2r(.7, df = nrow(mat$x) - 4)
    expected <- d[abs(d$r) >= rvalue_threshold, ]
    rownames(expected) <- NULL

    expect_equal(actual, expected)
})

context("conversions")

test_that("p2t finds the (1 - p/2)th quantile of a t-distribution", {
    p <- .1
    df <- 5
    expect_equal(p2t(p, df), stats::qt(1 - p/2, df))
})
