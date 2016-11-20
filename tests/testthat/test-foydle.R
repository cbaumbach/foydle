source("setup.R")

context("foydle")

test_that("happy path", {
    mat <- create_data()
    actual <- foydle(mat$x, mat$y, mat$z)
    expected <- foydle_lm(mat$x, mat$y, mat$z)
    expect_same_contents(actual, expected)
})

test_that("we can specify the (column) names of the result table", {
    mat <- create_data()
    actual <- foydle(mat$x, mat$y, mat$z, "data/out.txt", names = c("a", "b", "c", "r-value"))
    expected <- foydle_lm(mat$x, mat$y, mat$z)
    names(expected) <- c("a", "b", "c", "r-value")
    expect_same_contents(actual, expected)
})

test_that("we can select results below a p-value threshold", {
    mat <- create_data()
    actual <- foydle(mat$x, mat$y, mat$z, "data/out.txt", pvalue_threshold = .7)
    d <- foydle_lm(mat$x, mat$y, mat$z)
    rvalue_threshold <- p2r(.7, df = nrow(mat$x) - 4)
    expected <- d[abs(d$r) >= rvalue_threshold, ]
    rownames(expected) <- NULL
    expect_same_contents(actual, expected)
})

test_that("we can use multiple cores", {
    mat <- create_data()
    actual <- foydle(mat$x, mat$y, mat$z, "data/out.txt", cores = 2)
    expect_same_contents(actual, foydle_lm(mat$x, mat$y, mat$z))
})

test_that("we can suppress the return value", {
    mat <- create_data()
    expect_null(foydle(mat$x, mat$y, mat$z, "data/out.txt", with_return = FALSE))
    actual <- read.delim("data/out.txt", as.is = TRUE)
    expected <- foydle_lm(mat$x, mat$y, mat$z)
    expect_same_contents(actual, expected)
})

context("conversions")

test_that("p2t finds the (1 - p/2)th quantile of a t-distribution", {
    p <- .1
    df <- 5
    expect_equal(p2t(p, df), stats::qt(1 - p/2, df))
})
