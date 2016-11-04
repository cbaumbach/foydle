source("setup.R")

context("foydle")

test_that("happy path", {
    mat <- create_data()
    result <- read.delim(as.is = TRUE,
        foydle(mat$x, mat$y, mat$z, output_file = "data/out.txt"))
    expect_equal(result, foydle_lm(mat$x, mat$y, mat$z))
})

test_that("we can specify the column names for the result table", {
    mat <- create_data()
    result <- read.delim(as.is = TRUE,
        foydle(mat$x, mat$y, mat$z, output_file = "data/out.txt",
            colnames = c("a", "b", "c")))
    expect_equal(colnames(result), c("a", "b", "c", "r"))
})

test_that("we can select results below a p-value threshold", {
    mat <- create_data()
    result <- read.delim(as.is = TRUE,
        foydle(mat$x, mat$y, mat$z, output_file = "data/out.txt",
            pvalue_threshold = .7))
    expected <- local({
        d <- foydle_lm(mat$x, mat$y, mat$z)
        d <- with(d, d[abs(r) >= p2r(.7, df = nrow(mat$x) - 4), ])
        rownames(d) <- NULL
        d
    })
    expect_equal(result, expected)
})

context("conversions")

test_that("p2t finds the (1 - p/2)th quantile of a t-distribution", {
    p <- .1; df <- 5
    expect_equal(p2t(p, df), stats::qt(1 - p/2, df))
})
