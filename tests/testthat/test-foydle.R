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
