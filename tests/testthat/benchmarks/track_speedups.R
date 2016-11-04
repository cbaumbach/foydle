library(microbenchmark)

mat <- create_data(x = 20, y = 50, z = 100)

(data <- round(summary(microbenchmark(
    foydle(mat$x, mat$y, mat$z, "out.txt")))[-1]))

# min  lq mean median  uq max neval  comment
# 130 131  134    135 137 138   100  without optimization
