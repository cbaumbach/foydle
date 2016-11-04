library(microbenchmark)

mat <- create_data(x = 20, y = 50, z = 100)

(data <- round(summary(microbenchmark(
    foydle(mat$x, mat$y, mat$z, "out.txt")))[-1]))

# min  lq mean median  uq max neval  comment
# 130 131  134    135 137 138   100  without optimization
#  89  90   91     91  93  94   100  move computations out of inner loop
#  85  86   87     86  87  97   100  precompute YY and ZZ
#  83  85   86     85  87  88   100  center all matrix columns ahead of time
