library(microbenchmark)

mat <- create_data(x = 20, y = 50, z = 100)

(data <- round(summary(microbenchmark(
    foydle(mat$x, mat$y, mat$z, threshold = .05 / 20 / 50 / 100)))[-1]))

# min  lq mean median  uq max neval  comment
# 130 131  134    135 137 138   100  without optimization
#  89  90   91     91  93  94   100  move computations out of inner loop
#  85  86   87     86  87  97   100  precompute YY and ZZ
#  83  85   86     85  87  88   100  center all matrix columns ahead of time
#  80  81   81     81  81  86   100  swap ymat and zmat if ncol(zmat) > ncol(ymat)
#  71  71   71     71  71  74   100  precompute norms of X and YZ
# ---------------------------------------------------------------------------------
#  70  70   71     71  72  76   100  compile FORTRAN code with -O3 instead of -O2
#  69  70   70     70  70  74   100  also compile C code with -O3 instead of -O2
# ---------------------------------------------------------------------------------
# 150 151  155    152 154 274   100  with output file and return value
#  91  92   94     92  93 211   100  without output file
# 134 134  135    134 135 144   100  with output file and without return value
#  22  23   23     23  23  24   100  without output file and threshold .05
#  17  17   17     17  17  18   100  without output file and threshold .05 / 20 / 50 / 100
