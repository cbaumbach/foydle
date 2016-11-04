create_data <- function(x = 2, y = 3, z = 4) {
    list(x = create_matrix(paste0("x", seq_len(x))),
        y = create_matrix(paste0("y", seq_len(y))),
        z = create_matrix(paste0("z", seq_len(z))))
}

create_matrix <- function(colnames, nobs = 20) {
    matrix(rnorm(nobs * length(colnames)),
        ncol = length(colnames),
        dimnames = list(NULL, colnames))
}
