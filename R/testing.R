create_data <- function(x = 2, y = 3, z = 4, nobs = 20) {
    list(x = create_matrix(paste0("x", seq_len(x)), nobs),
        y = create_matrix(paste0("y", seq_len(y)), nobs),
        z = create_matrix(paste0("z", seq_len(z)), nobs))
}

create_matrix <- function(colnames, nobs) {
    matrix(stats::rnorm(nobs * length(colnames)),
        ncol = length(colnames),
        dimnames = list(NULL, colnames))
}
