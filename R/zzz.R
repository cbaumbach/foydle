#' @useDynLib foydle
NULL

.onUnload <- function(libpath) {
    library.dynam.unload("foydle", libpath)
}
