#' @keywords internal
#' @useDynLib misskappa, .registration = TRUE
#' @importFrom Rcpp sourceCpp
"_PACKAGE"

.onUnload <- function(libpath) {
  library.dynam.unload("misskappa", libpath)
}
