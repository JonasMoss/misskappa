#' @keywords internal
#' @useDynLib misskappa, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats influence
"_PACKAGE"

.onUnload <- function(libpath) {
  library.dynam.unload("misskappa", libpath)
}
