# Hand-written counterpart to Rcpp::compileAttributes output, so the
# package builds in one shot without invoking compileAttributes.

rcpp_kappa_raw <- function(x, method, weight_type, values, em_options) {
  .Call(`_misskappa_rcpp_kappa_raw`, x, method, weight_type, values, em_options)
}

rcpp_kappa_continuous <- function(x, method, weight_type) {
  .Call(`_misskappa_rcpp_kappa_continuous`, x, method, weight_type)
}

rcpp_kappa_counts <- function(x, weight_type, values) {
  .Call(`_misskappa_rcpp_kappa_counts`, x, weight_type, values)
}
