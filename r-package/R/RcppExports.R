# Hand-written counterpart to Rcpp::compileAttributes output, so the
# package builds in one shot without invoking compileAttributes.

rcpp_kappa_raw <- function(x, method, weight_type, values, em_options) {
  .Call(`_misskappa_rcpp_kappa_raw`, x, method, weight_type, values, em_options)
}

rcpp_kappa_continuous <- function(x, method, weight_type) {
  .Call(`_misskappa_rcpp_kappa_continuous`, x, method, weight_type)
}

rcpp_alpha_raw <- function(x, method, values, em_options) {
  .Call(`_misskappa_rcpp_alpha_raw`, x, method, values, em_options)
}

rcpp_kappa_counts <- function(x, weight_type, values) {
  .Call(`_misskappa_rcpp_kappa_counts`, x, weight_type, values)
}

rcpp_kappa_quadratic <- function(x, values, vcov_type, relative_kurtosis) {
  .Call(`_misskappa_rcpp_kappa_quadratic`, x, values, vcov_type, relative_kurtosis)
}

rcpp_kappa_quadratic_counts <- function(x, values, r_total) {
  .Call(`_misskappa_rcpp_kappa_quadratic_counts`, x, values, r_total)
}

rcpp_kappa_fiml_counts <- function(x, weight_type, values, r_total, em_options) {
  .Call(`_misskappa_rcpp_kappa_fiml_counts`, x, weight_type, values, r_total, em_options)
}

rcpp_fiml_louis_spectrum <- function(x, weight_type, values, em_options) {
  .Call(`_misskappa_rcpp_fiml_louis_spectrum`, x, weight_type, values, em_options)
}

rcpp_kappa_gwise_categorical <- function(x, distance_type, g, max_chance_tuples) {
  .Call(`_misskappa_rcpp_kappa_gwise_categorical`, x, distance_type, g, max_chance_tuples)
}

rcpp_kappa_gwise_continuous <- function(x, distance_type, g, max_chance_tuples) {
  .Call(`_misskappa_rcpp_kappa_gwise_continuous`, x, distance_type, g, max_chance_tuples)
}
