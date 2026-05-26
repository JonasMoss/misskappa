#' @description
#' Estimation and inference for weighted agreement coefficients, including
#' Fleiss' kappa, Conger's kappa (multirater Cohen's kappa), and the Brennan–Prediger coefficient.
#' Supports both rating-level and count-level data, with or without missing values.
#' Provides multiple estimators for incomplete data (available-case, inverse-probability weighting,
#' and maximum likelihood via an EM algorithm). Standard errors are computed using large-sample
#' theory (including Louis' method for the EM-based estimator).
#'
#' @importFrom Rcpp sourceCpp
#' @useDynLib misskappa, .registration = TRUE
#' @name misskappa-package
"_PACKAGE"
