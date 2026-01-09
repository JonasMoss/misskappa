#' @description
#' Estimation and inference for weighted agreement coefficients, including
#' Fleiss' kappa, Conger's kappa (multirater Cohen's kappa), and the Brennan–Prediger coefficient.
#' Supports both rating-level and count-level data, with or without missing values.
#' Missingness is assumed to be at random (MAR), and estimation is based on maximum likelihood
#' via the EM algorithm. Standard errors are computed using Louis' method, and the package
#' supports both asymptotic and bootstrap-based confidence intervals.
#'
#' @importFrom Rcpp sourceCpp
#' @useDynLib misskappa, .registration = TRUE
#' @name misskappa-package
"_PACKAGE"
