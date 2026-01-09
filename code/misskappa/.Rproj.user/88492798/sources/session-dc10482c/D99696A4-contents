#' Unified Interface for Inter-Rater Agreement Coefficients (Categorical Data)
#'
#' @description
#' Computes various weighted agreement coefficients (Conger's Kappa, Fleiss' Kappa,
#' and Brennan-Prediger) for multiple raters on categorical data.
#' This is the primary interface for analyzing raw `(subjects x raters)` rating data
#' where ratings are discrete categories.
#'
#' For continuous data, see `kappa_continuous()`. For count data, see `kappa_counts()`.
#'
#' @details
#' The `method` argument dispatches to different underlying C++ estimation engines.
#'
#' The `"quadratic"` method is a special case that applies quadratic weights based on
#' numeric scores assigned to the categories.
#'
#' @param x A numeric matrix of ratings (subjects-by-raters). Must contain integer
#'   category labels. `NA`s are permitted.
#' @param method A string specifying the estimation method: "ml", "ipw", "available", or "quadratic".
#' @param weight A string specifying the weighting scheme for disagreements.
#'   Supported values are `"identity"`, `"unweighted"`, `"linear"`, `"quadratic"`,
#'   `"ordinal"`, `"radical"`, `"ratio"`, `"circular"`, and `"bipolar"`.
#'   This argument is ignored and set to `"quadratic"` if `method = "quadratic"`.
#' @param values An optional numeric vector of values (scores) assigned to each category.
#'   If `NULL`, categories are inferred from the data and assigned sequential integer values.
#' @param ... Additional options passed to the specific estimation method. For
#'   `method = "ml"`, this can include `em_options = list(tol = ..., max_iter = ...)`
#'   to control the EM algorithm.
#'
#' @return A list containing the point estimates and the variance-covariance matrix.
#'
#' @export
#' @seealso kappa_continuous, kappa_counts
kappa_raw <- function(x,
                      method = c("ml", "ipw", "available", "gwet", "quadratic"),
                      weight = c(
                        "identity", "unweighted", "linear", "quadratic",
                        "ordinal", "radical", "ratio", "circular", "bipolar"
                      ),
                      values = NULL,
                      ...) {
  method <- match.arg(method)
  weight_arg <- match.arg(weight)

  if (!is.matrix(x) && !is.data.frame(x)) stop("'x' must be a matrix or data frame.")
  x_mat <- as.matrix(x)
  if (!is.numeric(x_mat)) stop("'x' must be numeric.")

  x_int <- x_mat
  storage.mode(x_int) <- "integer"

  user_options <- list(...)
  default_options <- list(em_options = list(tol = 1e-8, max_iter = 10000))
  final_options <- utils::modifyList(default_options, user_options, keep.null = TRUE)

  results_cpp <- unified_kappa_raw_rcpp(
    x_r = x_int,
    method = method,
    weight_type = weight_arg,
    values = values,
    options = final_options
  )

  estimates <- c(results_cpp$estimates)
  vcov <- results_cpp$vcov
  names(estimates) <- c("Conger", "Fleiss", "Brennan-Prediger")
  dimnames(vcov) <- list(names(estimates), names(estimates))

  return(list(estimates = estimates, vcov = vcov))
}

#' Inter-Rater Agreement Coefficients for Continuous Data
#'
#' @description
#' Computes Conger's Kappa and Fleiss' Kappa for multiple raters on continuous data.
#' This is equivalent to Lin's Concordance Correlation Coefficient and its generalization.
#'
#' @details
#' This function provides a unified interface for agreement on continuous measurements.
#' The `method` argument allows choosing between a fast parametric (`"quadratic"`)
#' approach and a more robust semi-parametric (`"available"`, `"ipw"`) approach.
#' The `"quadratic"` method is equivalent to the intraclass correlation. The other methods
#' use a U-statistic-based estimator for the disagreement measures.
#'
#' This function does not compute the Brennan-Prediger coefficient, as its definition
#' of chance agreement is not generally applicable to continuous data.
#'
#' @param x A numeric matrix of continuous ratings (subjects-by-raters). `NA`s are permitted.
#' @param method A string specifying the estimation method. Supports
#'   `"quadratic"`, `"available"`, and `"ipw"`.
#' @param weight A string specifying the weighting scheme. For `"available"` and `"ipw"`,
#'   supported values are `"identity"`, `"linear"`, `"quadratic"`,
#'   `"radical"`, and `"ratio"`. This argument is ignored for `method = "quadratic"`.
#'
#' @return A list containing the point estimates for Conger's Kappa and Fleiss' Kappa,
#'   and their 2x2 variance-covariance matrix.
#'
#' @export
#' @seealso kappa_raw, kappa_counts
kappa_continuous <- function(x,
                             method = c("quadratic", "available", "ipw", "gwet"),
                             weight = c(
                               "quadratic", "linear",
                               "radical", "ratio"
                             )) {
  method <- match.arg(method)
  weight <- match.arg(weight)

  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a matrix or data frame.")
  }

  x_mat <- as.matrix(x)

  if (!is.numeric(x_mat)) {
    stop("'x' must be numeric.")
  }

  results_cpp <- unified_kappa_continuous_rcpp(
    x_r = x_mat,
    method = method,
    weight_type = weight
  )

  estimates <- c(results_cpp$estimates)
  vcov <- results_cpp$vcov
  names(estimates) <- c("Conger", "Fleiss")
  dimnames(vcov) <- list(names(estimates), names(estimates))

  return(list(estimates = estimates, vcov = vcov))
}


#' Unified Interface for Inter-Rater Agreement (Counts Data)
#'
#' @description
#' Computes Fleiss' Kappa and the Brennan-Prediger coefficient from count data.
#' For Conger's Kappa, data must be in the raw (subjects x raters) format.
#'
#' @param x A numeric matrix of counts (subjects x categories).
#' @param method A string specifying the estimation method: "ml", "quadratic", or "available".
#' @param weight A string specifying the weighting scheme. Ignored if `method = "quadratic"`.
#' @param values An optional numeric vector of category values.
#' @param r Optional number of raters in the data set. If NULL, it is inferred from
#' the maximum row sum of `x`.
#' @param ... Additional options, e.g., `em_options = list(...)` for method "ml".
#'
#' @return A list containing the point estimates and the variance-covariance matrix.
#'
#' @export
#' @seealso kappa_raw, kappa_continuous
kappa_counts <- function(x,
                         method = c("ml", "quadratic", "available"),
                         weight = c(
                           "identity", "unweighted", "linear", "quadratic",
                           "ordinal", "radical", "ratio", "circular", "bipolar"
                         ),
                         values = NULL,
                         r = NULL,
                         ...) {
  method <- match.arg(method)
  weight <- match.arg(weight)
  if (method == "available" && weight %in% c("unweighted")) weight <- "identity"

  if (!is.matrix(x) && !is.data.frame(x)) stop("'x' must be a matrix or data frame.")
  x_mat <- as.matrix(x)
  if (!is.numeric(x_mat) || any(x_mat < 0, na.rm = TRUE) || any(is.na(x_mat))) {
    stop("'x' must be a matrix of non-negative integer counts.")
  }
  x_int <- round(x_mat)
  storage.mode(x_int) <- "integer"

  if (is.null(r)) {
    r <- max(rowSums(x_int, na.rm = TRUE))
  }

  user_options <- list(...)
  default_options <- list(em_options = list(tol = 1e-8, max_iter = 10000))
  final_options <- utils::modifyList(default_options, user_options, keep.null = TRUE)

  if (method == "ipw") {
    stop("The 'ipw' method is not applicable to counts data, as individual rater information is lost. Use 'np' instead.")
  }

  results_cpp <- unified_kappa_counts_rcpp(
    x_r = x_int, r = as.integer(r), method = method,
    weight_type = weight, values = values, options = final_options
  )

  estimates <- c(results_cpp$estimates)
  vcov <- results_cpp$vcov

  names(estimates) <- c("Fleiss", "Brennan-Prediger")
  dimnames(vcov) <- list(names(estimates), names(estimates))

  return(list(estimates = estimates, vcov = vcov))
}
