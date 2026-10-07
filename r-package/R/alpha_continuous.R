# Coefficient alpha for continuous items under ignorable missingness, via a
# normal FIML / EM covariance and a sandwich delta-method standard error. The
# estimator lives in the C++ normal-FIML backend; the R code keeps the
# user-facing validation and object shape.
#
# This is the {continuous, missing} cell of the alpha-under-missingness grid:
# the categorical multinomial-EM route lives in alpha_cat_fiml() /
# rcpp_alpha_raw(). The fitted moments and sandwich are checked in the test
# suite against a pure-R reference EM (tests/testthat/helper-normal-fiml-oracle.R)
# and lavaan's saturated h1 estimator, and the standard error against
# coefficientalpha (varphi = 0).

.normal_fiml_from_cpp <- function(raw, estimate_names, method, weight) {
  estimates <- as.numeric(raw$estimates)
  names(estimates) <- estimate_names
  vcov_mat <- raw$vcov
  dimnames(vcov_mat) <- list(estimate_names, estimate_names)
  psi <- raw$psi
  if (prod(dim(psi)) > 0L) colnames(psi) <- estimate_names
  structure(
    list(
      estimates = estimates,
      vcov = vcov_mat,
      psi = psi,
      method = method,
      weight = weight,
      moments = list(
        mu = as.numeric(raw$mu),
        Sigma = raw$Sigma,
        iterations = as.integer(raw$iterations),
        converged = isTRUE(raw$converged)
      )
    ),
    class = "misskappa_estimate"
  )
}

# ---- public estimator -------------------------------------------------------

#' Coefficient alpha for continuous items under missing data (normal FIML)
#'
#' @description
#' Backend for `alpha(estimator = "nt_fiml")`. Estimates Cronbach's
#' coefficient alpha for a battery of continuous items with missing entries,
#' via a multivariate-normal covariance fitted by full-information maximum
#' likelihood (the EM algorithm), valid under ignorable (MCAR or MAR)
#' missingness. The default covariance model is saturated. `em_options =
#' list(covariance = "variance_components")` fits the compound-symmetry
#' (fixed-rater) covariance model \eqn{\Sigma = \psi J + e I}. It treats items
#' as fixed, distinct from ten Hove's random-rater MLE-RE / ICC. The standard error
#' is the sandwich delta-method contraction of the alpha gradient with the
#' asymptotic covariance of the fitted moments.
#'
#' The fitted moments reproduce lavaan's saturated estimator and the standard
#' error matches `coefficientalpha::alpha(varphi = 0)`.
#'
#' This backend is normally reached through [alpha()], which checks that every
#' item is observed and, for the default saturated covariance model, every item
#' pair is jointly observed before dispatch. Direct callers should enforce the
#' same fixed-item condition. For `covariance = "variance_components"` this is
#' a connected item co-observation graph.
#'
#' @param x A subjects-by-items numeric matrix or data frame. `NA` marks
#'   missing entries. Rows that are entirely missing are dropped.
#' @param em_options Named list tuning the EM fit: `tol` (convergence
#'   tolerance on the moments, default `1e-8`), `max_iter` (default `10000`),
#'   `covariance`, and `mean`. `covariance = "saturated"` fits the
#'   unrestricted normal covariance (default). `"variance_components"` fits
#'   the compound-symmetry (fixed-rater) covariance model
#'   \eqn{\Sigma = \psi J + e I}. `mean = "free"` estimates one mean per item
#'   (default). `mean = "equal"` estimates a common mean and is supported with
#'   `covariance = "variance_components"`. `fd_h` is accepted for backward
#'   compatibility and ignored because the observed information is analytic.
#'   Pass any subset.
#'
#' @return An object of class `misskappa_estimate` carrying one coefficient
#'   named `alpha` and its asymptotic covariance. Additional fields: `moments`
#'   (the fitted `mu`, `Sigma`, EM `iterations`, and `converged` flag) and
#'   `psi` (per-subject influence-function rows). Methods: `print`, `coef`,
#'   `vcov`, `confint`, `as.data.frame`, and `stats::influence`.
#'
#' @references
#' Zhang, Z., & Yuan, K.-H. (2016). Robust coefficients alpha and omega and
#' confidence intervals with outlying observations and missing data.
#' *Educational and Psychological Measurement*, 76(3), 387-411.
#'
#' @keywords internal
alpha_continuous <- function(x, em_options = list()) {
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a matrix or data frame.")
  }
  X <- as.matrix(x)
  if (!is.numeric(X)) stop("'x' must be numeric.")
  storage.mode(X) <- "double"
  p <- ncol(X)
  if (p < 2L) stop("coefficient alpha requires at least two items.")

  R <- is.finite(X)
  keep <- rowSums(R) > 0L
  if (!all(keep)) {
    X <- X[keep, , drop = FALSE]
    R <- R[keep, , drop = FALSE]
  }
  if (any(colSums(R) == 0L)) {
    stop("every item must be observed for at least one subject.")
  }

  opt <- utils::modifyList(
    list(tol = 1e-8, max_iter = 10000L), em_options
  )
  covariance <- .normal_covariance_model(opt)
  mean_model <- .normal_mean_model(opt)
  .check_normal_mean_covariance(mean_model, covariance)
  opt$covariance <- covariance
  opt$mean <- mean_model

  fit <- .normal_fiml_from_cpp(
    rcpp_alpha_normal_fiml(X, opt),
    estimate_names = "alpha",
    method = if (covariance == "variance_components" && mean_model == "equal") {
      "alpha-continuous-fiml-exchangeable"
    } else if (covariance == "variance_components") {
      "alpha-continuous-fiml-variance-components"
    } else {
      "alpha-continuous-fiml"
    },
    weight = "score"
  )
  fit$covariance <- covariance
  fit$mean <- mean_model
  fit
}
