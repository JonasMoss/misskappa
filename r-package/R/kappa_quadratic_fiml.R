# Quadratic (Conger / Fleiss / Brennan-Prediger) kappa for continuous ratings
# under ignorable missingness, via a normal FIML / EM covariance and a
# delta-method standard error. This is the quadratic-kappa counterpart of
# alpha_continuous(): the two share the C++ normal EM, the analytic casewise
# scores, the observed information, and the sandwich contraction. The only
# coefficient-specific piece is the gradient (a pure-R reference copy lives in
# tests/testthat/helper-normal-fiml-oracle.R).
#
# Unlike alpha, the quadratic kappas depend on the rater means through
# t3 = sum_j (mu_j - mubar)^2, so the gradient carries a non-zero mean block.

.normal_covariance_model <- function(em_options) {
  cov <- em_options$covariance
  if (is.null(cov)) return("saturated")
  if (!is.character(cov) || length(cov) != 1L || is.na(cov)) {
    stop("em_options$covariance must be a single string.", call. = FALSE)
  }
  cov <- gsub("-", "_", tolower(trimws(cov)), fixed = TRUE)
  if (cov %in% c("saturated")) return("saturated")
  if (cov %in% c("variance_components")) return("variance_components")
  stop("unsupported em_options$covariance: ", em_options$covariance,
       call. = FALSE)
}

.normal_mean_model <- function(em_options) {
  mean <- em_options$mean
  if (is.null(mean)) return("free")
  if (!is.character(mean) || length(mean) != 1L || is.na(mean)) {
    stop("em_options$mean must be a single string.", call. = FALSE)
  }
  mean <- gsub("-", "_", tolower(trimws(mean)), fixed = TRUE)
  if (mean %in% c("free")) return("free")
  if (mean %in% c("equal", "common")) return("equal")
  stop("unsupported em_options$mean: ", em_options$mean, call. = FALSE)
}

.check_normal_mean_covariance <- function(mean, covariance) {
  if (mean == "equal" && covariance != "variance_components") {
    stop("em_options$mean = \"equal\" currently requires ",
         "em_options$covariance = \"variance_components\".",
         call. = FALSE)
  }
  invisible(NULL)
}

.quadratic_bp_chance <- function(x, values = NULL) {
  if (is.null(values)) {
    observed <- x[is.finite(x)]
    if (!length(observed)) stop("at least one rating must be observed.", call. = FALSE)
    rng <- range(observed)
    chance <- diff(rng)^2 / 6
  } else {
    if (!is.numeric(values)) stop("'values' must be numeric.", call. = FALSE)
    values <- as.numeric(values)
    if (length(values) < 2L || any(!is.finite(values))) {
      stop("'values' must contain at least two finite scores.", call. = FALSE)
    }
    C <- length(values)
    chance <- (2 / C^2) * (C * sum(values^2) - sum(values)^2)
  }
  if (!is.finite(chance) || chance <= 0) {
    stop("Brennan-Prediger's quadratic chance disagreement is zero.",
         call. = FALSE)
  }
  chance
}

# ---- public estimator -------------------------------------------------------

#' Quadratic kappa under missing data (normal FIML)
#'
#' @description
#' Backend for `kappa(estimator = "nt_fiml")`. Estimates the quadratically
#' weighted Conger, Fleiss, and Brennan-Prediger agreement coefficients for
#' continuous (or numerically scored) ratings with missing entries, via a
#' multivariate-normal covariance fitted by full-information maximum likelihood
#' (the EM algorithm), valid under ignorable (MCAR or MAR) missingness. The
#' default covariance model is saturated. `em_options = list(covariance =
#' "variance_components")` fits the compound-symmetry (fixed-rater) covariance
#' model \eqn{\Sigma = \psi J + e I}. It treats raters as fixed, distinct from
#' ten Hove's random-rater MLE-RE / ICC. The coefficients are smooth functions of
#' the fitted mean and covariance, and the standard error is a delta-method
#' contraction of their gradient with the asymptotic covariance of the fitted
#' moments. This is the quadratic-kappa counterpart of [alpha_continuous()] and
#' shares its EM, casewise-score, and sandwich machinery.
#'
#' This backend is normally reached through [kappa()], which checks that every
#' rater is observed and, for the default saturated covariance model, every
#' rater pair is jointly observed before dispatch. Direct callers should
#' enforce the same fixed-rater condition. For `covariance =
#' "variance_components"` this is a connected rater co-observation graph.
#'
#' @param x A subjects-by-raters numeric matrix or data frame. `NA` marks
#'   missing entries. Rows that are entirely missing are dropped.
#' @param values Optional finite numeric score vector defining the uniform
#'   Brennan-Prediger chance baseline. If `NULL`, the continuous uniform
#'   baseline on the observed score range is used.
#' @param em_options Named list tuning the EM fit: `tol` (default `1e-8`),
#'   `max_iter` (default `10000`), `covariance`, and `mean`. `covariance =
#'   "saturated"` fits the unrestricted normal FIML model.
#'   `"variance_components"` fits the compound-symmetry covariance model
#'   \eqn{\Sigma = \psi J + e I} with \eqn{\psi \ge 0}. `mean = "free"`
#'   estimates one mean per rater/item (default). `mean = "equal"` estimates a
#'   common mean and is supported with `covariance = "variance_components"`,
#'   giving the exchangeable normal model used by counts-format normal FIML.
#'   Pass any subset.
#'
#' @return An object of class `misskappa_estimate` carrying the `Conger`,
#'   `Fleiss`, and `Brennan-Prediger` coefficients and their 3x3 asymptotic
#'   covariance. Additional fields: `moments` (the fitted `mu`, `Sigma`, EM
#'   `iterations`, and `converged` flag) and `psi` (per-subject
#'   influence-function rows). Methods: `print`, `coef`, `vcov`, `confint`,
#'   `as.data.frame`, and `stats::influence`.
#'
#' @noRd
kappa_quadratic_fiml <- function(x, values = NULL, em_options = list()) {
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a matrix or data frame.")
  }
  X <- as.matrix(x)
  if (!is.numeric(X)) stop("'x' must be numeric.")
  storage.mode(X) <- "double"
  p <- ncol(X)
  if (p < 2L) stop("quadratic kappa requires at least two raters.")

  R <- is.finite(X)
  keep <- rowSums(R) > 0L
  if (!all(keep)) {
    X <- X[keep, , drop = FALSE]
    R <- R[keep, , drop = FALSE]
  }
  if (any(colSums(R) == 0L)) {
    stop("every rater must be observed for at least one subject.")
  }

  opt <- utils::modifyList(
    list(tol = 1e-8, max_iter = 10000L), em_options
  )
  covariance <- .normal_covariance_model(opt)
  mean_model <- .normal_mean_model(opt)
  .check_normal_mean_covariance(mean_model, covariance)
  opt$covariance <- covariance
  opt$mean <- mean_model
  bp_chance <- .quadratic_bp_chance(X, values)

  fit <- .normal_fiml_from_cpp(
    rcpp_kappa_quadratic_fiml(X, values, opt),
    estimate_names = c("Conger", "Fleiss", "Brennan-Prediger"),
    method = if (covariance == "variance_components" && mean_model == "equal") {
      "quadratic-fiml-exchangeable"
    } else if (covariance == "variance_components") {
      "quadratic-fiml-variance-components"
    } else {
      "quadratic-fiml"
    },
    weight = "quadratic"
  )
  fit$covariance <- covariance
  fit$mean <- mean_model
  fit$bp_chance <- bp_chance
  fit
}
