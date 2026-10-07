# Build a numeric item-score matrix for alpha(). With values = NULL, the input
# is already the score matrix. With values supplied, observed categories are
# sorted and mapped to the supplied scores, preserving the old alpha() scoring
# convention.
.alpha_score_matrix <- function(x, values = NULL) {
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a matrix or data frame.")
  }
  X <- as.matrix(x)
  if (!is.numeric(X)) stop("'x' must be numeric.")
  storage.mode(X) <- "double"
  if (ncol(X) < 2L) stop("coefficient alpha requires at least two items.")

  if (!is.null(values)) {
    if (!is.numeric(values)) stop("'values' must be numeric.")
    observed <- is.finite(X)
    if (!any(observed)) stop("all item responses are missing.")
    categories <- sort(unique(c(X[observed])))
    if (length(values) != length(categories)) {
      stop("Length of 'values' must equal the number of unique observed categories.")
    }
    scored <- matrix(NA_real_, nrow = nrow(X), ncol = ncol(X),
                     dimnames = dimnames(X))
    scored[observed] <- values[match(X[observed], categories)]
    X <- scored
  }
  X
}

# Surface the FIML null-space diagnostic: `null_frac` is, per coefficient,
# the fraction of the delta-method gradient lying in the truncated null space
# of the Louis information. Zero means the coefficient is an estimable
# function of the identified directions; larger values mean the saturated
# nuisance fit is selection-dependent along directions the coefficient
# touches, so warn once and keep the estimate.
.warn_null_frac <- function(null_frac, em_options = list(), threshold = 0.01) {
  if (is.null(null_frac) || length(null_frac) == 0L) return(invisible(NULL))
  # With flattening the posterior mode is unique, so the selection-dependence
  # half of the message no longer applies; the diagnostic stays available in
  # the returned object either way.
  flatten <- em_options[["flatten"]]
  if (!is.null(flatten) && is.finite(flatten) && flatten > 0) {
    return(invisible(NULL))
  }
  mx <- max(c(null_frac, 0), na.rm = TRUE)
  if (is.finite(mx) && mx > threshold) {
    warning(sprintf(paste0(
      "The saturated joint distribution is not uniquely identified from the ",
      "observed missing-data patterns (null-space gradient fraction %.3f). ",
      "The coefficient is still identified, but the point estimate depends ",
      "on which likelihood maximizer the EM selected, and the SE uses a ",
      "rank-truncated pseudo-inverse. Pass em_options = list(flatten = 0.1) ",
      "to select a unique interior posterior mode (point estimate essentially ",
      "unchanged; standard errors become conservative)."), mx), call. = FALSE)
  }
  invisible(NULL)
}

# Attach exact complete-data delete-one evaluations without changing the
# estimator object when the requested backend has a different point
# functional. The latter guard matters for structured likelihood estimators.
.attach_delete_one <- function(fit, out, coefficient_names) {
  if (is.null(out) || length(out$estimate) != length(coefficient_names)) {
    return(fit)
  }
  names(out$estimate) <- coefficient_names
  colnames(out$leave_one_out) <- coefficient_names
  if (!is.null(out$chance) && length(out$chance) > 0L) {
    names(out$chance) <- coefficient_names
    colnames(out$chance_leave_one_out) <- coefficient_names
  }
  target <- unname(fit$estimates[coefficient_names])
  scale <- max(1, abs(target), abs(out$estimate), na.rm = TRUE)
  if (any(!is.finite(target)) || any(!is.finite(out$estimate)) ||
      max(abs(target - out$estimate)) > 1e-7 * scale) {
    return(fit)
  }
  fit$delete_one <- out
  fit
}

.attach_alpha_delete_one <- function(fit, x) {
  if (nrow(x) < 3L || !all(is.finite(x))) return(fit)
  .attach_delete_one(fit, rcpp_alpha_delete_one(x), "alpha")
}

.attach_raw_kappa_delete_one <- function(fit, x, weight, values,
                                         em_options = list()) {
  if (nrow(x) < 3L || anyNA(x) || any(!is.finite(x))) return(fit)
  flatten <- em_options[["flatten"]]
  if (identical(fit$method, "fiml") && !is.null(flatten) && flatten > 0) {
    return(fit)
  }
  out <- rcpp_kappa_delete_one(x, weight, values)
  .attach_delete_one(fit, out, c("Conger", "Fleiss", "Brennan-Prediger"))
}

.attach_continuous_kappa_delete_one <- function(fit, x, weight) {
  if (nrow(x) < 3L || !all(is.finite(x))) return(fit)
  out <- rcpp_kappa_continuous_delete_one(x, weight)
  .attach_delete_one(fit, out, c("Conger", "Fleiss"))
}

.attach_counts_kappa_delete_one <- function(fit, x, weight, values,
                                            em_options = list()) {
  if (nrow(x) < 3L || anyNA(x) || any(!is.finite(x))) return(fit)
  flatten <- em_options[["flatten"]]
  if (identical(fit$method, "cat_fiml") && !is.null(flatten) && flatten > 0) {
    return(fit)
  }
  out <- rcpp_kappa_counts_delete_one(x, weight, values)
  .attach_delete_one(fit, out, c("Fleiss", "Brennan-Prediger"))
}

.alpha_from_cpp <- function(out, method, em_options = list()) {
  .warn_null_frac(out$null_frac, em_options)
  .new_estimate(out, "alpha", method = method, weight = NULL)
}

#' Coefficient alpha with missing data
#'
#' @description
#' Estimates Cronbach's coefficient alpha for a subjects-by-items matrix of
#' item scores with missing entries. Pick the estimator with `estimator`:
#' \itemize{
#'   \item `"pairwise"`: pairwise-available covariance moments under MCAR,
#'     with an influence-function sandwich standard error.
#'   \item `"cat_fiml"`: saturated-multinomial full-information maximum
#'     likelihood for finite-category items under ignorable missingness.
#'   \item `"nt_fiml"`: robust normal-theory FIML (saturated Gaussian
#'     covariance by EM, sandwich delta-method SE) for continuous items under
#'     ignorable missingness.
#'   \item `"nt_fiml_cs"`: the same normal FIML under the compound-symmetry
#'     covariance model \eqn{\Sigma = \psi J + e I} (fixed items, needs only a
#'     connected item co-observation graph). Not ten Hove's random-rater
#'     MLE-RE / ICC.
#' }
#'
#' For ordinal or otherwise categorical responses under `"pairwise"` or
#' `"nt_fiml"`, either pass the scored numeric matrix directly or pass integer
#' category codes with `values`. The observed categories are sorted and mapped
#' to `values` before estimation.
#'
#' All three estimators target fixed-item coefficient alpha. With missing
#' entries, every item must be observed at least once. The pairwise,
#' categorical-FIML, and default saturated normal-FIML estimators also require
#' every item pair to be jointly observed by at least one subject. The
#' structured normal-FIML compound-symmetry covariance model only requires a
#' connected item co-observation graph.
#'
#' @param x A subjects-by-items numeric matrix or data frame. `NA` and other
#'   non-finite values indicate missing entries. For `estimator = "cat_fiml"`,
#'   entries are integer category codes.
#' @param estimator One of `"pairwise"`, `"cat_fiml"`, `"nt_fiml"`, or
#'   `"nt_fiml_cs"` (the compound-symmetry covariance variant of `"nt_fiml"`).
#' @param values Optional numeric vector of scores for observed categories.
#'   When supplied, its length must equal the number of unique finite entries
#'   in `x`.
#' @param em_options Named list tuning the EM fit for the likelihood
#'   estimators (an error for `"pairwise"`). For `"nt_fiml"` / `"nt_fiml_cs"`:
#'   `tol` and `max_iter`. For `"cat_fiml"`: `tol`, `max_iter`,
#'   `prune_tol`, `start_alpha`, `info_rcond` (the relative
#'   eigenvalue cutoff used when inverting Louis' observed information), and
#'   `flatten` (total Dirichlet pseudo-mass spread over the complete pattern
#'   table. Any positive value makes the fitted table the unique interior
#'   posterior mode when the saturated likelihood is flat, at the cost of
#'   shrinking it toward uniform with weight `flatten / (n + flatten)`.
#'   `0`, the default, is strict ML). Flattening is a uniqueness device
#'   that leaves the point estimate essentially
#'   unchanged but makes the reported standard errors conservative
#'   (roughly 50\% too wide in calibration), so leave it at `0` unless a
#'   unique reproducible fit matters more than SE sharpness. Pass any
#'   subset.
#'
#' @return An object of class `misskappa_estimate` carrying one coefficient
#'   named `alpha` and its asymptotic covariance matrix. Methods: `print`,
#'   `coef`, `vcov`, `confint`, `as.data.frame`, and `stats::influence`.
#'   The object also carries a `psi` component (the n-by-1 matrix of
#'   per-subject influence functions).
#'
#' @examples
#' # Continuous item battery: the textual subscale of the Holzinger-Swineford
#' # (1939) data. Normal-theory FIML alpha, valid under ignorable missingness.
#' textual <- dat.holzinger1939[, c("x4", "x5", "x6")]
#' fit <- alpha(textual, estimator = "nt_fiml")
#' coef(fit)
#' confint(fit)
#'
#' # Compound-symmetry covariance variant: the same FIML under
#' # Sigma = psi J + e I, which needs only a connected item co-observation
#' # graph rather than every item pair jointly observed.
#' alpha(textual, estimator = "nt_fiml_cs")
#'
#' @examplesIf requireNamespace("psych", quietly = TRUE)
#' # Real item-level missing data: the Neuroticism scale of psych::bfi
#' # (2800 respondents, ~360 with at least one missing item). The FIML
#' # estimator is valid under ignorable missingness.
#' data(bfi, package = "psych")
#' N <- paste0("N", 1:5)
#' alpha(bfi[, N], estimator = "nt_fiml")
#'
#' @references
#' Moss, J. (in preparation). *Agreement coefficients with missing ratings.*
#' Working paper (kappa-missing); citation to be updated on preprint and
#' publication.
#'
#' @export
alpha <- function(x,
                  estimator = c("pairwise", "cat_fiml", "nt_fiml", "nt_fiml_cs"),
                  values = NULL,
                  em_options = list()) {
  estimator <- match.arg(estimator)
  em_options <- .check_em_options(em_options, switch(estimator,
    cat_fiml = c(.em_options_categorical, "flatten"),
    nt_fiml = , nt_fiml_cs = .em_options_normal,
    character(0)), estimator)
  is_cs <- estimator == "nt_fiml_cs"
  .check_pattern_identifiable(
    .pattern_observed(x), unit = "item",
    require = if (is_cs) "connected" else "complete",
    coefficient = if (is_cs) {
      "coefficient alpha under the compound-symmetry covariance model"
    } else {
      "coefficient alpha"
    })

  if (estimator == "pairwise") {
    X <- .alpha_score_matrix(x, values)
    out <- rcpp_alpha_available_continuous(X)
    fit <- .alpha_from_cpp(out, method = "pairwise")
    return(.attach_alpha_delete_one(fit, X))
  }

  if (estimator %in% c("nt_fiml", "nt_fiml_cs")) {
    X <- .alpha_score_matrix(x, values)
    em_options$covariance <- if (is_cs) "variance_components" else "saturated"
    fit <- alpha_continuous(X, em_options = em_options)
    fit$method <- estimator
    fit$covariance <- NULL    # the estimator name carries the structure
    fit$mean <- NULL
    if (!is_cs) fit <- .attach_alpha_delete_one(fit, X)
    return(fit)
  }

  fit <- alpha_cat_fiml(x, values = values, em_options = em_options)
  fit$method <- "cat_fiml"
  X <- .alpha_score_matrix(x, values)
  .attach_alpha_delete_one(fit, X)
}

#' Saturated categorical FIML coefficient alpha (internal)
#'
#' @description
#' Backend for `alpha(estimator = "cat_fiml")`. Fits the
#' saturated multinomial full-response distribution with EM, then maps that
#' distribution to the implied scored covariance matrix.
#'
#' This fixed-item functional requires every item to be observed and every item
#' pair to be jointly observed by at least one subject. The public [alpha()]
#' wrapper checks this before dispatch. Direct callers reach the C++ guard.
#'
#' @param x A subjects-by-items matrix of integer category codes. `NA`s
#'   indicate missing entries.
#' @param values Optional numeric vector of category scores. Defaults to the
#'   sorted unique observed categories.
#' @param em_options Named list of options for the categorical EM fit:
#'   `tol`, `max_iter`, `prune_tol`, `start_alpha`, `info_rcond`, `flatten`.
#'
#' @return An object of class `misskappa_estimate` carrying one coefficient
#'   named `alpha` and its asymptotic covariance matrix.
#'
#' @noRd
alpha_cat_fiml <- function(x, values = NULL, em_options = list()) {
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a matrix or data frame.")
  }
  x_mat <- as.matrix(x)
  if (!is.numeric(x_mat)) stop("'x' must be numeric.")
  storage.mode(x_mat) <- "integer"

  out <- rcpp_alpha_raw(
    x = x_mat,
    method = "fiml",
    values = values,
    em_options = em_options
  )

  .alpha_from_cpp(out, method = "alpha-cat-fiml", em_options = em_options)
}

# Build a numeric score matrix from raw ratings. With values = NULL the input
# is used as-is; with values supplied, observed categories are sorted and
# mapped to the supplied scores.
.score_matrix <- function(x, values = NULL) {
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a matrix or data frame.")
  }
  X <- as.matrix(x)
  if (!is.numeric(X)) stop("'x' must be numeric.")
  storage.mode(X) <- "double"
  if (!is.null(values)) {
    if (!is.numeric(values)) stop("'values' must be numeric.")
    observed <- is.finite(X)
    if (!any(observed)) stop("all ratings are missing.")
    categories <- sort(unique(c(X[observed])))
    if (length(values) != length(categories)) {
      stop("Length of 'values' must equal the number of unique observed categories.")
    }
    X[observed] <- values[match(X[observed], categories)]
  }
  X
}

# Internal: full categorical-kappa dispatch over every C++ method
# ("available", "ipw", "gwet", "fiml"). The public kappa() exposes only
# "ipw" / "cat_fiml"; simulations reach the dropped methods (available-case,
# Gwet) through this helper, which returns the same misskappa_estimate object.
estimate_kappa_raw <- function(x,
                               method = c("available", "ipw", "fiml", "gwet"),
                               weight = "identity",
                               values = NULL,
                               em_options = list()) {
  method <- match.arg(method)
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a matrix or data frame.")
  }
  x_mat <- as.matrix(x)
  if (!is.numeric(x_mat)) stop("'x' must be numeric.")
  storage.mode(x_mat) <- "integer"

  out <- rcpp_kappa_raw(x = x_mat, method = method, weight_type = weight,
                        values = values, em_options = em_options)
  .warn_null_frac(out$null_frac, em_options)
  fit <- .new_estimate(out, c("Conger", "Fleiss", "Brennan-Prediger"),
                       method = method, weight = weight)
  .attach_raw_kappa_delete_one(fit, x_mat, weight, values, em_options)
}

#' Multi-weight saturated categorical FIML kappa (internal)
#'
#' @description
#' Fits the saturated multinomial FIML model once, then maps the fitted
#' distribution to several agreement-weight functionals. This is intended for
#' simulations and validation code that need nominal, linear, and quadratic
#' Cat-FIML coefficients for the same data matrix. The public [kappa()] API
#' deliberately keeps one `weight` per call.
#'
#' This fixed-rater functional requires every rater to be observed and every
#' rater pair to be jointly observed by at least one subject. Otherwise the
#' requested saturated pairwise agreement functional is not identified.
#'
#' @param x A subjects-by-raters matrix of integer category codes. `NA`s
#'   indicate missing entries.
#' @param weights Character vector containing any of `"identity"`,
#'   `"nominal"`, `"linear"`, or `"quadratic"`. `"nominal"` is an alias for
#'   `"identity"`.
#' @param values Optional numeric category scores for metric weights.
#' @param em_options Named list of categorical EM options.
#'
#' @return A named list of `misskappa_estimate` objects.
#'
#' @noRd
estimate_kappa_fiml_multi <- function(x,
                                      weights = c("identity", "linear", "quadratic"),
                                      values = NULL,
                                      em_options = list()) {
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a matrix or data frame.")
  }
  x_mat <- as.matrix(x)
  if (!is.numeric(x_mat)) stop("'x' must be numeric.")
  storage.mode(x_mat) <- "integer"

  choices <- c("identity", "nominal", "linear", "quadratic")
  weights <- match.arg(weights, choices, several.ok = TRUE)
  weights <- unique(ifelse(weights == "nominal", "identity", weights))

  raw <- rcpp_kappa_fiml_multi(
    x = x_mat,
    weight_types = weights,
    values = values,
    em_options = em_options
  )

  out <- lapply(seq_along(raw), function(i) {
    if (i == 1L) .warn_null_frac(raw[[i]]$null_frac, em_options)
    .new_estimate(raw[[i]], c("Conger", "Fleiss", "Brennan-Prediger"),
                  method = "fiml", weight = weights[[i]])
  })
  names(out) <- weights
  out
}

#' Weighted agreement coefficients with missing data
#'
#' @description
#' Estimates Conger / Fleiss / Brennan-Prediger weighted agreement
#' coefficients for raw ratings with missing entries. Pick the estimator with
#' `estimator`:
#' \itemize{
#'   \item `"ipw"`: inverse-probability-weighted moment estimator for
#'     categorical ratings under MCAR, consistent under non-exchangeable raters
#'     with rater-varying observation rates.
#'   \item `"cat_fiml"`: saturated-multinomial full-information maximum
#'     likelihood for categorical ratings under ignorable missingness.
#'   \item `"pairwise"`: pairwise-available moment estimator for the
#'     quadratically weighted (scored) coefficient under MCAR.
#'   \item `"nt_fiml"`: robust normal-theory FIML (saturated covariance) for
#'     the quadratically weighted coefficient under ignorable missingness.
#'   \item `"nt_fiml_cs"`: the same FIML under the compound-symmetry
#'     covariance model \eqn{\Sigma = \psi J + e I} (fixed raters, needs only a
#'     connected rater co-observation graph).
#' }
#'
#' `"ipw"` and `"cat_fiml"` take integer category codes and the nominal,
#' linear, or quadratic `weight`. `"pairwise"`, `"nt_fiml"`, and `"nt_fiml_cs"`
#' are quadratic by construction: they treat `x` as numeric scores (pass scores
#' directly, or integer codes with `values`), and `weight` defaults to
#' `"quadratic"` for them.
#'
#' Brennan-Prediger compares the observed disagreement with a uniform draw from
#' the rating scale, so it is reported only when the scale is known. The
#' categorical estimators take the scale from the categories (the observed
#' codes, or `values`). The scored estimators `"nt_fiml"` and `"nt_fiml_cs"`
#' report it when `values` declares a discrete scale, or for vector-valued
#' ratings when `ranges` declares a continuous one. `"pairwise"` reports Conger
#' and Fleiss only.
#'
#' Two generalizations are reached through the same call. Set `g > 2` for the
#' Frechet / Hubert g-wise (multirater) kernels with `"ipw"` or `"cat_fiml"`,
#' which report the Conger-type (distinct rater combinations) and Fleiss-type
#' coefficients. Pass a subjects-by-raters-by-features array to score
#' vector-valued ratings. With `"ipw"` or `"pairwise"` and no feature metric the
#' loss is component-separable (per-feature weights via `feature_weights`, or
#' `loss = "rms"`). Pass a full feature metric `W` (or use `"nt_fiml"` /
#' `"nt_fiml_cs"`) for the squared vector loss `(x - y)' W (x - y)`:
#' `"pairwise"` is the moment estimator and `"nt_fiml"` / `"nt_fiml_cs"` the
#' saturated / block compound-symmetry normal FIML
#' (`Sigma = J_R \%x\% Psi + I_R \%x\% E`).
#'
#' For scalar fixed-rater coefficients with missing entries, the observed-data
#' pattern must identify the requested functional. For ordinary pairwise kappa
#' (`g = 2`, and all quadratic requests), every rater must be observed at least
#' once and every rater pair must be jointly observed by at least one subject
#' (`"nt_fiml_cs"` needs only a connected co-observation graph). For
#' non-quadratic `g > 2`, every requested rater g-tuple must be jointly
#' observed by at least one subject. Violations are reported as
#' non-identification before the estimator runs.
#'
#' @param x A subjects-by-raters matrix or data frame of scalar ratings, or a
#'   subjects-by-raters-by-features 3-D array of vector-valued ratings (the
#'   array form is auto-detected). Integer category codes for `"ipw"` /
#'   `"cat_fiml"`. Numeric scores for `"pairwise"` / `"nt_fiml"` /
#'   `"nt_fiml_cs"`. `NA` marks missing entries.
#' @param estimator One of `"ipw"`, `"cat_fiml"`, `"pairwise"`, `"nt_fiml"`,
#'   or `"nt_fiml_cs"`. Vector-valued ratings support `"ipw"` and `"pairwise"`
#'   (component-separable) and, with the squared loss, `"pairwise"`,
#'   `"nt_fiml"`, and `"nt_fiml_cs"`. For non-quadratic `g > 2` use `"ipw"`
#'   (which equals the complete-data estimator when nothing is missing) or
#'   `"cat_fiml"`.
#' @param weight Weighting scheme: `"nominal"` (the default for `"ipw"` /
#'   `"cat_fiml"`), `"linear"`, or `"quadratic"` (the only, and default, choice
#'   for the scored estimators). `"hubert"` is the all-raters-equal multirater
#'   kernel for `g > 2`. For vector-valued ratings `weight` selects the
#'   component loss (`"nominal"` -> Hamming, `"linear"` -> L1, `"quadratic"` ->
#'   squared).
#' @param g Arity of the multirater disagreement kernel. `g = 2` (the default)
#'   is ordinary pairwise kappa. `g > 2` uses the Frechet / Hubert g-wise
#'   family. `weight = "quadratic"` is g-invariant, so `g` is ignored there and
#'   the closed form is used. Not used for vector-valued ratings.
#' @param values Optional numeric vector of category scores. The sorted
#'   distinct observed ratings are mapped to `values` in order, so its length
#'   must equal their number. For `"ipw"` / `"cat_fiml"` the scores
#'   parameterize the linear and quadratic weights (default: the codes
#'   themselves). For the scored estimators they convert codes to scores and
#'   declare the discrete rating scale that Brennan-Prediger needs. Not used
#'   for vector-valued ratings.
#' @param em_options Named list of EM options for the likelihood estimators
#'   (an error for `"ipw"` and `"pairwise"`). For `"cat_fiml"`: `tol`,
#'   `max_iter`, `prune_tol`, `start_alpha`, `info_rcond`, and `flatten` (total
#'   Dirichlet pseudo-mass spread over the complete pattern table. Any positive
#'   value selects the unique interior posterior mode when the saturated
#'   likelihood is flat, shrinking the fitted table toward uniform with weight
#'   `flatten / (n + flatten)`. `0`, the default, is strict ML). Flattening is
#'   a uniqueness device that leaves the point estimate essentially unchanged
#'   but makes the reported standard errors conservative (roughly 50\% too wide
#'   in calibration). For `"nt_fiml"` / `"nt_fiml_cs"`: `tol` and `max_iter`.
#'   Pass any subset.
#' @param ... Mode-specific arguments. Anything not listed here is an error.
#'   For vector-valued ratings: `W`, a features-by-features symmetric
#'   positive-semidefinite feature metric for the squared loss (with
#'   `"pairwise"`, `"nt_fiml"`, `"nt_fiml_cs"`; the FIML estimators default to
#'   the identity); `ranges`, a features-by-2 matrix of `(lo, hi)` score bounds
#'   declaring a continuous rating scale for Brennan-Prediger (with
#'   `"nt_fiml"` / `"nt_fiml_cs"`); and `feature_weights` / `loss` for the
#'   component-separable `"ipw"` / `"pairwise"` path (non-negative per-feature
#'   weights, and `loss = "rms"` for the root-mean-square loss). For `g > 2`:
#'   `max_chance_tuples`, the cap on the number of chance tuples enumerated
#'   (default `5e6`).
#'
#' @details
#' The R interface exposes the main weighting schemes used in current
#' agreement work: `"nominal"`, `"linear"`, and `"quadratic"` (plus the
#' `"hubert"` g-wise kernel for `g > 2`). Older Krippendorff / irrCAC-style
#' categorical schemes (`"ordinal"`, `"radical"`, `"ratio"`, `"circular"`,
#' `"bipolar"`) remain in the C++ library for validation and parity work.
#'
#' @references
#' Moss, J. (in preparation). *Agreement coefficients with missing ratings.*
#' Working paper (kappa-missing); citation to be updated on preprint and
#' publication.
#'
#' Moss, J. (in preparation). *Quadratically weighted agreement coefficients.*
#' Working paper (quadratic); citation to be updated on preprint and
#' publication.
#'
#' Krippendorff, K. (2011). Computing Krippendorff's alpha-reliability.
#'
#' Gwet, K. L. (2019). irrCAC: Computing chance-corrected agreement
#' coefficients among raters.
#'
#' @return An object of class `misskappa_estimate` carrying the named
#'   coefficient estimates, their asymptotic covariance matrix, and the number
#'   of subjects `n`. Methods: `print`, `coef`, `vcov`, `confint`,
#'   `as.data.frame`, and `stats::influence`. The object also carries a `psi`
#'   component (the n-by-K matrix of per-subject influence functions,
#'   satisfying `vcov == crossprod(psi) / n^2`) for building custom contrasts
#'   or joint tests.
#'
#' @examples
#' # Categorical ratings with missing entries (Gwet 2014): the
#' # inverse-probability-weighted estimator under MCAR, with linear weights.
#' kappa(dat.gwet2014, estimator = "ipw", weight = "linear")
#'
#' # Scored coefficient (quadratic loss) via the pairwise-available moment
#' # estimator. Treats the columns as numeric scores.
#' kappa(dat.zapf2016, estimator = "pairwise")
#'
#' # Normal-theory FIML. Declaring the 1-5 rating scale with `values` adds
#' # Brennan-Prediger.
#' kappa(dat.zapf2016, estimator = "nt_fiml", values = 1:5)
#'
#' # Compound-symmetry covariance variant: the fixed-rater variance-component
#' # model. Consistent and slightly more efficient than "nt_fiml" when the
#' # structure holds, and identified by a connected (not necessarily complete)
#' # rater co-observation graph.
#' kappa(dat.zapf2016, estimator = "nt_fiml_cs")
#'
#' # Planned-incomplete design (synthetic): some rater pairs are never
#' # co-observed, so the saturated estimators are not identified, but the
#' # compound-symmetry estimator needs only a connected co-observation graph.
#' kappa(dat.tenhove2025, estimator = "nt_fiml_cs")
#'
#' # g-wise (multirater) nominal kappa over triples of raters.
#' kappa(dat.gwet2014, estimator = "ipw", weight = "nominal", g = 3)
#'
#' # Vector-valued ratings: a subjects-by-raters-by-features array. Each of the
#' # four expert observers classifies crackles at six chest sites, so a rating
#' # is a six-vector. The default Hamming loss counts disagreeing sites.
#' kappa(dat.vanbelle2019[, 1:4, ], estimator = "ipw")
#'
#' # Full feature-metric quadratic vector agreement under the block
#' # compound-symmetry normal FIML, with a feature metric W over the six sites.
#' sites <- dat.vanbelle2019[, 1:4, ]
#' kappa(sites, estimator = "nt_fiml_cs", W = diag(dim(sites)[3]))
#'
#' @export
kappa <- function(x,
                  estimator = c("ipw", "cat_fiml", "pairwise", "nt_fiml",
                                "nt_fiml_cs"),
                  weight = c("nominal", "linear", "quadratic", "hubert"),
                  g = 2L,
                  values = NULL,
                  em_options = list(),
                  ...) {
  estimator <- match.arg(estimator)
  weight_supplied <- !missing(weight)
  weight <- match.arg(weight)
  dots <- list(...)
  em_allowed <- switch(estimator,
    cat_fiml = c(.em_options_categorical, "flatten"),
    nt_fiml = , nt_fiml_cs = .em_options_normal,
    character(0))

  ## Vector-valued ratings: subjects-by-raters-by-features array.
  if (is.array(x) && length(dim(x)) == 3L) {
    if (!missing(g)) {
      stop("'g' does not apply to vector-valued ratings.", call. = FALSE)
    }
    if (!is.null(values)) {
      stop("'values' does not apply to vector-valued ratings; use 'ranges' ",
           "to declare the score range.", call. = FALSE)
    }
    em_options <- .check_em_options(em_options, em_allowed, estimator)
    return(.kappa_vector_dispatch(x, estimator, weight, weight_supplied,
                                  em_options, dots))
  }

  if (!is.numeric(g) || length(g) != 1L || !is.finite(g) || g < 2 ||
      g != round(g)) {
    stop("'g' must be an integer >= 2.")
  }
  g <- as.integer(g)
  if (weight == "hubert" && g == 2L) {
    stop('weight = "hubert" is only defined for g > 2; use weight = ',
         '"nominal" for g = 2.')
  }

  ## Quadratic weighting is a moment functional and therefore g-invariant: any
  ## g collapses to the same coefficient, so always take the closed form
  ## (the g = 2 path) and never the combinatorial g-wise enumeration.
  if (weight == "quadratic" || g == 2L) {
    .check_dots(dots, character(0), "scalar ratings with g = 2")
    em_options <- .check_em_options(em_options, em_allowed, estimator)
    ## Every public scalar g = 2 path surfaces a fixed-rater coefficient, so
    ## the saturated functional needs a complete rater co-observation graph.
    ## The compound-symmetry normal FIML only needs a connected one; the
    ## reduced information check in C++ remains the final guard.
    is_cs <- estimator == "nt_fiml_cs"
    .check_pattern_identifiable(
      .pattern_observed(x), unit = "rater",
      require = if (is_cs) "connected" else "complete",
      coefficient = if (is_cs) {
        "Conger's kappa under the compound-symmetry covariance model"
      } else {
        "Conger's kappa (and Brennan-Prediger)"
      },
      arity = 2L)
    return(.kappa_scalar_g2(x, estimator, weight, weight_supplied, values,
                            em_options))
  }

  ## g > 2, non-quadratic: the Frechet / Hubert g-wise family.
  .check_dots(dots, "max_chance_tuples", "g-wise kappa")
  if (estimator == "cat_fiml") em_allowed <- .em_options_categorical
  em_options <- .check_em_options(em_options, em_allowed, estimator)
  .check_pattern_identifiable(
    .pattern_observed(x), unit = "rater",
    coefficient = sprintf("%s-wise kappa", g),
    arity = g)
  .kappa_scalar_gwise(x, estimator, weight, g, values, em_options, dots)
}

# Ordinary pairwise (g = 2) kappa, also reached by any quadratic request
# because quadratic kappa does not depend on g.
.kappa_scalar_g2 <- function(x, estimator, weight, weight_supplied, values,
                             em_options) {
  if (estimator %in% c("pairwise", "nt_fiml", "nt_fiml_cs")) {
    if (weight_supplied && weight != "quadratic") {
      stop(sprintf('estimator = "%s" requires weight = "quadratic".', estimator))
    }
    X <- .score_matrix(x, values)
    if (estimator == "pairwise") {
      fit <- kappa_continuous(X, method = "available", weight = "quadratic")
    } else {
      em_options$covariance <- if (estimator == "nt_fiml_cs") {
        "variance_components"
      } else {
        "saturated"
      }
      fit <- kappa_quadratic_fiml(X, values = values, em_options = em_options)
      fit$covariance <- NULL    # the estimator name carries the structure
      fit$mean <- NULL
      # Brennan-Prediger needs a declared scale (see ?kappa).
      if (is.null(values)) fit <- .drop_coefficients(fit, "Brennan-Prediger")
    }
    fit$method <- estimator
    fit$weight <- "quadratic"
    return(fit)
  }

  cpp_weight <- if (weight == "nominal") "identity" else weight
  fit <- estimate_kappa_raw(
    x,
    method = if (estimator == "ipw") "ipw" else "fiml",
    weight = cpp_weight,
    values = values,
    em_options = em_options
  )
  fit$method <- estimator
  fit$weight <- weight
  fit
}

# g > 2, non-quadratic. Maps the public weight to a g-wise kernel and the
# public estimator to the backend method, then delegates to kappa_gwise().
.kappa_scalar_gwise <- function(x, estimator, weight, g, values, em_options,
                                dots) {
  map <- switch(weight,
    nominal = list(distance = "nominal", family = "categorical"),
    hubert  = list(distance = "hubert",  family = "categorical"),
    linear  = list(distance = "absolute", family = "continuous"),
    stop(sprintf(
      'weight = "%s" has no g-wise (g > 2) kernel; use "nominal", "linear", "hubert", or "quadratic".',
      weight))
  )

  method <- switch(estimator,
    ipw      = "ipw",
    cat_fiml = if (map$family == "categorical") "fiml" else stop(
      'estimator = "cat_fiml" has no continuous g-wise variant; use ',
      'weight = "nominal"/"hubert", or estimator = "ipw".', call. = FALSE),
    stop(sprintf(paste0(
      'estimator = "%s" is a quadratic estimator; for g > 2 with weight = ',
      '"%s" use "ipw" (equal to the complete-data estimator when no ratings ',
      'are missing) or "cat_fiml".'), estimator, weight), call. = FALSE)
  )

  max_chance_tuples <- if (!is.null(dots$max_chance_tuples)) {
    dots$max_chance_tuples
  } else {
    5000000L
  }
  xx <- if (map$family == "continuous") .score_matrix(x, values) else x

  fit <- kappa_gwise(xx, distance = map$distance, method = method, g = g,
                     em_options = em_options,
                     max_chance_tuples = max_chance_tuples)
  fit$method <- estimator
  fit$weight <- weight
  fit
}

# Vector-valued ratings (subjects-by-raters-by-features array).
.kappa_vector_dispatch <- function(x, estimator, weight, weight_supplied,
                                   em_options, dots) {
  W <- dots$W
  ranges <- dots$ranges
  # Full-matrix-W quadratic family: the squared vector loss (x - y)' W (x - y).
  # The normal-FIML covariance models ("nt_fiml" saturated, "nt_fiml_cs" block
  # compound symmetry) report Conger and Fleiss, plus Brennan-Prediger when
  # `ranges` declares the score range.
  if (estimator %in% c("nt_fiml", "nt_fiml_cs")) {
    .check_dots(dots, c("W", "ranges"),
                sprintf('vector-valued estimator = "%s"', estimator))
    if (weight_supplied && weight != "quadratic") {
      stop(sprintf('estimator = "%s" requires weight = "quadratic".', estimator))
    }
    fit <- kappa_vector_quadratic(x, method = estimator, W = W,
                                  em_options = em_options, ranges = ranges)
    if (is.null(ranges)) fit <- .drop_coefficients(fit, "Brennan-Prediger")
    fit$method <- estimator
    fit$weight <- "quadratic"
    return(fit)
  }
  if (estimator == "pairwise" && !is.null(W)) {
    .check_dots(dots, "W",
                'vector-valued estimator = "pairwise" with a feature metric W')
    if (weight_supplied && weight != "quadratic") {
      stop('a full feature metric W requires weight = "quadratic".')
    }
    fit <- kappa_vector_quadratic(x, method = "pairwise", W = W)
    fit$method <- estimator
    fit$weight <- "quadratic"
    return(fit)
  }
  if (!estimator %in% c("pairwise", "ipw")) {
    stop(sprintf(
      'vector-valued ratings support estimator = "ipw", "pairwise", "nt_fiml", or "nt_fiml_cs"; got "%s".',
      estimator), call. = FALSE)
  }
  # Component-separable (diagonal-weight) path.
  .check_dots(dots, c("feature_weights", "loss"),
              sprintf('component-separable vector estimator = "%s"', estimator))
  loss <- if (!is.null(dots$loss)) {
    dots$loss
  } else {
    switch(weight,
      nominal   = "hamming",
      linear    = "absolute",
      quadratic = "squared",
      stop(sprintf(
        'weight = "%s" has no vector component loss; use "nominal", "linear", or "quadratic", or pass loss = "rms".',
        weight))
    )
  }
  fit <- kappa_vector(x, method = estimator, loss = loss,
                      feature_weights = dots$feature_weights)
  fit$method <- estimator
  fit$weight <- if (!is.null(dots$loss)) dots$loss else weight
  fit
}

#' Closed-form quadratic-loss kappa estimator
#'
#' @description
#' Closed-form moment-based estimator for the quadratically-weighted
#' Conger / Fleiss / Brennan-Prediger family. Treats categorical ratings
#' as numeric scores and uses per-rater means and covariances.
#'
#' Use this when the quadratic weighting is the right loss and a parametric
#' moment-based estimate is preferred over the U-statistic version
#' (`kappa(x, estimator = "pairwise")`). The two
#' agree under the usual asymptotics. The closed form can be faster on
#' large `n`.
#'
#' With missing ratings, every rater must be observed at least once and every
#' rater pair must be jointly observed by at least one subject.
#'
#' @param x A subjects-by-raters numeric matrix of category scores. `NA`
#'   marks missing entries.
#' @param values Length-C numeric vector of category scores. The quadratic
#'   loss is `(values[i] - values[j])^2 / (max - min)^2`.
#'   The covariance is the empirical covariance of the reduced quadratic
#'   moment summaries from their row-wise estimating equations.
#'
#' @return A `misskappa_estimate` object with `Conger`, `Fleiss`,
#'   `Brennan-Prediger` coefficients, the 3x3 vcov, and per-subject
#'   influence-function rows.
#'
#' @noRd
kappa_quadratic <- function(x, values) {
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a matrix or data frame.")
  }
  x_mat <- as.matrix(x)
  storage.mode(x_mat) <- "double"
  if (!is.numeric(values)) stop("'values' must be numeric.")

  out <- rcpp_kappa_quadratic(
    x = x_mat,
    values = as.numeric(values)
  )
  .new_estimate(out, c("Conger", "Fleiss", "Brennan-Prediger"),
                method = "quadratic", weight = "quadratic")
}

#' Convert rater-identified ratings to counts format
#'
#' @description
#' Aggregates a subjects-by-raters matrix of categorical ratings into the
#' subjects-by-categories counts format consumed by [kappa_counts()]: each
#' row is one subject, each column one category, and `out[i, k]` is the number
#' of raters who assigned subject `i` to category `k`. Missing entries (`NA`
#' or, for numeric input, any non-finite value) are simply not counted, so a
#' subject rated by a subset of raters contributes a row whose total is the
#' number of raters who actually scored it.
#'
#' Aggregating to counts **discards rater identity**, which is exactly the
#' exchangeable-raters representation. Feeding the result to
#' `kappa_counts(estimator = "cat_fiml")` therefore fits the exchangeable
#' count-FIML: the maximum-likelihood estimator one obtains by imposing rater
#' exchangeability on the saturated categorical FIML of
#' `kappa(estimator = "cat_fiml")`. Under genuinely non-exchangeable raters
#' (rater-specific marginals), prefer the rater-identified
#' `kappa(x, estimator = "cat_fiml")` instead, which keeps per-rater
#' information.
#'
#' @param x A subjects-by-raters matrix or data frame of category labels
#'   (integer codes, characters, or factors). `NA` (and, for numeric
#'   input, any non-finite value) marks a missing rating.
#' @param categories Optional vector giving the full set of categories and
#'   their column order in the output. Defaults to the sorted unique observed
#'   labels. Supply it to pin a known category set (including categories that
#'   happen to be unobserved, which then appear as all-zero columns) or to
#'   control column order. Observed labels outside `categories` are an error.
#'
#' @return An integer subjects-by-categories matrix whose columns are named by
#'   `categories` and whose row names are carried over from `x`. Suitable as
#'   the `x` argument of [kappa_counts()].
#'
#' @seealso [kappa_counts()] for the count-format estimators.
#'   [kappa()] for the rater-identified (non-exchangeable) estimators.
#'
#' @examples
#' # Gwet (2014): 20 subjects, five raters, categories {0, 1, 2, 3}, with
#' # missing ratings. Aggregate to counts, then fit the exchangeable
#' # count-FIML (equivalent to exchangeability-constrained cat_fiml).
#' counts <- ratings_to_counts(dat.gwet2014)
#' head(counts)
#' kappa_counts(counts, estimator = "cat_fiml")
#'
#' # The Fleiss-Cuzick moment estimator is available from the same counts.
#' kappa_counts(counts, estimator = "fleiss_cuzick")
#'
#' @export
ratings_to_counts <- function(x, categories = NULL) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x)) {
    stop("'x' must be a subjects-by-raters matrix or data frame.")
  }
  if (is.numeric(x)) x[!is.finite(x)] <- NA

  observed <- x[!is.na(x)]
  if (length(observed) == 0L) {
    stop("'x' has no observed (non-missing) ratings.")
  }

  if (is.null(categories)) {
    categories <- sort(unique(as.vector(observed)))
  } else {
    extra <- setdiff(unique(as.vector(observed)), categories)
    if (length(extra) > 0L) {
      shown <- extra[seq_len(min(5L, length(extra)))]
      stop("Observed labels absent from 'categories': ",
           paste(shown, collapse = ", "),
           if (length(extra) > 5L) ", ..." else "")
    }
  }
  C <- length(categories)

  idx <- matrix(match(x, categories), nrow = nrow(x))  # NA ratings stay NA
  counts <- matrix(0L, nrow = nrow(x), ncol = C,
                   dimnames = list(rownames(x), as.character(categories)))
  for (i in seq_len(nrow(x))) {
    row_idx <- idx[i, !is.na(idx[i, ])]
    if (length(row_idx) > 0L) counts[i, ] <- tabulate(row_idx, nbins = C)
  }
  counts
}

#' Weighted agreement coefficients for counts-format input
#'
#' @description
#' Estimates Fleiss and Brennan-Prediger weighted agreement coefficients
#' for counts-format data: each row is one subject, each column is a
#' category, and `x[i, k]` is the number of raters who assigned subject `i`
#' to category `k`. Row sums (number of raters per subject) need not be
#' uniform across subjects.
#'
#' The count-format moment estimator follows the Fleiss-Cuzick unequal-judges
#' convention: observed row disagreement is weighted by `r_i - 1`, and chance
#' disagreement uses the pooled rating-token margin. This differs from the
#' unit-weighted distribution/count convention used by some software. The
#' unit-weighted comparator is kept only in the C++ API for validation studies.
#'
#' Rows with fewer than two observed ratings contain no within-subject pair
#' information and receive zero observed-disagreement weight. At least one row
#' with two or more ratings is required.
#'
#' Counts data discards rater identity, so all counts estimators impose
#' exchangeability. The categorical count-FIML model treats full count
#' compositions as exchangeable iid rater draws. The compound-symmetry
#' estimators (`"nt_fiml_cs"` and the moment estimator `"fleiss_cuzick"`) fit
#' the corresponding common-mean, compound-symmetry scored model
#' \eqn{\Sigma = \psi J + e I}. For non-exchangeable raters, use rater-identified
#' data and `kappa()` instead.
#'
#' The quadratic compound-symmetry likelihood estimator (`"nt_fiml_cs"`) works
#' **directly on the counts**: each subject contributes only the sufficient
#' statistics \eqn{(m_i, \sum_c n_{ic} v_c, \sum_c n_{ic} v_c^2)} of its observed
#' ratings, so no realizing subjects-by-raters table is constructed. The
#' coefficients are the `r_total`-invariant maps
#' \eqn{\kappa_{\mathrm{Fleiss}} = \psi/(\psi+e)} and
#' \eqn{\kappa_{\mathrm{BP}} = 1 - 2e/c_{\mathrm{BP}}}, and every route returns
#' robust (sandwich) standard errors.
#'
#' To obtain counts from a rater-identified subjects-by-raters matrix, use
#' [ratings_to_counts()]. Chaining it into `estimator = "cat_fiml"` fits the
#' exchangeable categorical count-FIML.
#'
#' Conger is not reported (raters are not identified in this input format),
#' and neither IPW nor Gwet are meaningful (per-rater observation rates are
#' aggregated away by the counts representation).
#'
#' @param x A subjects-by-categories non-negative integer matrix.
#' @param estimator One of `"fleiss_cuzick"` (the Fleiss-Cuzick count-format
#'   moment estimator, the available-case / moment member of the
#'   compound-symmetry family), `"cat_fiml"` (EM over the composition simplex
#'   with multivariate hypergeometric weights for completing partial counts), or
#'   `"nt_fiml_cs"` (the exchangeable common-mean compound-symmetry normal FIML,
#'   quadratic weight only, the efficient ignorable-missingness MLE). All agree
#'   when every row of `x` sums to `r_total`. With partial counts (some
#'   `r_i < r_total`) `"cat_fiml"` and `"nt_fiml_cs"` are the efficient routes and
#'   `"fleiss_cuzick"` is the consistent-under-MCAR moment estimator. For a
#'   complete-case (listwise) analysis, drop the incomplete rows yourself and
#'   pass the result to any estimator.
#' @param weight Weighting scheme: `"nominal"` (default), `"linear"`, or
#'   `"quadratic"`. `"nt_fiml_cs"` is quadratic only, and its default.
#' @param values Optional length-C numeric vector of category scores used
#'   by the metric weightings. Defaults to `1:C`.
#' @param r_total Total number of raters per subject. Defaults to the
#'   maximum observed row sum. Required for `"cat_fiml"` when rows have varying
#'   totals and the maximum observed total is not the design total.
#'   `"fleiss_cuzick"` and `"nt_fiml_cs"` are `r_total`-invariant and ignore it.
#' @param em_options Named list of EM/optimiser options (an error for
#'   `"fleiss_cuzick"`). For `estimator = "cat_fiml"`: `tol`, `max_iter`,
#'   `prune_tol`, `start_alpha`, `info_rcond`. The raw-data `flatten` option
#'   does not apply there: count-data categorical FIML fits a low-dimensional
#'   composition model, not the saturated joint. For `estimator =
#'   "nt_fiml_cs"`: `tol` and `max_iter` tune the optimiser.
#'
#' @references
#' Moss, J. (in preparation). *Quadratically weighted agreement coefficients.*
#' Working paper (quadratic); citation to be updated on preprint and
#' publication.
#'
#' Krippendorff, K. (2011). Computing Krippendorff's alpha-reliability.
#'
#' Gwet, K. L. (2019). irrCAC: Computing chance-corrected agreement
#' coefficients among raters.
#'
#' @return A `misskappa_estimate` object with `Fleiss` and
#'   `Brennan-Prediger` coefficients, the 2x2 vcov, the number of subjects
#'   `n`, and a `psi` component of per-subject influence rows.
#'
#' @examples
#' # Fleiss (1971) psychiatric diagnoses: 30 subjects, 6 raters, 5 categories,
#' # in counts format (one row per subject, one column per category).
#' kappa_counts(dat.fleiss1971, estimator = "fleiss_cuzick")
#'
#' # Treat the categories as ordered scores with a quadratic loss.
#' kappa_counts(dat.fleiss1971, estimator = "fleiss_cuzick", weight = "quadratic")
#'
#' # Exchangeable common-mean compound-symmetry FIML, fit directly on the
#' # per-subject counts (quadratic weight).
#' kappa_counts(dat.fleiss1971, estimator = "nt_fiml_cs")
#'
#' @export
kappa_counts <- function(x,
                         estimator = c("fleiss_cuzick", "cat_fiml", "nt_fiml_cs"),
                         weight = c("nominal", "linear", "quadratic"),
                         values = NULL,
                         r_total = NULL,
                         em_options = list()) {
  estimator <- match.arg(estimator)
  weight_supplied <- !missing(weight)
  weight <- match.arg(weight)
  if (estimator == "nt_fiml_cs") {
    if (weight_supplied && weight != "quadratic") {
      stop('estimator = "nt_fiml_cs" requires weight = "quadratic".',
           call. = FALSE)
    }
    weight <- "quadratic"
  }
  em_options <- .check_em_options(em_options, switch(estimator,
    cat_fiml = .em_options_categorical,
    nt_fiml_cs = .em_options_normal,
    character(0)), estimator)
  cpp_weight <- if (weight == "nominal") "identity" else weight
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a matrix or data frame.")
  }
  x_mat <- as.matrix(x)
  if (!is.numeric(x_mat)) stop("'x' must be numeric.")
  storage.mode(x_mat) <- "integer"

  if (estimator == "nt_fiml_cs") {
    # Direct-on-counts exchangeable CS quadratic kappa (no realizing table;
    # see kappa_counts_cs.R).
    return(.kappa_counts_cs(x_mat, values, em_options))
  }
  out <- if (estimator == "fleiss_cuzick") {
    rcpp_kappa_counts(x = x_mat, weight_type = cpp_weight, values = values)
  } else {
    if (is.null(r_total)) r_total <- as.integer(max(rowSums(x_mat)))
    rcpp_kappa_fiml_counts(
      x = x_mat,
      weight_type = cpp_weight,
      values = values,
      r_total = as.integer(r_total),
      em_options = em_options
    )
  }
  fit <- .new_estimate(out, c("Fleiss", "Brennan-Prediger"),
                       method = estimator, weight = weight)
  .attach_counts_kappa_delete_one(fit, x_mat, cpp_weight, values, em_options)
}

#' Weighted agreement coefficients for continuous ratings
#'
#' @description
#' Estimates Conger / Fleiss weighted agreement coefficients for raw
#' real-valued ratings under missing data. The MCAR moment routes
#' (`"available"`, `"ipw"`, `"gwet"`) cover any continuous loss kernel. The
#' normal-FIML route (`"fiml"`, quadratic kernel only) fits a multivariate-
#' normal covariance by EM and is valid under ignorable (MCAR or MAR)
#' missingness. The default covariance model is saturated. `em_options =
#' list(covariance = "variance_components")` fits the compound-symmetry
#' (fixed-rater) covariance model. Missing entries are encoded as `NA` (or any
#' non-finite value).
#'
#' For the saturated missing-data routes here, every rater must be observed at
#' least once and every rater pair must be jointly observed by at least one
#' subject. The compound-symmetry normal-FIML covariance model only requires
#' a connected rater co-observation graph.
#'
#' @param x A subjects-by-raters numeric matrix.
#' @param method One of `"available"`, `"ipw"`, `"gwet"`, or `"fiml"`.
#'   `"fiml"` requires `weight = "quadratic"` and dispatches to
#'   [kappa_quadratic_fiml()].
#' @param weight Continuous loss kernel: `"identity"` (binary 0/1),
#'   `"linear"`, `"quadratic"`, `"radical"`, or `"ratio"`. The data range
#'   `[min(x), max(x)]` is used to parameterise the kernel. Only
#'   `"quadratic"` is available with `method = "fiml"`.
#' @param values Optional finite numeric score vector used only for
#'   `method = "fiml"` to define Brennan-Prediger's uniform chance baseline.
#'   If `NULL`, the continuous uniform baseline on the observed score range is
#'   used.
#' @param em_options Used only when `method = "fiml"`: named list with `tol`,
#'   `max_iter`, `covariance`, and `mean`. `covariance = "saturated"` fits the
#'   unrestricted covariance (default). `"variance_components"` fits
#'   the compound-symmetry (fixed-rater) covariance model
#'   \eqn{\Sigma = \psi J + e I}. `mean = "free"` estimates one mean per
#'   rater/item (default). `mean = "equal"` estimates a common mean and
#'   requires `covariance = "variance_components"`.
#'
#' @return A `misskappa_estimate` object carrying named coefficients
#'   (`Conger`, `Fleiss`, `Brennan-Prediger`) and the 3x3 asymptotic
#'   covariance matrix for `method = "fiml"`. Other continuous methods carry
#'   `Conger` and `Fleiss`.
#'
#' @noRd
kappa_continuous <- function(x,
                             method = c("available", "ipw", "gwet", "fiml"),
                             weight = c("quadratic", "linear", "identity",
                                        "unweighted", "radical", "ratio"),
                             values = NULL,
                             em_options = list()) {
  method <- match.arg(method)
  weight <- match.arg(weight)

  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a matrix or data frame.")
  }
  x_mat <- as.matrix(x)
  if (!is.numeric(x_mat)) stop("'x' must be numeric.")
  storage.mode(x_mat) <- "double"

  if (method == "fiml") {
    if (weight != "quadratic") {
      stop("method = \"fiml\" is only defined for the quadratic kernel; the ",
           "normal-FIML coefficients are functions of the mean and covariance ",
           "only under quadratic weighting.")
    }
    return(kappa_quadratic_fiml(x_mat, values = values, em_options = em_options))
  }

  out <- rcpp_kappa_continuous(
    x = x_mat,
    method = method,
    weight_type = weight
  )
  fit <- .new_estimate(out, c("Conger", "Fleiss"), method = method,
                       weight = weight)
  .attach_continuous_kappa_delete_one(fit, x_mat, weight)
}

#' Component-separable vector agreement coefficients
#'
#' @description
#' Internal pilot for pairwise vector-valued ratings with component-wise
#' missingness. Input is a rectangular subjects-by-raters-by-features array.
#' The estimator forms diagonal-weighted component-loss moments and reports
#' Conger and Fleiss coefficients.
#'
#' With IPW, each positive-weight rater-feature cell must be observed at
#' least once. A vector coefficient is undefined when no positive-weight rater
#' pairs are jointly observed within subjects.
#'
#' @param x Numeric array with dimensions subjects, raters, features.
#' @param method `"pairwise"` for observed component-pair moments or `"ipw"`
#'   for inverse-probability weighting by rater-feature and
#'   rater-pair-feature observation rates.
#' @param loss Component-separable loss: `"hamming"`, `"absolute"`,
#'   `"squared"`, or `"rms"`.
#' @param feature_weights Non-negative diagonal feature weights. Defaults to
#'   equal weights.
#'
#' @return A `misskappa_estimate` object with `Conger` and `Fleiss`
#'   coefficients and a 2x2 influence-function covariance matrix.
#'
#' @noRd
kappa_vector <- function(x,
                         method = c("pairwise", "ipw"),
                         loss = c("hamming", "absolute", "squared", "rms"),
                         feature_weights = NULL) {
  method <- match.arg(method)
  loss <- match.arg(loss)
  dims <- dim(x)
  if (length(dims) != 3L) {
    stop("'x' must be a subjects-by-raters-by-features array.")
  }
  if (dims[1L] < 1L || dims[2L] < 2L || dims[3L] < 1L) {
    stop("'x' must have at least one subject, two raters, and one feature.")
  }
  if (!is.numeric(x)) stop("'x' must be numeric.")
  storage.mode(x) <- "double"

  n <- dims[1L]
  R <- dims[2L]
  p <- dims[3L]
  if (is.null(feature_weights)) {
    feature_weights <- rep(1, p)
  }
  if (!is.numeric(feature_weights) || length(feature_weights) != p ||
      any(!is.finite(feature_weights)) || any(feature_weights < 0) ||
      !any(feature_weights > 0)) {
    stop("'feature_weights' must be a non-negative numeric vector of length ",
         "dim(x)[3] with at least one positive entry.")
  }

  flat <- matrix(NA_real_, nrow = n, ncol = R * p)
  for (r in seq_len(R)) {
    for (l in seq_len(p)) {
      flat[, (r - 1L) * p + l] <- x[, r, l]
    }
  }

  out <- rcpp_kappa_vector(
    x = flat,
    features = as.integer(p),
    method = method,
    loss_type = loss,
    feature_weights = as.numeric(feature_weights)
  )

  .new_estimate(out, c("Conger", "Fleiss"),
                method = paste0("vector_", method), weight = loss,
                loss = loss, feature_weights = as.numeric(feature_weights))
}

#' G-wise agreement coefficients for rectangular ratings
#'
#' @description
#' Estimates g-wise Conger-type and Fleiss-type agreement coefficients using
#' multirater disagreement kernels (the Frechet / Hubert family). The
#' `"complete"` estimator requires every entry observed. `"ipw"` and `"fiml"`
#' admit missing entries (FIML is categorical only) when every requested rater
#' g-tuple is jointly observed by at least one subject.
#'
#' @param x A subjects-by-raters matrix or data frame.
#' @param distance Multirater disagreement kernel. `"nominal"` uses Frechet
#'   mode disagreement for categorical ratings. `"absolute"` uses median
#'   absolute deviation. `"quadratic"` uses mean squared deviation. `"hubert"`
#'   uses all-raters-equal disagreement.
#' @param method `"complete"` (complete-data), `"ipw"` (MCAR inverse-probability
#'   weighting), or `"fiml"` (categorical full-information ML under ignorable
#'   missingness, not defined for the continuous kernels).
#' @param g Arity of the multirater distance. Defaults to all raters
#'   (`ncol(x)`). Must be between 2 and `ncol(x)`.
#' @param em_options Named list of EM options used by `method = "fiml"`.
#' @param max_chance_tuples Maximum number of direct `n^g` item tuples to
#'   evaluate before stopping for continuous distances, or finite category
#'   tuples for categorical distances.
#'
#' @return A `misskappa_estimate` object with `Conger` and `Fleiss`
#'   coefficients and a 2x2 influence-function covariance matrix.
#'
#' @noRd
kappa_gwise <- function(x,
                        distance = c("nominal", "absolute", "quadratic", "hubert"),
                        method = c("complete", "ipw", "fiml"),
                        g = NULL,
                        em_options = list(),
                        max_chance_tuples = 5000000L) {
  distance <- match.arg(distance)
  method <- match.arg(method)
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a matrix or data frame.")
  }
  x_mat <- as.matrix(x)
  if (!is.numeric(x_mat)) stop("'x' must be numeric.")
  if (is.null(g)) g <- ncol(x_mat)
  if (!is.numeric(g) || length(g) != 1L || !is.finite(g) ||
      g < 2 || g > ncol(x_mat)) {
    stop("'g' must be an integer between 2 and ncol(x).")
  }
  if (!is.numeric(max_chance_tuples) || length(max_chance_tuples) != 1L ||
      !is.finite(max_chance_tuples) || max_chance_tuples < 1 ||
      max_chance_tuples > .Machine$integer.max) {
    stop("'max_chance_tuples' must be a positive integer.")
  }

  categorical <- distance %in% c("nominal", "hubert")
  if (method == "fiml" && !categorical) {
    stop('method = "fiml" is only available for categorical g-wise distances.')
  }
  if (method == "complete" && any(!is.finite(x_mat))) {
    stop("the complete-data g-wise estimator requires complete, finite ratings.")
  }

  if (categorical) {
    cats <- sort(unique(c(x_mat)))
    x_indexed <- matrix(match(x_mat, cats) - 1L, nrow = nrow(x_mat), ncol = ncol(x_mat))
    storage.mode(x_indexed) <- "integer"
    out <- rcpp_kappa_gwise_categorical(
      x = x_indexed,
      distance_type = distance,
      method = method,
      g = as.integer(g),
      max_chance_tuples = as.integer(max_chance_tuples),
      em_options = em_options
    )
  } else {
    storage.mode(x_mat) <- "double"
    out <- rcpp_kappa_gwise_continuous(
      x = x_mat,
      distance_type = distance,
      method = method,
      g = as.integer(g),
      max_chance_tuples = as.integer(max_chance_tuples)
    )
  }

  .new_estimate(out, c("Conger", "Fleiss"),
                method = paste0("gwise_", method), weight = distance,
                distance = distance, g = as.integer(g))
}
