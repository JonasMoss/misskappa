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

.alpha_from_cpp <- function(out, method, em_options = list()) {
  estimates <- as.numeric(out$estimates)
  names(estimates) <- "alpha"
  vcov_mat <- out$vcov
  dimnames(vcov_mat) <- list(names(estimates), names(estimates))
  psi_mat <- out$psi
  if (prod(dim(psi_mat)) > 0L) colnames(psi_mat) <- names(estimates)
  .warn_null_frac(out$null_frac, em_options)

  structure(
    list(
      estimates = estimates,
      vcov = vcov_mat,
      psi = psi_mat,
      null_frac = out$null_frac,
      method = method,
      weight = "score"
    ),
    class = "misskappa_estimate"
  )
}

#' Coefficient alpha with missing data
#'
#' @description
#' Estimates Cronbach's coefficient alpha for a subjects-by-items matrix of
#' item scores with missing entries. Pick the estimator with `estimator`:
#' \itemize{
#'   \item `"pairwise"` --- pairwise-available covariance moments under MCAR,
#'     with an influence-function sandwich standard error;
#'   \item `"cat_fiml"` --- saturated-multinomial full-information maximum
#'     likelihood for finite-category items under ignorable missingness;
#'   \item `"nt_fiml"` --- robust normal-theory FIML (saturated Gaussian
#'     covariance by EM, sandwich delta-method SE) for continuous items under
#'     ignorable missingness.
#' }
#'
#' For ordinal or otherwise categorical responses under `"pairwise"` or
#' `"nt_fiml"`, either pass the scored numeric matrix directly or pass integer
#' category codes with `values`; the observed categories are sorted and mapped
#' to `values` before estimation.
#'
#' All three estimators target fixed-item coefficient alpha. With missing
#' entries, every item must be observed at least once and every item pair must
#' be jointly observed by at least one subject. If that complete pairwise
#' co-observation condition fails, alpha is not identified from the observed
#' missing-data pattern and the fit errors before running the estimator.
#'
#' @param x A subjects-by-items numeric matrix or data frame; `NA` and other
#'   non-finite values indicate missing entries. For `estimator = "cat_fiml"`,
#'   entries are integer category codes.
#' @param estimator One of `"pairwise"`, `"cat_fiml"`, or `"nt_fiml"`.
#' @param values Optional numeric vector of scores for observed categories.
#'   When supplied, its length must equal the number of unique finite entries
#'   in `x`.
#' @param em_options Named list tuning the EM fit for the likelihood
#'   estimators. For `"nt_fiml"`: `tol`, `max_iter`; `fd_h` is accepted for
#'   backward compatibility and ignored. For `"cat_fiml"`: `tol`, `max_iter`,
#'   `prune_tol`, `start_alpha`, `info_rcond` (the relative
#'   eigenvalue cutoff used when inverting Louis' observed information), and
#'   `flatten` (total Dirichlet pseudo-mass spread over the complete pattern
#'   table; any positive value makes the fitted table the unique interior
#'   posterior mode when the saturated likelihood is flat, at the cost of
#'   shrinking it toward uniform with weight `flatten / (n + flatten)`;
#'   `0`, the default, is strict ML). Flattening is a uniqueness device,
#'   not an inference upgrade: it leaves the point estimate essentially
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
#' @examplesIf requireNamespace("psych", quietly = TRUE)
#' # Real item-level missing data: the Neuroticism scale of psych::bfi
#' # (2800 respondents, ~360 with at least one missing item). The FIML
#' # estimator is valid under ignorable missingness.
#' data(bfi, package = "psych")
#' N <- paste0("N", 1:5)
#' alpha(bfi[, N], estimator = "nt_fiml")
#'
#' @export
alpha <- function(x,
                  estimator = c("pairwise", "cat_fiml", "nt_fiml"),
                  values = NULL,
                  em_options = list()) {
  estimator <- match.arg(estimator)
  .check_pattern_identifiable(.pattern_observed(x), unit = "item",
                              coefficient = "coefficient alpha")

  if (estimator == "pairwise") {
    X <- .alpha_score_matrix(x, values)
    out <- rcpp_alpha_available_continuous(X)
    return(.alpha_from_cpp(out, method = "pairwise"))
  }

  if (estimator == "nt_fiml") {
    X <- .alpha_score_matrix(x, values)
    fit <- alpha_continuous(X, em_options = em_options)
    fit$method <- "nt_fiml"
    return(fit)
  }

  fit <- alpha_cat_fiml(x, values = values, em_options = em_options)
  fit$method <- "cat_fiml"
  fit
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
#' wrapper checks this before dispatch; direct callers reach the C++ guard.
#'
#' @param x A subjects-by-items matrix of integer category codes; `NA`s
#'   indicate missing entries.
#' @param values Optional numeric vector of category scores. Defaults to the
#'   sorted unique observed categories.
#' @param em_options Named list of options for the categorical EM fit:
#'   `tol`, `max_iter`, `prune_tol`, `start_alpha`, `info_rcond`, `flatten`.
#'
#' @return An object of class `misskappa_estimate` carrying one coefficient
#'   named `alpha` and its asymptotic covariance matrix.
#'
#' @keywords internal
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
  estimates <- as.numeric(out$estimates)
  names(estimates) <- c("Conger", "Fleiss", "Brennan-Prediger")
  vcov_mat <- out$vcov
  dimnames(vcov_mat) <- list(names(estimates), names(estimates))
  psi_mat <- out$psi
  if (prod(dim(psi_mat)) > 0L) colnames(psi_mat) <- names(estimates)
  null_frac <- out$null_frac
  if (length(null_frac) == length(estimates)) names(null_frac) <- names(estimates)
  .warn_null_frac(null_frac, em_options)
  structure(list(estimates = estimates, vcov = vcov_mat, psi = psi_mat,
                 null_frac = null_frac, method = method, weight = weight),
            class = "misskappa_estimate")
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
#' rater pair to be jointly observed by at least one subject; otherwise the
#' requested saturated pairwise agreement functional is not identified.
#'
#' @param x A subjects-by-raters matrix of integer category codes; `NA`s
#'   indicate missing entries.
#' @param weights Character vector containing any of `"identity"`,
#'   `"nominal"`, `"linear"`, or `"quadratic"`. `"nominal"` is an alias for
#'   `"identity"`.
#' @param values Optional numeric category scores for metric weights.
#' @param em_options Named list of categorical EM options.
#'
#' @return A named list of `misskappa_estimate` objects.
#'
#' @keywords internal
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
    estimates <- as.numeric(raw[[i]]$estimates)
    names(estimates) <- c("Conger", "Fleiss", "Brennan-Prediger")
    vcov_mat <- raw[[i]]$vcov
    dimnames(vcov_mat) <- list(names(estimates), names(estimates))
    psi_mat <- raw[[i]]$psi
    if (prod(dim(psi_mat)) > 0L) colnames(psi_mat) <- names(estimates)
    null_frac <- raw[[i]]$null_frac
    if (length(null_frac) == length(estimates)) names(null_frac) <- names(estimates)
    if (i == 1L) .warn_null_frac(null_frac, em_options)
    structure(
      list(estimates = estimates, vcov = vcov_mat, psi = psi_mat,
           null_frac = null_frac, method = "fiml", weight = weights[[i]]),
      class = "misskappa_estimate"
    )
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
#'   \item `"ipw"` --- inverse-probability-weighted moment estimator for
#'     categorical ratings under MCAR, consistent under non-exchangeable raters
#'     with rater-varying observation rates;
#'   \item `"cat_fiml"` --- saturated-multinomial full-information maximum
#'     likelihood for categorical ratings under ignorable missingness;
#'   \item `"pairwise"` --- pairwise-available moment estimator for the
#'     quadratically weighted (scored) coefficient under MCAR;
#'   \item `"nt_fiml"` --- robust normal-theory FIML for the quadratically
#'     weighted coefficient under ignorable missingness.
#' }
#'
#' `"ipw"` and `"cat_fiml"` take integer category codes and the nominal,
#' linear, or quadratic `weight`.
#' `"pairwise"` and `"nt_fiml"` are quadratic by construction: they treat `x`
#' as numeric scores (pass scores directly, or integer codes with `values`)
#' and require `weight = "quadratic"`.
#'
#' Two generalizations are reached through the same call. Set `g > 2` for the
#' Frechet / Hubert g-wise (multirater) kernels, which report the Conger-type
#' (distinct rater combinations) and Fleiss-type coefficients. Pass a
#' subjects-by-raters-by-features array to score vector-valued ratings with a
#' component-separable loss. See the arguments for the estimator / weight
#' combinations supported in each mode.
#'
#' For scalar fixed-rater coefficients with missing entries, the observed-data
#' pattern must identify the requested functional. For ordinary pairwise kappa
#' (`g = 2`, and all quadratic requests), every rater must be observed at least
#' once and every rater pair must be jointly observed by at least one subject.
#' For non-quadratic `g > 2`, every requested rater g-tuple must be jointly
#' observed by at least one subject. Violations are reported as
#' non-identification before the estimator runs.
#'
#' @param x A subjects-by-raters matrix or data frame of scalar ratings, or a
#'   subjects-by-raters-by-features 3-D array of vector-valued ratings (the
#'   array form is auto-detected). Integer category codes for `"ipw"` /
#'   `"cat_fiml"`; numeric scores for `"pairwise"` / `"nt_fiml"`. `NA` marks
#'   missing entries.
#' @param estimator One of `"ipw"`, `"cat_fiml"`, `"pairwise"`, `"nt_fiml"`.
#'   Vector-valued ratings support `"pairwise"` and `"ipw"` only; continuous
#'   g-wise kernels (`weight = "linear"`, `g > 2`) support `"pairwise"` and
#'   `"ipw"` only.
#' @param weight Weighting scheme for `"ipw"` / `"cat_fiml"`: `"nominal"` (the
#'   default), `"linear"`, or `"quadratic"`. `"pairwise"` / `"nt_fiml"` require
#'   `"quadratic"`. For `g > 2` only `"nominal"`, `"linear"`, `"hubert"`, and
#'   `"quadratic"` are defined, where `"hubert"` is the all-raters-equal
#'   multirater kernel (g > 2 only). For vector-valued ratings `weight` selects
#'   the component loss (`"nominal"` -> Hamming, `"linear"` -> L1,
#'   `"quadratic"` -> squared).
#' @param g Arity of the multirater disagreement kernel. `g = 2` (the default)
#'   is ordinary pairwise kappa; `g > 2` uses the Frechet / Hubert g-wise
#'   family. `weight = "quadratic"` is g-invariant, so `g` is ignored there and
#'   the cheap closed form is used. Ignored for vector-valued ratings.
#' @param values Optional numeric vector of category scores used by the metric
#'   weightings (and by the quadratic / continuous g-wise estimators to map
#'   integer codes to scores). Defaults to the sorted unique observed
#'   categories.
#' @param em_options Named list of EM options for the likelihood estimators
#'   (`"cat_fiml"`, `"nt_fiml"`). For `"cat_fiml"`: `tol`, `max_iter`,
#'   `prune_tol`, `start_alpha`, `info_rcond`, and `flatten` (total Dirichlet
#'   pseudo-mass spread over the complete pattern table; any positive value
#'   selects the unique interior posterior mode when the saturated likelihood
#'   is flat, shrinking the fitted table toward uniform with weight
#'   `flatten / (n + flatten)`; `0`, the default, is strict ML). Flattening
#'   is a uniqueness device, not an inference upgrade: it leaves the point
#'   estimate essentially unchanged but makes the reported standard errors
#'   conservative (roughly 50\% too wide in calibration). Pass any
#'   subset.
#' @param ... Mode-specific extras: `feature_weights` and `loss` (e.g.
#'   `"rms"`) for vector-valued ratings, and `max_chance_tuples` (cap on the
#'   number of chance tuples enumerated) for `g > 2`.
#'
#' @details
#' The public R API deliberately exposes only the main weighting schemes used in
#' current agreement work: `"nominal"`, `"linear"`, and `"quadratic"` (plus the
#' `"hubert"` g-wise kernel for `g > 2`). Older Krippendorff / irrCAC-style
#' categorical schemes (`"ordinal"`, `"radical"`, `"ratio"`, `"circular"`,
#' `"bipolar"`) remain in the internal C++ layer and Rcpp glue for validation
#' and parity work, but are unsupported and not part of the recommended R
#' interface.
#'
#' @references
#' Krippendorff, K. (2011). Computing Krippendorff's alpha-reliability.
#'
#' Gwet, K. L. (2019). irrCAC: Computing chance-corrected agreement
#' coefficients among raters.
#'
#' @return An object of class `misskappa_estimate` carrying the named
#'   coefficient estimates and their asymptotic covariance matrix. Methods:
#'   `print`, `coef`, `vcov`, `confint`, `as.data.frame`, and
#'   `stats::influence`. The object also carries a `psi` component --- the
#'   n-by-K matrix of per-subject influence functions, satisfying
#'   `vcov == crossprod(psi) / n^2` --- for power users who want to build
#'   their own contrasts or joint tests.
#'
#' @examples
#' # Categorical ratings with missing entries (Gwet 2014): the
#' # inverse-probability-weighted estimator under MCAR, with linear weights.
#' kappa(dat.gwet2014, estimator = "ipw", weight = "linear")
#'
#' # Scored coefficient (quadratic loss) via the pairwise-available moment
#' # estimator; treats the columns as numeric scores.
#' kappa(dat.zapf2016, estimator = "pairwise")
#'
#' # g-wise (multirater) nominal kappa over triples of raters.
#' kappa(dat.gwet2014, estimator = "ipw", weight = "nominal", g = 3)
#'
#' # Vector-valued ratings: a subjects-by-raters-by-features array. Each of the
#' # four expert observers classifies crackles at six chest sites, so a rating
#' # is a six-vector; the default Hamming loss counts disagreeing sites.
#' kappa(dat.vanbelle2019[, 1:4, ], estimator = "ipw")
#'
#' @export
kappa <- function(x,
                  estimator = c("ipw", "cat_fiml", "pairwise", "nt_fiml"),
                  weight = c("nominal", "linear", "quadratic", "hubert"),
                  g = 2L,
                  values = NULL,
                  em_options = list(),
                  ...) {
  estimator <- match.arg(estimator)
  weight_supplied <- !missing(weight)
  weight <- match.arg(weight)
  dots <- list(...)

  ## Vector-valued ratings: subjects-by-raters-by-features array.
  if (is.array(x) && length(dim(x)) == 3L) {
    return(.kappa_vector_dispatch(x, estimator, weight, dots))
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
  ## g collapses to the same coefficient, so always take the cheap closed form
  ## (the g = 2 path) and never the combinatorial g-wise enumeration.
  if (weight == "quadratic" || g == 2L) {
    ## Every public scalar g = 2 path surfaces a fixed-rater coefficient
    ## (Conger, and Brennan-Prediger for the raw estimators), so the saturated
    ## functional needs a complete rater co-observation graph.
    .check_pattern_identifiable(
      .pattern_observed(x), unit = "rater",
      coefficient = "Conger's kappa (and Brennan-Prediger)",
      arity = 2L)
    return(.kappa_scalar_g2(x, estimator, weight, weight_supplied, values,
                            em_options))
  }

  ## g > 2, non-quadratic: the Frechet / Hubert g-wise family.
  .check_pattern_identifiable(
    .pattern_observed(x), unit = "rater",
    coefficient = sprintf("%s-wise kappa", g),
    arity = g)
  .kappa_scalar_gwise(x, estimator, weight, g, values, em_options, dots)
}

# Ordinary pairwise (g = 2) kappa: the historical kappa() body, also reached by
# any quadratic request because quadratic kappa does not depend on g.
.kappa_scalar_g2 <- function(x, estimator, weight, weight_supplied, values,
                             em_options) {
  if (estimator %in% c("pairwise", "nt_fiml")) {
    if (weight_supplied && weight != "quadratic") {
      stop(sprintf('estimator = "%s" requires weight = "quadratic".', estimator))
    }
    X <- .score_matrix(x, values)
    fit <- kappa_continuous(
      X,
      method = if (estimator == "pairwise") "available" else "fiml",
      weight = "quadratic",
      em_options = em_options
    )
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
    pairwise = "complete",
    cat_fiml = if (map$family == "categorical") "fiml" else stop(
      'estimator = "cat_fiml" has no continuous g-wise variant; use ',
      'weight = "nominal"/"hubert", or estimator = "ipw".'),
    nt_fiml  = stop('estimator = "nt_fiml" is only available for ',
                    'weight = "quadratic" (where g is ignored).')
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

# Vector-valued ratings (subjects-by-raters-by-features array). Maps the public
# weight to a component loss and delegates to kappa_vector(); `g` is ignored.
.kappa_vector_dispatch <- function(x, estimator, weight, dots) {
  if (!estimator %in% c("pairwise", "ipw")) {
    stop(sprintf(
      'vector-valued ratings support estimator = "pairwise" or "ipw"; got "%s".',
      estimator))
  }
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
  fit$weight <- loss
  fit
}

#' @export
print.misskappa_estimate <- function(x, digits = 4, level = 0.95,
                                     transform = c("none", "fisher"), ...) {
  transform <- match.arg(transform)
  hdr <- sprintf("misskappa: estimator=%s, weight=%s", x$method, x$weight)
  if (!is.null(x$g) && x$g > 2L) hdr <- sprintf("%s, g=%d", hdr, x$g)
  cat(hdr, "\n", sep = "")
  se <- sqrt(diag(x$vcov))
  ci <- stats::confint(x, level = level, transform = transform)
  tab <- data.frame(estimate = x$estimates, se = se,
                    lower = ci[, 1L], upper = ci[, 2L],
                    row.names = names(x$estimates))
  print(round(tab, digits = digits), ...)
  invisible(x)
}

#' @export
coef.misskappa_estimate <- function(object, ...) object$estimates

#' @export
vcov.misskappa_estimate <- function(object, ...) object$vcov

#' Confidence intervals for misskappa coefficients
#'
#' @description
#' Wald confidence intervals for the coefficients of a `misskappa_estimate`,
#' on the natural scale (`transform = "none"`) or via the variance-stabilising
#' Fisher z transform (`transform = "fisher"`). The Fisher interval is built on
#' the `atanh` scale with a delta-method standard error
#' \eqn{\mathrm{se}_z = \mathrm{se} / (1 - \hat\theta^2)} and back-transformed
#' with `tanh`, so it always lies in \eqn{(-1, 1)} and tends to have better
#' small-sample coverage near the upper boundary.
#'
#' @param object A `misskappa_estimate` object.
#' @param parm Optional subset of coefficients (names or indices); defaults to
#'   all.
#' @param level Confidence level.
#' @param transform Either `"none"` (natural-scale Wald interval) or `"fisher"`
#'   (delta-method interval on the `atanh` scale, back-transformed with
#'   `tanh`). Coefficients with \eqn{|\hat\theta| \ge 1} yield `NA` limits
#'   under `"fisher"`.
#' @param ... Unused; present for S3 generic conformance.
#'
#' @return A two-column numeric matrix of lower and upper limits, one row per
#'   coefficient.
#'
#' @examples
#' fit <- kappa(dat.gwet2014, estimator = "ipw")
#' confint(fit)                       # natural-scale Wald interval
#' confint(fit, transform = "fisher") # Fisher z interval; always within (-1, 1)
#'
#' @export
confint.misskappa_estimate <- function(object, parm = NULL, level = 0.95,
                                       transform = c("none", "fisher"), ...) {
  transform <- match.arg(transform)
  est <- object$estimates
  se <- sqrt(diag(object$vcov))
  z <- stats::qnorm((1 + level) / 2)

  if (transform == "fisher") {
    interior <- abs(est) < 1
    if (!all(interior)) {
      warning("Fisher transform requires |estimate| < 1; ",
              "returning NA limits for boundary coefficients.")
    }
    g <- ifelse(interior, atanh(ifelse(interior, est, 0)), NA_real_)
    se_z <- ifelse(interior, se / (1 - est^2), NA_real_)
    lo <- tanh(g - z * se_z)
    hi <- tanh(g + z * se_z)
  } else {
    lo <- est - z * se
    hi <- est + z * se
  }

  ci <- cbind(lo, hi)
  colnames(ci) <- c(sprintf("%.1f %%", 100 * (1 - level) / 2),
                    sprintf("%.1f %%", 100 * (1 + level) / 2))
  rownames(ci) <- names(est)
  if (!is.null(parm)) ci <- ci[parm, , drop = FALSE]
  ci
}

#' @export
as.data.frame.misskappa_estimate <- function(x, row.names = NULL,
                                             optional = FALSE, ...) {
  data.frame(
    coefficient = names(x$estimates),
    estimate = unname(x$estimates),
    se = sqrt(diag(x$vcov)),
    row.names = row.names,
    stringsAsFactors = FALSE
  )
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
#' agree under the usual asymptotics; the closed form can be faster on
#' large `n`.
#'
#' With missing ratings, every rater must be observed at least once and every
#' rater pair must be jointly observed by at least one subject.
#'
#' @param x A subjects-by-raters numeric matrix of category scores; `NA`
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
#' @keywords internal
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
  estimates <- as.numeric(out$estimates)
  names(estimates) <- c("Conger", "Fleiss", "Brennan-Prediger")
  vcov_mat <- out$vcov
  dimnames(vcov_mat) <- list(names(estimates), names(estimates))
  psi_mat <- out$psi
  if (prod(dim(psi_mat)) > 0L) colnames(psi_mat) <- names(estimates)
  structure(
    list(
      estimates = estimates,
      vcov = vcov_mat,
      psi = psi_mat,
      method = "quadratic",
      weight = "quadratic"
    ),
    class = "misskappa_estimate"
  )
}

#' Weighted agreement coefficients for counts-format input
#'
#' @description
#' Estimates Fleiss and Brennan-Prediger weighted agreement coefficients
#' for counts-format data: each row is one subject; each column is a
#' category; `x[i, k]` is the number of raters who assigned subject `i`
#' to category `k`. Row sums (number of raters per subject) need not be
#' uniform across subjects.
#'
#' The count-format moment estimator follows the Fleiss-Cuzick unequal-judges
#' convention: observed row disagreement is weighted by `r_i - 1`, and chance
#' disagreement uses the pooled rating-token margin. This differs from the
#' unit-weighted distribution/count convention used by some software; the
#' unit-weighted comparator is kept only in the C++ API for validation studies.
#'
#' Rows with fewer than two observed ratings contain no within-subject pair
#' information and receive zero observed-disagreement weight; at least one row
#' with two or more ratings is required.
#'
#' Counts data assumes **exchangeable iid raters** drawn from a single
#' category distribution. For non-exchangeable raters, use rater-
#' identified data and `kappa(x, estimator = "cat_fiml")` instead.
#'
#' Conger is not reported (raters are not identified in this input format),
#' and neither IPW nor Gwet are meaningful (per-rater observation rates are
#' aggregated away by the counts representation).
#'
#' @param x A subjects-by-categories non-negative integer matrix.
#' @param estimator Either `"fleiss_cuzick"` (the Fleiss-Cuzick count-format
#'   moment estimator) or `"cat_fiml"` (EM over the composition simplex with
#'   multivariate hypergeometric weights for completing partial counts). The
#'   two agree exactly when every row of `x` sums to `r_total`; with partial
#'   counts (some `r_i < r_total`) the `"cat_fiml"` estimator can be more
#'   efficient.
#' @param weight Weighting scheme: `"nominal"` (default), `"linear"`, or
#'   `"quadratic"`.
#' @param values Optional length-C numeric vector of category scores used
#'   by the metric weightings. Defaults to `1:C`.
#' @param r_total Total number of raters per subject. Defaults to the
#'   maximum observed row sum. Required for `estimator = "cat_fiml"` when rows
#'   have varying totals.
#' @param em_options Named list of EM options for `estimator = "cat_fiml"`:
#'   `tol`, `max_iter`, `prune_tol`, `start_alpha`, `info_rcond`. The
#'   raw-data `flatten` option does not apply here: the count-data FIML fits
#'   a low-dimensional composition model, not the saturated joint.
#'
#' @details
#' Older Krippendorff / irrCAC-style categorical schemes (`"ordinal"`,
#' `"radical"`, `"ratio"`, `"circular"`, `"bipolar"`) remain available only
#' through unsupported internal helpers and the C++ `misskappa::loss` factories;
#' they are retained for validation and parity work, not as part of the
#' recommended R interface.
#'
#' @references
#' Krippendorff, K. (2011). Computing Krippendorff's alpha-reliability.
#'
#' Gwet, K. L. (2019). irrCAC: Computing chance-corrected agreement
#' coefficients among raters.
#'
#' @return A `misskappa_estimate` object with `Fleiss` and
#'   `Brennan-Prediger` coefficients, the 2x2 vcov, a `psi` component, and
#'   the registered `stats::influence` method.
#'
#' @examples
#' # Fleiss (1971) psychiatric diagnoses: 30 subjects, 6 raters, 5 categories,
#' # in counts format (one row per subject, one column per category).
#' kappa_counts(dat.fleiss1971, estimator = "fleiss_cuzick")
#'
#' # Treat the categories as ordered scores with a quadratic loss.
#' kappa_counts(dat.fleiss1971, estimator = "fleiss_cuzick", weight = "quadratic")
#'
#' @export
kappa_counts <- function(x,
                         estimator = c("fleiss_cuzick", "cat_fiml"),
                         weight = c("nominal", "linear", "quadratic"),
                         values = NULL,
                         r_total = NULL,
                         em_options = list()) {
  # Backward-compatible local alias; not documented as a public option.
  if (identical(estimator, "pairwise")) estimator <- "fleiss_cuzick"
  estimator <- match.arg(estimator)
  weight <- match.arg(weight)
  cpp_weight <- if (weight == "nominal") "identity" else weight
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a matrix or data frame.")
  }
  x_mat <- as.matrix(x)
  if (!is.numeric(x_mat)) stop("'x' must be numeric.")
  storage.mode(x_mat) <- "integer"

  if (estimator == "fleiss_cuzick") {
    out <- rcpp_kappa_counts(x = x_mat, weight_type = cpp_weight, values = values)
  } else {
    if (is.null(r_total)) r_total <- as.integer(max(rowSums(x_mat)))
    out <- rcpp_kappa_fiml_counts(
      x = x_mat,
      weight_type = cpp_weight,
      values = values,
      r_total = as.integer(r_total),
      em_options = em_options
    )
  }

  estimates <- as.numeric(out$estimates)
  names(estimates) <- c("Fleiss", "Brennan-Prediger")
  vcov_mat <- out$vcov
  dimnames(vcov_mat) <- list(names(estimates), names(estimates))
  psi_mat <- out$psi
  if (prod(dim(psi_mat)) > 0L) colnames(psi_mat) <- names(estimates)

  structure(
    list(
      estimates = estimates,
      vcov = vcov_mat,
      psi = psi_mat,
      method = estimator,
      weight = weight
    ),
    class = "misskappa_estimate"
  )
}

#' Weighted agreement coefficients for continuous ratings
#'
#' @description
#' Estimates Conger / Fleiss weighted agreement coefficients for raw
#' real-valued ratings under missing data. The MCAR moment routes
#' (`"available"`, `"ipw"`, `"gwet"`) cover any continuous loss kernel; the
#' normal-FIML route (`"fiml"`, quadratic kernel only) fits a saturated
#' multivariate-normal covariance by EM and is valid under ignorable (MCAR or
#' MAR) missingness. Missing entries are encoded as `NA` (or any non-finite
#' value).
#'
#' Brennan-Prediger is not reported for continuous data; the chance-
#' disagreement baseline requires a finite number of categories.
#'
#' For all missing-data routes here, every rater must be observed at least once
#' and every rater pair must be jointly observed by at least one subject.
#'
#' @param x A subjects-by-raters numeric matrix.
#' @param method One of `"available"`, `"ipw"`, `"gwet"`, or `"fiml"`.
#'   `"fiml"` requires `weight = "quadratic"` and dispatches to
#'   [kappa_quadratic_fiml()].
#' @param weight Continuous loss kernel: `"identity"` (binary 0/1),
#'   `"linear"`, `"quadratic"`, `"radical"`, or `"ratio"`. The data range
#'   `[min(x), max(x)]` is used to parameterise the kernel. Only
#'   `"quadratic"` is available with `method = "fiml"`.
#' @param em_options Used only when `method = "fiml"`: named list with `tol`
#'   and `max_iter`. `fd_h` is accepted for backward compatibility and ignored.
#'
#' @return A `misskappa_estimate` object carrying named coefficients
#'   (`Conger`, `Fleiss`) and the 2x2 asymptotic covariance matrix.
#'
#' @keywords internal
kappa_continuous <- function(x,
                             method = c("available", "ipw", "gwet", "fiml"),
                             weight = c("quadratic", "linear", "identity",
                                        "unweighted", "radical", "ratio"),
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
           "Conger and Fleiss kappas are functions of the mean and covariance ",
           "only under quadratic weighting.")
    }
    return(kappa_quadratic_fiml(x_mat, em_options = em_options))
  }

  out <- rcpp_kappa_continuous(
    x = x_mat,
    method = method,
    weight_type = weight
  )

  estimates <- as.numeric(out$estimates)
  names(estimates) <- c("Conger", "Fleiss")
  vcov_mat <- out$vcov
  dimnames(vcov_mat) <- list(names(estimates), names(estimates))
  psi_mat <- out$psi
  if (prod(dim(psi_mat)) > 0L) colnames(psi_mat) <- names(estimates)

  structure(
    list(
      estimates = estimates,
      vcov = vcov_mat,
      psi = psi_mat,
      method = method,
      weight = weight
    ),
    class = "misskappa_estimate"
  )
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
#' @keywords internal
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

  estimates <- as.numeric(out$estimates)
  names(estimates) <- c("Conger", "Fleiss")
  vcov_mat <- out$vcov
  dimnames(vcov_mat) <- list(names(estimates), names(estimates))
  psi_mat <- out$psi
  if (prod(dim(psi_mat)) > 0L) colnames(psi_mat) <- names(estimates)

  structure(
    list(
      estimates = estimates,
      vcov = vcov_mat,
      psi = psi_mat,
      method = paste0("vector_", method),
      weight = loss,
      loss = loss,
      feature_weights = as.numeric(feature_weights)
    ),
    class = "misskappa_estimate"
  )
}

#' G-wise agreement coefficients for rectangular ratings
#'
#' @description
#' Estimates g-wise Conger-type and Fleiss-type agreement coefficients using
#' multirater disagreement kernels (the Frechet / Hubert family). The
#' `"complete"` estimator requires every entry observed; `"ipw"` and `"fiml"`
#' admit missing entries (FIML is categorical only) when every requested rater
#' g-tuple is jointly observed by at least one subject.
#'
#' @param x A subjects-by-raters matrix or data frame.
#' @param distance Multirater disagreement kernel. `"nominal"` uses Frechet
#'   mode disagreement for categorical ratings; `"absolute"` uses median
#'   absolute deviation; `"quadratic"` uses mean squared deviation; `"hubert"`
#'   uses all-raters-equal disagreement.
#' @param method `"complete"` (complete-data), `"ipw"` (MCAR inverse-probability
#'   weighting), or `"fiml"` (categorical full-information ML under ignorable
#'   missingness; not defined for the continuous kernels).
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
#' @keywords internal
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

  estimates <- as.numeric(out$estimates)
  names(estimates) <- c("Conger", "Fleiss")
  vcov_mat <- out$vcov
  dimnames(vcov_mat) <- list(names(estimates), names(estimates))
  psi_mat <- out$psi
  if (prod(dim(psi_mat)) > 0L) colnames(psi_mat) <- names(estimates)

  structure(
    list(
      estimates = estimates,
      vcov = vcov_mat,
      psi = psi_mat,
      method = paste0("gwise_", method),
      weight = distance,
      distance = distance,
      g = as.integer(g)
    ),
    class = "misskappa_estimate"
  )
}
