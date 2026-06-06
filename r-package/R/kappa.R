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

.alpha_from_cpp <- function(out, method) {
  estimates <- as.numeric(out$estimates)
  names(estimates) <- "alpha"
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
      weight = "score"
    ),
    class = "misskappa_estimate"
  )
}

#' Coefficient alpha with missing data
#'
#' @description
#' Estimates Cronbach's coefficient alpha for a subjects-by-items matrix of
#' item scores with missing entries. Missingness is handled either by
#' pairwise-available covariance moments under MCAR (`method = "available"`)
#' or by saturated full-information maximum likelihood under ignorable
#' missingness (`method = "fiml"`). For FIML, `type` selects the response
#' model: `"normal"` fits a saturated multivariate-normal covariance by EM
#' (continuous items), and `"categorical"` fits the saturated multinomial
#' full-response distribution by EM and maps it to the implied scored
#' covariance (finite-category items).
#'
#' For ordinal or otherwise categorical responses under `"available"` or
#' `type = "normal"`, either pass the scored numeric matrix directly or pass
#' integer category codes with `values`; the observed categories are sorted
#' and mapped to `values` before estimation.
#'
#' @param x A subjects-by-items numeric matrix or data frame; `NA` and other
#'   non-finite values indicate missing entries. For
#'   `method = "fiml", type = "categorical"`, entries are integer category
#'   codes.
#' @param method One of `"available"` (pairwise-available covariance moments
#'   with an influence-function sandwich SE) or `"fiml"` (saturated EM under
#'   ignorable missingness).
#' @param type For `method = "fiml"`, the response model: `"normal"` (saturated
#'   Gaussian covariance via EM with a sandwich delta-method SE) or
#'   `"categorical"` (saturated multinomial EM mapped to the scored
#'   covariance). Ignored when `method = "available"`.
#' @param values Optional numeric vector of scores for observed categories.
#'   When supplied, its length must equal the number of unique finite entries
#'   in `x`.
#' @param em_options Named list tuning the EM fit. For `type = "normal"`:
#'   `tol`, `max_iter`, `fd_h`. For `type = "categorical"`: `tol`, `max_iter`,
#'   `prune_tol`, `start_alpha`, `info_rcond` (the relative eigenvalue cutoff
#'   used when inverting Louis' observed information). Pass any subset.
#'
#' @return An object of class `misskappa_estimate` carrying one coefficient
#'   named `alpha` and its asymptotic covariance matrix. Methods: `print`,
#'   `coef`, `vcov`, `confint`, `as.data.frame`, and, when available,
#'   `influence`.
#'
#' @examples
#' set.seed(1)
#' n <- 400L; p <- 5L
#' L <- chol(0.3 + 0.7 * diag(p))
#' x <- matrix(rnorm(n * p), n, p) %*% L
#' x[matrix(runif(n * p) < 0.15, n, p)] <- NA
#' fit <- alpha(x, method = "fiml")
#' coef(fit)
#' confint(fit)
#'
#' @export
alpha <- function(x,
                  method = c("available", "fiml"),
                  type = c("normal", "categorical"),
                  values = NULL,
                  em_options = list()) {
  method <- match.arg(method)
  type <- match.arg(type)

  if (method == "available") {
    X <- .alpha_score_matrix(x, values)
    out <- rcpp_alpha_available_continuous(X)
    return(.alpha_from_cpp(out, method = "alpha-available"))
  }

  if (type == "normal") {
    X <- .alpha_score_matrix(x, values)
    return(alpha_continuous(X, em_options = em_options))
  }

  alpha_cat_fiml(x, values = values, em_options = em_options)
}

#' Saturated categorical FIML coefficient alpha (internal)
#'
#' @description
#' Backend for `alpha(method = "fiml", type = "categorical")`. Fits the
#' saturated multinomial full-response distribution with EM, then maps that
#' distribution to the implied scored covariance matrix.
#'
#' @param x A subjects-by-items matrix of integer category codes; `NA`s
#'   indicate missing entries.
#' @param values Optional numeric vector of category scores. Defaults to the
#'   sorted unique observed categories.
#' @param em_options Named list of options for the categorical EM fit:
#'   `tol`, `max_iter`, `prune_tol`, `start_alpha`, `info_rcond`.
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

  .alpha_from_cpp(out, method = "alpha-cat-fiml")
}

#' Weighted agreement coefficients with missing data
#'
#' @description
#' Estimates Cohen / Fleiss / Conger / Brennan-Prediger weighted agreement
#' coefficients for raw categorical ratings, supporting missing entries
#' under MCAR (`"available"`, `"ipw"`, `"gwet"`) or MAR (`"fiml"`).
#'
#' @param x A subjects-by-raters matrix of integer category codes; `NA`s
#'   indicate missing entries.
#' @param method One of `"available"`, `"ipw"`, `"fiml"`, `"gwet"`.
#' @param weight Weighting scheme: `"identity"` (the default; equivalent to
#'   `"unweighted"`), `"linear"`, `"quadratic"`, `"ordinal"`, `"radical"`,
#'   `"ratio"`, `"circular"`, or `"bipolar"`.
#' @param values Optional numeric vector of category scores used by the
#'   metric weightings. Defaults to the sorted unique observed categories.
#' @param em_options Named list of options for `method = "fiml"`:
#'   `tol`, `max_iter`, `prune_tol`, `start_alpha`, `info_rcond`.
#'   `info_rcond` is the relative eigenvalue cutoff used when inverting
#'   Louis' observed information. Pass any subset.
#'
#' @return An object of class `misskappa_estimate` carrying the named
#'   coefficient estimates and the 3x3 asymptotic covariance matrix.
#'   Methods: `print`, `coef`, `vcov`, `confint`, `as.data.frame`.
#'
#' @export
kappa <- function(x,
                  method = c("available", "ipw", "fiml", "gwet"),
                  weight = c("identity", "unweighted", "linear", "quadratic",
                             "ordinal", "radical", "ratio", "circular", "bipolar"),
                  values = NULL,
                  em_options = list()) {
  method <- match.arg(method)
  weight <- match.arg(weight)

  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a matrix or data frame.")
  }
  x_mat <- as.matrix(x)
  if (!is.numeric(x_mat)) stop("'x' must be numeric.")
  storage.mode(x_mat) <- "integer"

  out <- rcpp_kappa_raw(
    x = x_mat,
    method = method,
    weight_type = weight,
    values = values,
    em_options = em_options
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
      method = method,
      weight = weight
    ),
    class = "misskappa_estimate"
  )
}

#' @export
print.misskappa_estimate <- function(x, digits = 4, level = 0.95,
                                     transform = c("none", "fisher"), ...) {
  transform <- match.arg(transform)
  cat(sprintf("misskappa: method=%s, weight=%s\n", x$method, x$weight))
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

#' Per-subject influence-function matrix
#'
#' @description
#' Returns the `n x K` matrix of per-subject influence functions for the
#' coefficient estimates, where `n` is the number of subjects and `K` is
#' the number of coefficients. Estimators that do not expose influence
#' functions return `NULL`.
#'
#' When non-null, the influence-function matrix satisfies
#' `vcov(object) == (1 / n^2) * crossprod(influence(object))` up to
#' floating-point noise. Stack the columns from independent fits on the
#' same data with `joint_vcov()` to build a joint asymptotic covariance
#' across estimators, weight schemes, or rater pairs.
#'
#' @param model A `misskappa_estimate` object.
#' @param ... Unused; present for S3 generic conformance.
#'
#' @return A numeric matrix of dimension `n x K` with column names matching
#'   `names(coef(model))`, or `NULL` if the estimator does not expose
#'   influence functions.
#'
#' @export
influence.misskappa_estimate <- function(model, ...) {
  psi <- model$psi
  if (is.null(psi) || prod(dim(psi)) == 0L) return(NULL)
  psi
}

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
#' (`kappa(x, method = "available", weight = "quadratic")`). The two
#' agree under the usual asymptotics; the closed form can be faster on
#' large `n`.
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
#' Counts data assumes **exchangeable iid raters** drawn from a single
#' category distribution. For non-exchangeable raters, use rater-
#' identified data and `kappa(x, method = "fiml")` instead.
#'
#' Conger is not reported (raters are not identified in this input format),
#' and neither IPW nor Gwet are meaningful (per-rater observation rates are
#' aggregated away by the counts representation).
#'
#' @param x A subjects-by-categories non-negative integer matrix.
#' @param method Either `"available"` (moment-based; pooled over observed
#'   pair counts) or `"fiml"` (EM over the composition simplex with multi-
#'   variate hypergeometric weights for completing partial counts). The
#'   two agree exactly when every row of `x` sums to `r_total`; with
#'   partial counts (some `r_i < r_total`) the FIML estimator can be more
#'   efficient.
#' @param weight Weighting scheme: `"identity"` (default; equivalent to
#'   `"unweighted"`), `"linear"`, `"quadratic"`, `"ordinal"`, `"radical"`,
#'   `"ratio"`, `"circular"`, or `"bipolar"`.
#' @param values Optional length-C numeric vector of category scores used
#'   by the metric weightings. Defaults to `1:C`.
#' @param r_total Total number of raters per subject. Defaults to the
#'   maximum observed row sum. Required for `method = "fiml"` when rows
#'   have varying totals.
#' @param em_options Named list of EM options for `method = "fiml"`:
#'   `tol`, `max_iter`, `prune_tol`, `start_alpha`, `info_rcond`.
#'
#' @return A `misskappa_estimate` object with `Fleiss` and
#'   `Brennan-Prediger` coefficients and the 2x2 vcov.
#'
#' @export
kappa_counts <- function(x,
                         method = c("available", "fiml"),
                         weight = c("identity", "unweighted", "linear", "quadratic",
                                    "ordinal", "radical", "ratio", "circular", "bipolar"),
                         values = NULL,
                         r_total = NULL,
                         em_options = list()) {
  method <- match.arg(method)
  weight <- match.arg(weight)
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a matrix or data frame.")
  }
  x_mat <- as.matrix(x)
  if (!is.numeric(x_mat)) stop("'x' must be numeric.")
  storage.mode(x_mat) <- "integer"

  if (method == "available") {
    out <- rcpp_kappa_counts(x = x_mat, weight_type = weight, values = values)
  } else {
    if (is.null(r_total)) r_total <- as.integer(max(rowSums(x_mat)))
    out <- rcpp_kappa_fiml_counts(
      x = x_mat,
      weight_type = weight,
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
      method = method,
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
#' @param x A subjects-by-raters numeric matrix.
#' @param method One of `"available"`, `"ipw"`, `"gwet"`, or `"fiml"`.
#'   `"fiml"` requires `weight = "quadratic"` and dispatches to
#'   [kappa_quadratic_fiml()].
#' @param weight Continuous loss kernel: `"identity"` (binary 0/1),
#'   `"linear"`, `"quadratic"`, `"radical"`, or `"ratio"`. The data range
#'   `[min(x), max(x)]` is used to parameterise the kernel. Only
#'   `"quadratic"` is available with `method = "fiml"`.
#' @param em_options Used only when `method = "fiml"`: named list with `tol`,
#'   `max_iter`, and `fd_h`.
#'
#' @return A `misskappa_estimate` object carrying named coefficients
#'   (`Conger`, `Fleiss`) and the 2x2 asymptotic covariance matrix.
#'
#' @export
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

#' G-wise agreement coefficients for complete rectangular ratings
#'
#' @description
#' Estimates closed-data g-wise Cohen-type and Fleiss-type agreement
#' coefficients using multirater disagreement kernels. This is the
#' Frechet / Hubert family for complete rectangular designs: every subject
#' must have the same set of observed raters, and missing values are not
#' supported.
#'
#' @param x A complete subjects-by-raters matrix or data frame.
#' @param distance Multirater disagreement kernel. `"nominal"` uses Frechet
#'   mode disagreement for categorical ratings; `"absolute"` uses median
#'   absolute deviation; `"quadratic"` uses mean squared deviation; `"hubert"`
#'   uses all-raters-equal disagreement.
#' @param g Arity of the multirater distance. Defaults to all raters
#'   (`ncol(x)`). Must be between 2 and `ncol(x)`.
#' @param max_chance_tuples Maximum number of direct `n^g` item tuples to
#'   evaluate before stopping for continuous distances, or finite category
#'   tuples for categorical distances.
#'
#' @return A `misskappa_estimate` object with `Cohen` and `Fleiss`
#'   coefficients and a 2x2 influence-function covariance matrix.
#'
#' @export
kappa_gwise <- function(x,
                        distance = c("nominal", "absolute", "quadratic", "hubert"),
                        g = NULL,
                        max_chance_tuples = 5000000L) {
  distance <- match.arg(distance)
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

  if (distance %in% c("nominal", "hubert")) {
    if (any(!is.finite(x_mat))) {
      stop("'x' must be complete and finite for g-wise categorical distances.")
    }
    cats <- sort(unique(c(x_mat)))
    x_indexed <- matrix(match(x_mat, cats) - 1L, nrow = nrow(x_mat), ncol = ncol(x_mat))
    storage.mode(x_indexed) <- "integer"
    out <- rcpp_kappa_gwise_categorical(
      x = x_indexed,
      distance_type = distance,
      g = as.integer(g),
      max_chance_tuples = as.integer(max_chance_tuples)
    )
  } else {
    storage.mode(x_mat) <- "double"
    out <- rcpp_kappa_gwise_continuous(
      x = x_mat,
      distance_type = distance,
      g = as.integer(g),
      max_chance_tuples = as.integer(max_chance_tuples)
    )
  }

  estimates <- as.numeric(out$estimates)
  names(estimates) <- c("Cohen", "Fleiss")
  vcov_mat <- out$vcov
  dimnames(vcov_mat) <- list(names(estimates), names(estimates))
  psi_mat <- out$psi
  if (prod(dim(psi_mat)) > 0L) colnames(psi_mat) <- names(estimates)

  structure(
    list(
      estimates = estimates,
      vcov = vcov_mat,
      psi = psi_mat,
      method = "gwise",
      weight = distance,
      distance = distance,
      g = as.integer(g)
    ),
    class = "misskappa_estimate"
  )
}
