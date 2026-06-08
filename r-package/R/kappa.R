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
#' @param x A subjects-by-items numeric matrix or data frame; `NA` and other
#'   non-finite values indicate missing entries. For `estimator = "cat_fiml"`,
#'   entries are integer category codes.
#' @param estimator One of `"pairwise"`, `"cat_fiml"`, or `"nt_fiml"`.
#' @param values Optional numeric vector of scores for observed categories.
#'   When supplied, its length must equal the number of unique finite entries
#'   in `x`.
#' @param em_options Named list tuning the EM fit for the likelihood
#'   estimators. For `"nt_fiml"`: `tol`, `max_iter`, `fd_h`. For `"cat_fiml"`:
#'   `tol`, `max_iter`, `prune_tol`, `start_alpha`, `info_rcond` (the relative
#'   eigenvalue cutoff used when inverting Louis' observed information). Pass
#'   any subset.
#'
#' @return An object of class `misskappa_estimate` carrying one coefficient
#'   named `alpha` and its asymptotic covariance matrix. Methods: `print`,
#'   `coef`, `vcov`, `confint`, `as.data.frame`. The object also carries a
#'   `psi` component (the n-by-1 matrix of per-subject influence functions).
#'
#' @examples
#' # Continuous item battery: the textual subscale of the Holzinger-Swineford
#' # (1939) data. Normal-theory FIML alpha, valid under ignorable missingness.
#' textual <- as.matrix(dat.holzinger1939[, c("x4", "x5", "x6")])
#' fit <- alpha(textual, estimator = "nt_fiml")
#' coef(fit)
#' confint(fit)
#'
#' @examplesIf requireNamespace("psych", quietly = TRUE)
#' # Real item-level missing data: the Neuroticism scale of psych::bfi
#' # (2800 respondents, ~360 with at least one missing item). The FIML
#' # estimator is valid under ignorable missingness.
#' data(bfi, package = "psych")
#' alpha(as.matrix(bfi[, paste0("N", 1:5)]), estimator = "nt_fiml")
#'
#' @export
alpha <- function(x,
                  estimator = c("pairwise", "cat_fiml", "nt_fiml"),
                  values = NULL,
                  em_options = list()) {
  estimator <- match.arg(estimator)

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
  structure(list(estimates = estimates, vcov = vcov_mat, psi = psi_mat,
                 method = method, weight = weight),
            class = "misskappa_estimate")
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
#' `"ipw"` and `"cat_fiml"` take integer category codes and any `weight`.
#' `"pairwise"` and `"nt_fiml"` are quadratic by construction: they treat `x`
#' as numeric scores (pass scores directly, or integer codes with `values`)
#' and require `weight = "quadratic"`.
#'
#' @param x A subjects-by-raters matrix or data frame. Integer category codes
#'   for `"ipw"` / `"cat_fiml"`; numeric scores for `"pairwise"` / `"nt_fiml"`.
#'   `NA` marks missing entries.
#' @param estimator One of `"ipw"`, `"cat_fiml"`, `"pairwise"`, `"nt_fiml"`.
#' @param weight Weighting scheme for `"ipw"` / `"cat_fiml"`: `"nominal"` (the
#'   default), `"linear"`, `"quadratic"`, `"ordinal"`, `"radical"`, `"ratio"`,
#'   `"circular"`, or `"bipolar"`. `"pairwise"` / `"nt_fiml"` require
#'   `"quadratic"`.
#' @param values Optional numeric vector of category scores used by the metric
#'   weightings (and by the quadratic estimators to map integer codes to
#'   scores). Defaults to the sorted unique observed categories.
#' @param em_options Named list of EM options for the likelihood estimators
#'   (`"cat_fiml"`, `"nt_fiml"`). Pass any subset.
#'
#' @return An object of class `misskappa_estimate` carrying the named
#'   coefficient estimates and their asymptotic covariance matrix. Methods:
#'   `print`, `coef`, `vcov`, `confint`, `as.data.frame`. The object also
#'   carries a `psi` component --- the n-by-K matrix of per-subject influence
#'   functions, satisfying `vcov == crossprod(psi) / n^2` --- for power users
#'   who want to build their own contrasts or joint tests.
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
#' @export
kappa <- function(x,
                  estimator = c("ipw", "cat_fiml", "pairwise", "nt_fiml"),
                  weight = c("nominal", "linear", "quadratic", "ordinal",
                             "radical", "ratio", "circular", "bipolar"),
                  values = NULL,
                  em_options = list()) {
  estimator <- match.arg(estimator)
  weight_supplied <- !missing(weight)
  weight <- match.arg(weight)

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

#' @export
print.misskappa_estimate <- function(x, digits = 4, level = 0.95,
                                     transform = c("none", "fisher"), ...) {
  transform <- match.arg(transform)
  cat(sprintf("misskappa: estimator=%s, weight=%s\n", x$method, x$weight))
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
#' identified data and `kappa(x, estimator = "cat_fiml")` instead.
#'
#' Conger is not reported (raters are not identified in this input format),
#' and neither IPW nor Gwet are meaningful (per-rater observation rates are
#' aggregated away by the counts representation).
#'
#' @param x A subjects-by-categories non-negative integer matrix.
#' @param estimator Either `"pairwise"` (moment-based; pooled over observed
#'   pair counts) or `"cat_fiml"` (EM over the composition simplex with
#'   multivariate hypergeometric weights for completing partial counts). The
#'   two agree exactly when every row of `x` sums to `r_total`; with partial
#'   counts (some `r_i < r_total`) the `"cat_fiml"` estimator can be more
#'   efficient.
#' @param weight Weighting scheme: `"nominal"` (default), `"linear"`,
#'   `"quadratic"`, `"ordinal"`, `"radical"`, `"ratio"`, `"circular"`, or
#'   `"bipolar"`.
#' @param values Optional length-C numeric vector of category scores used
#'   by the metric weightings. Defaults to `1:C`.
#' @param r_total Total number of raters per subject. Defaults to the
#'   maximum observed row sum. Required for `estimator = "cat_fiml"` when rows
#'   have varying totals.
#' @param em_options Named list of EM options for `estimator = "cat_fiml"`:
#'   `tol`, `max_iter`, `prune_tol`, `start_alpha`, `info_rcond`.
#'
#' @return A `misskappa_estimate` object with `Fleiss` and
#'   `Brennan-Prediger` coefficients and the 2x2 vcov.
#'
#' @examples
#' # Fleiss (1971) psychiatric diagnoses: 30 subjects, 6 raters, 5 categories,
#' # in counts format (one row per subject, one column per category).
#' kappa_counts(dat.fleiss1971, estimator = "pairwise")
#'
#' # Treat the categories as ordered scores with a quadratic loss.
#' kappa_counts(dat.fleiss1971, estimator = "pairwise", weight = "quadratic")
#'
#' @export
kappa_counts <- function(x,
                         estimator = c("pairwise", "cat_fiml"),
                         weight = c("nominal", "linear", "quadratic",
                                    "ordinal", "radical", "ratio", "circular", "bipolar"),
                         values = NULL,
                         r_total = NULL,
                         em_options = list()) {
  estimator <- match.arg(estimator)
  weight <- match.arg(weight)
  cpp_weight <- if (weight == "nominal") "identity" else weight
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a matrix or data frame.")
  }
  x_mat <- as.matrix(x)
  if (!is.numeric(x_mat)) stop("'x' must be numeric.")
  storage.mode(x_mat) <- "integer"

  if (estimator == "pairwise") {
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
#' @keywords internal
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
