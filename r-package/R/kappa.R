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
#'   `tol`, `max_iter`, `prune_tol`, `start_alpha`. Pass any subset.
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

  structure(
    list(
      estimates = estimates,
      vcov = vcov_mat,
      method = method,
      weight = weight
    ),
    class = "misskappa_estimate"
  )
}

#' @export
print.misskappa_estimate <- function(x, digits = 4, ...) {
  cat(sprintf("misskappa: method=%s, weight=%s\n", x$method, x$weight))
  se <- sqrt(diag(x$vcov))
  tab <- data.frame(estimate = x$estimates, se = se, row.names = names(x$estimates))
  print(round(tab, digits = digits), ...)
  invisible(x)
}

#' @export
coef.misskappa_estimate <- function(object, ...) object$estimates

#' @export
vcov.misskappa_estimate <- function(object, ...) object$vcov

#' @export
confint.misskappa_estimate <- function(object, parm = NULL, level = 0.95, ...) {
  est <- object$estimates
  se <- sqrt(diag(object$vcov))
  z <- stats::qnorm((1 + level) / 2)
  lo <- est - z * se
  hi <- est + z * se
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
#'
#' @return A `misskappa_estimate` object with `Conger`, `Fleiss`,
#'   `Brennan-Prediger` coefficients and the 3x3 vcov.
#'
#' @export
kappa_quadratic <- function(x, values) {
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a matrix or data frame.")
  }
  x_mat <- as.matrix(x)
  storage.mode(x_mat) <- "double"
  if (!is.numeric(values)) stop("'values' must be numeric.")

  out <- rcpp_kappa_quadratic(x = x_mat, values = as.numeric(values))
  estimates <- as.numeric(out$estimates)
  names(estimates) <- c("Conger", "Fleiss", "Brennan-Prediger")
  vcov_mat <- out$vcov
  dimnames(vcov_mat) <- list(names(estimates), names(estimates))
  structure(
    list(
      estimates = estimates,
      vcov = vcov_mat,
      method = "quadratic",
      weight = "quadratic"
    ),
    class = "misskappa_estimate"
  )
}

#' Closed-form quadratic-loss kappa for counts-format input
#'
#' @description
#' Counts-format counterpart of `kappa_quadratic()`. Useful when only the
#' per-subject category counts are available (and not the rater-level data).
#'
#' @param x A subjects-by-categories non-negative integer matrix.
#' @param values Length-C numeric vector of category scores.
#' @param r_total Total number of raters per subject.
#'
#' @return A `misskappa_estimate` object with `Fleiss` and
#'   `Brennan-Prediger` coefficients and the 2x2 vcov.
#'
#' @export
kappa_quadratic_counts <- function(x, values, r_total) {
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a matrix or data frame.")
  }
  x_mat <- as.matrix(x)
  storage.mode(x_mat) <- "integer"
  if (!is.numeric(values)) stop("'values' must be numeric.")
  if (!is.numeric(r_total) || length(r_total) != 1L || r_total < 2) {
    stop("'r_total' must be an integer >= 2.")
  }

  out <- rcpp_kappa_quadratic_counts(
    x = x_mat, values = as.numeric(values), r_total = as.integer(r_total)
  )
  estimates <- as.numeric(out$estimates)
  names(estimates) <- c("Fleiss", "Brennan-Prediger")
  vcov_mat <- out$vcov
  dimnames(vcov_mat) <- list(names(estimates), names(estimates))
  structure(
    list(
      estimates = estimates,
      vcov = vcov_mat,
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
#'   `tol`, `max_iter`, `prune_tol`, `start_alpha`.
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

  structure(
    list(
      estimates = estimates,
      vcov = vcov_mat,
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
#' real-valued ratings under MCAR missingness (`"available"`, `"ipw"`,
#' `"gwet"`). Missing entries are encoded as `NA` (or any non-finite
#' value).
#'
#' Brennan-Prediger is not reported for continuous data; the chance-
#' disagreement baseline requires a finite number of categories.
#'
#' @param x A subjects-by-raters numeric matrix.
#' @param method One of `"available"`, `"ipw"`, `"gwet"`.
#' @param weight Continuous loss kernel: `"identity"` (binary 0/1),
#'   `"linear"`, `"quadratic"`, `"radical"`, or `"ratio"`. The data range
#'   `[min(x), max(x)]` is used to parameterise the kernel.
#'
#' @return A `misskappa_estimate` object carrying named coefficients
#'   (`Conger`, `Fleiss`) and the 2x2 asymptotic covariance matrix.
#'
#' @export
kappa_continuous <- function(x,
                             method = c("available", "ipw", "gwet"),
                             weight = c("quadratic", "linear", "identity",
                                        "unweighted", "radical", "ratio")) {
  method <- match.arg(method)
  weight <- match.arg(weight)

  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a matrix or data frame.")
  }
  x_mat <- as.matrix(x)
  if (!is.numeric(x_mat)) stop("'x' must be numeric.")
  storage.mode(x_mat) <- "double"

  out <- rcpp_kappa_continuous(
    x = x_mat,
    method = method,
    weight_type = weight
  )

  estimates <- as.numeric(out$estimates)
  names(estimates) <- c("Conger", "Fleiss")
  vcov_mat <- out$vcov
  dimnames(vcov_mat) <- list(names(estimates), names(estimates))

  structure(
    list(
      estimates = estimates,
      vcov = vcov_mat,
      method = method,
      weight = weight
    ),
    class = "misskappa_estimate"
  )
}
