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
