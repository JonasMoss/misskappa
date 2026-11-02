# The misskappa_estimate object: constructor, shared argument checks, and the
# S3 methods (print, coef, vcov, confint, as.data.frame).

# Build a misskappa_estimate from a backend result `out` carrying `estimates`,
# `vcov`, `psi` (n x K influence rows), and optionally the FIML diagnostics
# `null_frac` and `se_cutoff_ratio`. Extra named
# fields (moments, bp_chance, g, W, ...) are appended; NULL extras are skipped.
.new_estimate <- function(out, coefficients, method, weight, ...) {
  K <- length(coefficients)
  estimates <- stats::setNames(as.numeric(out$estimates), coefficients)
  vcov_mat <- as.matrix(out$vcov)
  dimnames(vcov_mat) <- list(coefficients, coefficients)
  psi <- out$psi
  if (!is.null(psi) && prod(dim(psi)) > 0L) {
    psi <- as.matrix(psi)
    colnames(psi) <- coefficients
  }
  fit <- list(
    estimates = estimates,
    vcov = vcov_mat,
    psi = psi,
    n = if (!is.null(psi) && nrow(psi) > 0L) nrow(psi) else NA_integer_,
    method = method,
    weight = weight
  )
  if (length(out$null_frac) == K) {
    fit$null_frac <- stats::setNames(as.numeric(out$null_frac), coefficients)
  }
  if (length(out$se_cutoff_ratio) == K) {
    fit$se_cutoff_ratio <- stats::setNames(as.numeric(out$se_cutoff_ratio),
                                           coefficients)
  }
  extra <- list(...)
  extra <- extra[!vapply(extra, is.null, logical(1))]
  fit[names(extra)] <- extra
  undefined <- !is.finite(estimates)
  if (any(undefined)) {
    # A zero chance disagreement (every rating in one category, or constant
    # scores) leaves the coefficient undefined; its zero "variance" is not a
    # standard error, so blank it rather than report false precision.
    fit$vcov[undefined, ] <- NA_real_
    fit$vcov[, undefined] <- NA_real_
    warning(sprintf(paste0(
      "%s undefined: the chance disagreement is zero, as when every rating ",
      "falls in one category."),
      paste(coefficients[undefined], collapse = ", ")), call. = FALSE)
  }
  structure(fit, class = "misskappa_estimate")
}

# Remove coefficients from a fit, keeping estimates, vcov, influence rows, the
# FIML diagnostics, and the delete-one summaries consistent.
.drop_coefficients <- function(fit, drop) {
  keep <- !names(fit$estimates) %in% drop
  if (all(keep)) return(fit)
  fit$estimates <- fit$estimates[keep]
  fit$vcov <- fit$vcov[keep, keep, drop = FALSE]
  if (!is.null(fit$psi) && prod(dim(fit$psi)) > 0L) {
    fit$psi <- fit$psi[, keep, drop = FALSE]
  }
  if (!is.null(fit$null_frac)) fit$null_frac <- fit$null_frac[keep]
  if (!is.null(fit$se_cutoff_ratio)) {
    fit$se_cutoff_ratio <- fit$se_cutoff_ratio[keep]
  }
  d1 <- fit$delete_one
  if (!is.null(d1)) {
    d1$estimate <- d1$estimate[keep]
    d1$leave_one_out <- d1$leave_one_out[, keep, drop = FALSE]
    if (length(d1$chance) == length(keep)) {
      d1$chance <- d1$chance[keep]
      d1$chance_leave_one_out <- d1$chance_leave_one_out[, keep, drop = FALSE]
    }
    fit$delete_one <- d1
  }
  if ("Brennan-Prediger" %in% drop) fit$bp_chance <- NULL
  fit
}

# Error on `...` arguments that do not apply in the current mode, instead of
# silently ignoring a typo or an argument meant for another input shape.
.check_dots <- function(dots, allowed, context) {
  if (!length(dots)) return(invisible(NULL))
  nm <- names(dots)
  if (is.null(nm) || any(!nzchar(nm))) {
    stop("arguments passed through '...' must be named.", call. = FALSE)
  }
  bad <- setdiff(nm, allowed)
  if (length(bad)) {
    stop(sprintf(
      "argument%s not used by %s: %s.%s",
      if (length(bad) > 1L) "s" else "", context,
      paste(bad, collapse = ", "),
      if (length(allowed)) {
        sprintf(" Allowed here: %s.", paste(allowed, collapse = ", "))
      } else {
        ""
      }), call. = FALSE)
  }
  invisible(NULL)
}

# EM tuning options read by each likelihood estimator.
.em_options_categorical <- c("tol", "max_iter", "prune_tol", "start_alpha",
                             "info_rcond")
.em_options_normal <- c("tol", "max_iter")

.check_em_options <- function(em_options, allowed, estimator) {
  if (!is.list(em_options)) {
    stop("'em_options' must be a named list.", call. = FALSE)
  }
  if (!length(em_options)) return(em_options)
  if (!length(allowed)) {
    stop(sprintf('em_options does not apply to estimator = "%s".', estimator),
         call. = FALSE)
  }
  nm <- names(em_options)
  if (is.null(nm) || any(!nzchar(nm))) {
    stop("'em_options' must be a named list.", call. = FALSE)
  }
  bad <- setdiff(nm, allowed)
  if (length(bad)) {
    stop(sprintf('unknown em_options for estimator = "%s": %s. Allowed: %s.',
                 estimator, paste(bad, collapse = ", "),
                 paste(allowed, collapse = ", ")), call. = FALSE)
  }
  em_options
}

#' @export
print.misskappa_estimate <- function(x, digits = 4, level = 0.95,
                                     method = c("wald", "jackknife"),
                                     transform = c("none", "fisher",
                                                   "log_complement"), ...) {
  method <- match.arg(method)
  transform <- match.arg(transform)
  hdr <- sprintf("misskappa: estimator = %s", x$method)
  if (!is.null(x$weight)) hdr <- sprintf("%s, weight = %s", hdr, x$weight)
  if (!is.null(x$g) && x$g > 2L) hdr <- sprintf("%s, g = %d", hdr, x$g)
  if (!is.null(x$n) && !is.na(x$n)) hdr <- sprintf("%s, n = %d", hdr, x$n)
  cat(hdr, "\n", sep = "")
  se <- sqrt(diag(x$vcov))
  ci <- stats::confint(x, level = level, method = method,
                       transform = transform)
  tab <- data.frame(estimate = x$estimates, se = se,
                    lower = ci[, 1L], upper = ci[, 2L],
                    row.names = names(x$estimates))
  print(round(tab, digits = digits), ...)
  if (!is.null(x$fallback)) {
    cat(sprintf('Note: estimator = "%s" failed (%s); this is the "%s" fallback.\n',
                x$fallback$from, sub("[.]$", "", x$fallback$reason), x$method))
  }
  if (isFALSE(x$moments$converged)) {
    cat("Note: the EM did not converge; the estimates may be unreliable.\n")
  }
  invisible(x)
}

#' @export
coef.misskappa_estimate <- function(object, ...) object$estimates

#' @export
vcov.misskappa_estimate <- function(object, ...) object$vcov

#' Confidence intervals for misskappa coefficients
#'
#' @description
#' Wald or exact-delete-one jackknife confidence intervals for a
#' `misskappa_estimate`. Wald intervals use the estimator covariance and a
#' normal critical value. Jackknife intervals transform every complete-data
#' leave-one-subject-out estimate directly, use the ordinary delete-one
#' jackknife standard error, and use a Student t critical value with `n - 1`
#' degrees of freedom.
#'
#' Besides the natural and Fisher scales, `transform = "log_complement"` uses
#' the upper-bound-only transform `-log(1 - theta)`, which stays defined for
#' coefficients below -1. At an exact boundary (perfect agreement), a
#' degenerate transformed jackknife is replaced by a one-sided hybrid: a
#' Clopper--Pearson upper limit for the probability of any item disagreement is
#' divided by a delete-one t lower limit for chance disagreement. When that
#' fallback is unavailable the interval is returned as `NA` with a warning.
#'
#' @section Jackknife availability:
#' The exact delete-one summaries behind `method = "jackknife"` are attached
#' automatically to fits with at least three subjects in these cases:
#' \itemize{
#'   \item [agreement()] on complete scalar ratings: `"ipw"` and `"cat_fiml"`
#'     (`g = 2`), and `"pairwise"`.
#'   \item [agreement_counts()]: `"pairwise"` for any row totals, and
#'     `"cat_fiml"` when every row has the same total.
#' }
#' Fits on incomplete data, `g > 2`, vector-valued ratings, the
#' compound-symmetry estimators, and categorical FIML with `flatten > 0` support
#' only `method = "wald"`.
#'
#' @param object A `misskappa_estimate` object.
#' @param parm Optional subset of coefficients (names or indices), defaults to
#'   all.
#' @param level Confidence level.
#' @param method Either `"wald"` or `"jackknife"`. See the section on
#'   jackknife availability.
#' @param transform One of `"none"`, `"fisher"`, or `"log_complement"`.
#'   Wald uses a delta-method standard error on the selected scale; jackknife
#'   transforms the exact delete-one estimates directly.
#' @param ... Unused, present for S3 generic conformance.
#'
#' @return A two-column numeric matrix of lower and upper limits, one row per
#'   coefficient.
#'
#' @references
#' Efron, B., and Tibshirani, R. J. (1993). *An Introduction to the
#' Bootstrap*. Chapman and Hall, Chapter 11.
#'
#' Clopper, C. J., and Pearson, E. S. (1934). The use of confidence or
#' fiducial limits illustrated in the case of the binomial. *Biometrika*,
#' 26(4), 404--413. \doi{10.1093/biomet/26.4.404}.
#'
#' @examples
#' fit <- agreement(dat.gwet2014, estimator = "ipw")
#' confint(fit)                       # natural-scale Wald interval
#' confint(fit, transform = "fisher") # Fisher z interval, always within (-1, 1)
#' complete_fit <- agreement_counts(dat.fleiss1971)
#' confint(complete_fit, method = "jackknife", transform = "fisher")
#'
#' @export
confint.misskappa_estimate <- function(object, parm = NULL, level = 0.95,
                                       method = c("wald", "jackknife"),
                                       transform = c("none", "fisher",
                                                     "log_complement"), ...) {
  method <- match.arg(method)
  transform <- match.arg(transform)
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("'level' must be a single number strictly between zero and one.",
         call. = FALSE)
  }
  est <- object$estimates

  forward <- switch(transform,
    none = function(x) x,
    fisher = function(x) {
      out <- x
      out[] <- NA_real_
      interior <- is.finite(x) & abs(x) < 1
      out[interior] <- atanh(x[interior])
      out
    },
    log_complement = function(x) {
      out <- x
      out[] <- NA_real_
      interior <- is.finite(x) & x < 1
      out[interior] <- -log1p(-x[interior])
      out
    }
  )
  inverse <- switch(transform,
    none = function(x) x,
    fisher = tanh,
    log_complement = function(x) 1 - exp(-x)
  )

  fallback <- rep(FALSE, length(est))
  if (method == "wald") {
    if (transform == "fisher" && any(!is.finite(est) | abs(est) >= 1)) {
      warning("Fisher transform requires |estimate| < 1; returning NA limits ",
              "for boundary coefficients.", call. = FALSE)
    }
    if (transform == "log_complement" && any(!is.finite(est) | est >= 1)) {
      warning("Log-complement transform requires estimate < 1; returning NA ",
              "limits for boundary coefficients.", call. = FALSE)
    }
    se <- sqrt(diag(object$vcov))
    critical <- stats::qnorm((1 + level) / 2)
    center <- forward(est)
    derivative <- switch(transform,
      none = rep(1, length(est)),
      fisher = ifelse(abs(est) < 1, 1 / (1 - est^2), NA_real_),
      log_complement = ifelse(est < 1, 1 / (1 - est), NA_real_)
    )
    transformed_se <- se * derivative
    lo <- inverse(center - critical * transformed_se)
    hi <- inverse(center + critical * transformed_se)
  } else {
    delete_one <- object$delete_one
    if (is.null(delete_one) || is.null(delete_one$leave_one_out)) {
      stop("Exact delete-one summaries are unavailable for this fit. ",
           "Jackknife intervals need complete data and a supported estimator; ",
           "see ?confint.misskappa_estimate.", call. = FALSE)
    }
    loo <- as.matrix(delete_one$leave_one_out)
    n <- nrow(loo)
    if (n < 3L || ncol(loo) != length(est)) {
      stop("Malformed exact delete-one summary in the fitted object.",
           call. = FALSE)
    }
    transformed <- forward(loo)
    center <- forward(est)
    transformed_mean <- colMeans(transformed, na.rm = FALSE)
    transformed_se <- sqrt((n - 1) / n *
      colSums(sweep(transformed, 2L, transformed_mean, "-")^2))
    critical <- stats::qt((1 + level) / 2, df = n - 1)
    lo <- inverse(center - critical * transformed_se)
    hi <- inverse(center + critical * transformed_se)

    invalid <- !is.finite(lo) | !is.finite(hi) |
      !is.finite(transformed_se) | transformed_se <= sqrt(.Machine$double.eps)
    has_boundary <-
      length(delete_one$item_disagreement) == n &&
      length(delete_one$chance) == length(est) &&
      identical(dim(delete_one$chance_leave_one_out), c(n, length(est)))
    can_fallback <- invalid & has_boundary
    if (any(can_fallback)) {
      events <- sum(delete_one$item_disagreement > sqrt(.Machine$double.eps))
      alpha_tail <- (1 - level) / 2
      probability_upper <- if (events == n) 1 else
        stats::qbeta(1 - alpha_tail, events + 1, n - events)
      chance_loo <- as.matrix(delete_one$chance_leave_one_out)
      chance_mean <- colMeans(chance_loo)
      chance_se <- sqrt((n - 1) / n *
        colSums(sweep(chance_loo, 2L, chance_mean, "-")^2))
      chance_lower <- delete_one$chance - critical * chance_se
      usable <- can_fallback & is.finite(chance_lower) & chance_lower > 0
      lo[usable] <- 1 - probability_upper / chance_lower[usable]
      hi[usable] <- 1
      fallback[usable] <- TRUE
      invalid[usable] <- FALSE
    }
    if (any(invalid)) {
      lo[invalid] <- NA_real_
      hi[invalid] <- NA_real_
      warning("The transformed jackknife is undefined or degenerate and no ",
              "kappa boundary fallback is available; returning NA limits.",
              call. = FALSE)
    }
  }

  ci <- cbind(lo, hi)
  colnames(ci) <- c(sprintf("%.1f %%", 100 * (1 - level) / 2),
                    sprintf("%.1f %%", 100 * (1 + level) / 2))
  rownames(ci) <- names(est)
  if (!is.null(parm)) {
    ci <- ci[parm, , drop = FALSE]
    fallback <- fallback[parm]
  }
  attr(ci, "method") <- method
  attr(ci, "transform") <- transform
  attr(ci, "boundary_fallback") <- stats::setNames(fallback, rownames(ci))
  ci
}

#' Coerce a misskappa fit to a data frame
#'
#' @description
#' One row per coefficient with the estimator, the estimate, its standard
#' error, and confidence limits from [confint()][confint.misskappa_estimate].
#' Fits from several calls can be stacked with `rbind()`.
#'
#' @param x A `misskappa_estimate` object.
#' @param row.names,optional Passed for S3 conformance. `row.names` sets the
#'   row names of the result.
#' @param level Confidence level for the `lower` and `upper` columns.
#' @param ... Further arguments to [confint()][confint.misskappa_estimate],
#'   such as `method` or `transform`.
#'
#' @return A data frame with columns `estimator`, `coefficient`, `estimate`,
#'   `se`, `lower`, and `upper`.
#'
#' @examples
#' fit <- agreement(dat.gwet2014, estimator = "ipw")
#' as.data.frame(fit)
#'
#' @export
as.data.frame.misskappa_estimate <- function(x, row.names = NULL,
                                             optional = FALSE, level = 0.95,
                                             ...) {
  ci <- stats::confint(x, level = level, ...)
  data.frame(
    estimator = rep(x$method, length(x$estimates)),
    coefficient = names(x$estimates),
    estimate = unname(x$estimates),
    se = unname(sqrt(diag(x$vcov))),
    lower = unname(ci[, 1L]),
    upper = unname(ci[, 2L]),
    row.names = row.names,
    stringsAsFactors = FALSE
  )
}
