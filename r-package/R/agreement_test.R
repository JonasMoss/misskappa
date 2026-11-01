# Public equality-test front door. It stands on the exposed primitives (coef,
# vcov, and the fit$psi influence-function component). Arbitrary contrasts are
# built from fit$psi directly (see the comparing-coefficients article).

.fit_labels <- function(fits) {
  labs <- names(fits)
  if (is.null(labs)) labs <- rep("", length(fits))
  blank <- !nzchar(labs)
  labs[blank] <- paste0("fit", seq_along(fits))[blank]
  labs
}

# Wald test that one coefficient is equal across >= 2 fits (or = `value` for a
# single fit). Returns an htest. Covariance assembly is the only thing that
# varies: single fit -> its vcov; paired -> influence-function cross-terms;
# independent -> block-diagonal of the per-fit vcovs.
.coef_equality_test <- function(fits, coef, paired, value, labels) {
  if (!all(vapply(fits, inherits, logical(1), "misskappa_estimate")))
    stop("all inputs must be 'misskappa_estimate' objects.", call. = FALSE)
  G <- length(fits)
  if (G > 2L && !identical(value, 0)) {
    stop("'value' must be 0 when comparing more than two fits.", call. = FALSE)
  }
  for (i in seq_len(G))
    if (!coef %in% names(stats::coef(fits[[i]])))
      stop(sprintf('coefficient "%s" not found in %s.', coef, labels[i]),
           call. = FALSE)
  beta <- vapply(fits, function(f) stats::coef(f)[[coef]], numeric(1))

  if (G == 1L) {
    se2  <- stats::vcov(fits[[1L]])[coef, coef]
    stat <- (beta - value)^2 / se2
    df   <- 1L
    method <- sprintf("One-sample Wald test of %s", coef)
    null_value <- stats::setNames(value, coef)
  } else {
    if (paired) {
      V <- .paired_vcov(fits, coef, labels)
    } else {
      V <- diag(vapply(fits, function(f) stats::vcov(f)[coef, coef], numeric(1)), G)
    }
    L     <- cbind(-1, diag(G - 1L))               # "all equal" contrast
    delta <- as.numeric(L %*% beta) - value
    stat  <- as.numeric(crossprod(delta, qr.solve(L %*% V %*% t(L), delta)))
    df    <- G - 1L
    method <- sprintf("%s Wald test of equal %s across %d fits",
                      if (paired) "Paired" else "Independent-sample",
                      coef, G)
    null_value <- if (G == 2L) {
      stats::setNames(value, sprintf("difference in %s", coef))
    } else {
      NULL
    }
  }

  structure(list(
    statistic = stats::setNames(stat, "X-squared"),
    parameter = stats::setNames(df, "df"),
    p.value   = stats::pchisq(stat, df, lower.tail = FALSE),
    estimate  = stats::setNames(beta, labels),
    null.value = null_value,
    alternative = if (is.null(null_value)) NULL else "two.sided",
    method    = method,
    data.name = paste(labels, collapse = ", ")
  ), class = "htest")
}

# Covariance of one coefficient across fits on the same subjects, from the
# per-subject influence rows. Each fit's influence rows are scaled by its own
# subject count, so Cov(theta_j, theta_k) = sum_i psi_ij psi_ik / (n_j n_k).
# Rows are matched by subject id (agreement() names them by row name or row
# number). A subject missing from a fit, because it had no ratings there,
# contributes zero influence to that fit. Fits without ids must have the same
# number of rows and are matched by position.
.paired_vcov <- function(fits, coef, labels) {
  psis <- lapply(seq_along(fits), function(i) {
    m <- fits[[i]]$psi
    if (is.null(m) || prod(dim(m)) == 0L)
      stop(sprintf("paired = TRUE needs influence functions; %s has none.",
                   labels[i]), call. = FALSE)
    m[, coef, drop = FALSE]
  })
  n <- vapply(psis, nrow, numeric(1))
  ids <- lapply(psis, rownames)
  if (any(vapply(ids, is.null, logical(1)))) {
    if (length(unique(n)) != 1L)
      stop("paired fits must have the same number of subjects; got ",
           paste(n, collapse = ", "), ".", call. = FALSE)
    psi <- do.call(cbind, psis)
  } else {
    all_ids <- unique(unlist(ids))
    psi <- vapply(psis, function(m) {
      out <- numeric(length(all_ids))
      out[match(rownames(m), all_ids)] <- m[, 1L]
      out
    }, numeric(length(all_ids)))
    psi <- matrix(psi, nrow = length(all_ids))
  }
  psi <- sweep(psi, 2L, n, "/")
  crossprod(psi)
}

#' Test equality of an agreement coefficient across fits
#'
#' @description
#' Wald test that a coefficient is equal across two or more fits, or equal to
#' `value` for a single fit. With more than two fits it tests that all are
#' equal (`G - 1` degrees of freedom).
#'
#' Fits are independent samples by default, so their variances add. Use
#' `paired = TRUE` for fits on the same subjects (two estimators, two
#' attributes, two rater groups, ...): the dependence then comes from the
#' per-subject influence rows, matched by the row names of the rating data, or
#' by row number when there are none.
#'
#' @param ... Fits from [agreement()] or [agreement_counts()]. Names label the
#'   output.
#' @param coef Coefficient to compare. Defaults to the first one (`"Conger"`
#'   for [agreement()], `"Fleiss"` for [agreement_counts()]).
#' @param paired `FALSE` (default) for independent samples, `TRUE` for fits on
#'   the same subjects.
#' @param value Null value for the coefficient (one fit) or the difference
#'   (two fits). Default `0`.
#'
#' @return An `htest` object.
#'
#' @examples
#' # Two studies with different subjects.
#' agreement_test(gwet  = agreement(dat.gwet2014),
#'                klein = agreement(dat.klein2018))
#'
#' # Two rater pairs on the same 50 subjects.
#' agreement_test(AB = agreement(dat.zapf2016[, 1:2], estimator = "pairwise"),
#'                CD = agreement(dat.zapf2016[, 3:4], estimator = "pairwise"),
#'                paired = TRUE)
#'
#' @export
agreement_test <- function(..., coef = NULL, paired = FALSE, value = 0) {
  fits <- list(...)
  if (length(fits) < 1L) stop("agreement_test() needs at least one fit.",
                              call. = FALSE)
  if (is.null(coef) && inherits(fits[[1L]], "misskappa_estimate")) {
    coef <- names(fits[[1L]]$estimates)[1L]
  }
  .coef_equality_test(fits, coef = coef, paired = paired, value = value,
                      labels = .fit_labels(fits))
}
