# Public equality-test front doors. These stand on the exposed primitives
# (coef, vcov, and the fit$psi influence-function component); the general
# linear-hypothesis engine (joint_vcov / wald_test) stays internal for the
# rare arbitrary-contrast case.

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
    stop("all inputs must be 'misskappa_estimate' objects.")
  G <- length(fits)
  if (G > 2L && !identical(value, 0)) {
    stop("'value' must be 0 when comparing more than two fits.", call. = FALSE)
  }
  for (i in seq_len(G))
    if (!coef %in% names(stats::coef(fits[[i]])))
      stop(sprintf('coefficient "%s" not found in %s.', coef, labels[i]))
  beta <- vapply(fits, function(f) stats::coef(f)[[coef]], numeric(1))

  if (G == 1L) {
    se2  <- stats::vcov(fits[[1L]])[coef, coef]
    stat <- (beta - value)^2 / se2
    df   <- 1L
    method <- sprintf("One-sample Wald test of %s", coef)
    null_value <- stats::setNames(value, coef)
  } else {
    if (paired) {
      psis <- lapply(seq_len(G), function(i) {
        m <- fits[[i]]$psi
        if (is.null(m) || prod(dim(m)) == 0L)
          stop(sprintf("paired = TRUE needs influence functions; %s has none.",
                       labels[i]), call. = FALSE)
        m[, coef]
      })
      if (length(unique(lengths(psis))) != 1L)
        stop("paired fits must have the same number of subjects; got ",
             paste(lengths(psis), collapse = ", "), ".")
      ids <- lapply(fits, function(f) rownames(f$psi))  # checked when present
      if (!any(vapply(ids, is.null, logical(1))) && length(unique(ids)) != 1L)
        stop("paired fits are not row-aligned (their fit$psi subject ids differ).",
             call. = FALSE)
      psi <- do.call(cbind, psis)
      n   <- nrow(psi)
      V   <- crossprod(psi) / (as.numeric(n) * as.numeric(n))
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

#' Test equality of a kappa agreement coefficient across fits
#'
#' @description
#' Wald test that a chosen kappa coefficient is equal across two or more
#' [agreement()] fits (or, with a single fit, equal to `value`). Pass the fits to
#' `...`, naming them to label the output.
#'
#' By default the fits are treated as independent samples, so their variances
#' add. Set `paired = TRUE` when the fits are computed on the **same subjects in
#' the same row order** (two estimators, two weight schemes, two rated
#' attributes, rater pairs, ...). The dependence is then taken from the
#' per-subject influence functions (`fit$psi`). With more than two fits this is
#' the joint hypothesis that all are equal (`G - 1` degrees of freedom).
#'
#' @param ... One or more `misskappa_estimate` objects from [agreement()] or
#'   [agreement_counts()] (a single fit tests `coef = value`). Names become labels
#'   in the output.
#' @param coef Coefficient to compare. Defaults to the first coefficient of the
#'   first fit (`"Conger"` for [agreement()], `"Fleiss"` for [agreement_counts()]).
#' @param paired `FALSE` (default) for independent samples, `TRUE` for fits on
#'   the same subjects in the same row order.
#' @param value Null value: the coefficient (single fit) or the difference
#'   between two fits (default `0`).
#'
#' @return An `htest` object.
#'
#' @examples
#' # Paired: are two rater pairs equally consistent on the same 50 subjects
#' # (Zapf 2016)? The dependence is read from the per-subject influence
#' # functions.
#' agreement_test(AB = agreement(dat.zapf2016[, 1:2], estimator = "pairwise"),
#'            CD = agreement(dat.zapf2016[, 3:4], estimator = "pairwise"),
#'            paired = TRUE)
#'
#' # Independent samples (different subjects in two studies): variances add.
#' agreement_test(gwet  = agreement(dat.gwet2014,  estimator = "ipw"),
#'            klein = agreement(dat.klein2018, estimator = "ipw"))
#'
#' @export
agreement_test <- function(..., coef = NULL, paired = FALSE, value = 0) {
  fits <- list(...)
  if (length(fits) < 1L) stop("agreement_test() needs at least one fit.")
  if (is.null(coef) && inherits(fits[[1L]], "misskappa_estimate")) {
    coef <- names(fits[[1L]]$estimates)[1L]
  }
  .coef_equality_test(fits, coef = coef, paired = paired, value = value,
                      labels = .fit_labels(fits))
}
