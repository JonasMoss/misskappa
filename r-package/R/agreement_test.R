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
  for (i in seq_len(G))
    if (!coef %in% names(stats::coef(fits[[i]])))
      stop(sprintf('coefficient "%s" not found in %s.', coef, labels[i]))
  beta <- vapply(fits, function(f) stats::coef(f)[[coef]], numeric(1))

  if (G == 1L) {
    se2  <- stats::vcov(fits[[1L]])[coef, coef]
    stat <- (beta - value)^2 / se2
    df   <- 1L
    method <- sprintf("One-sample Wald test that %s = %g", coef, value)
  } else {
    if (paired) {
      psis <- lapply(seq_len(G), function(i) {
        m <- fits[[i]]$psi
        if (is.null(m) || prod(dim(m)) == 0L)
          stop(sprintf("paired = TRUE needs influence functions; %s has none.",
                       labels[i]))
        m[, coef]
      })
      if (length(unique(lengths(psis))) != 1L)
        stop("paired fits must have the same number of subjects; got ",
             paste(lengths(psis), collapse = ", "), ".")
      ids <- lapply(fits, function(f) rownames(f$psi))   # alignment check, if ids present
      if (!any(vapply(ids, is.null, logical(1))) && length(unique(ids)) != 1L)
        stop("paired fits are not row-aligned (their fit$psi subject ids differ).")
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
    method <- sprintf("%s test of equal %s across %d fits",
                      if (paired) "Paired (dependent)" else "Independent-sample",
                      coef, G)
  }

  structure(list(
    statistic = stats::setNames(stat, "X-squared"),
    parameter = stats::setNames(df, "df"),
    p.value   = stats::pchisq(stat, df, lower.tail = FALSE),
    estimate  = stats::setNames(beta, labels),
    method    = method,
    data.name = paste(labels, collapse = ", ")
  ), class = "htest")
}

#' Test equality of a kappa agreement coefficient across fits
#'
#' @description
#' Wald test that a chosen kappa coefficient is equal across two or more
#' [kappa()] fits (or, with a single fit, equal to `value`). Pass the fits to
#' `...`; name them to label the output.
#'
#' Set `paired = TRUE` when the fits are computed on the **same subjects in the
#' same row order** (two estimators, two weight schemes, two rated attributes,
#' rater pairs, ...); the dependence is taken from the per-subject influence
#' functions (`fit$psi`). Set `paired = FALSE` for independent samples, where
#' the variances add. With more than two fits this is the joint hypothesis that
#' all are equal (`G - 1` degrees of freedom).
#'
#' @param ... Two or more `misskappa_estimate` objects from [kappa()] (a single
#'   fit tests `coef = value`). Names become labels in the output.
#' @param coef Coefficient to compare; default `"Conger"`.
#' @param paired `TRUE` (default) for same-subject fits, `FALSE` for independent
#'   samples.
#' @param value Null value: the coefficient (single fit) or the pairwise
#'   difference (default `0`).
#'
#' @return An `htest` object.
#'
#' @examples
#' # Paired: are two rater pairs equally consistent on the same 50 subjects
#' # (Zapf 2016)? The dependence is read from the per-subject influence
#' # functions, so paired = TRUE (the default) is appropriate.
#' kappa_test(AB = kappa(dat.zapf2016[, 1:2], estimator = "pairwise"),
#'            CD = kappa(dat.zapf2016[, 3:4], estimator = "pairwise"),
#'            coef = "Conger")
#'
#' # Independent samples (different subjects in two studies): variances add.
#' kappa_test(gwet  = kappa(dat.gwet2014,  estimator = "ipw"),
#'            klein = kappa(dat.klein2018, estimator = "ipw"),
#'            coef = "Conger", paired = FALSE)
#'
#' @export
kappa_test <- function(..., coef = "Conger", paired = TRUE, value = 0) {
  fits <- list(...)
  if (length(fits) < 1L) stop("kappa_test() needs at least one fit.")
  .coef_equality_test(fits, coef = coef, paired = paired, value = value,
                      labels = .fit_labels(fits))
}

#' Test equality of coefficient alpha across fits
#'
#' @description
#' Wald test that coefficient alpha is equal across two or more [alpha()] fits
#' (or, with a single fit, equal to `value`). See [kappa_test()] for the
#' `paired` semantics; `alpha` is the only coefficient, so there is no `coef`
#' argument.
#'
#' @param ... Two or more `misskappa_estimate` objects from [alpha()]. Names
#'   become labels in the output.
#' @param paired `TRUE` (default) for same-subject fits, `FALSE` for independent
#'   samples.
#' @param value Null value for the difference (or, with one fit, for alpha).
#'   Default `0`.
#'
#' @return An `htest` object.
#'
#' @examples
#' # Do the three Holzinger-Swineford (1939) subscales have equal reliability?
#' # The same students take all three, so the fits are paired (G-way, df = 2).
#' subs <- list(visual  = c("x1", "x2", "x3"),
#'              textual = c("x4", "x5", "x6"),
#'              speed   = c("x7", "x8", "x9"))
#' fits <- lapply(subs, function(v)
#'   alpha(dat.holzinger1939[, v], estimator = "nt_fiml"))
#' do.call(alpha_test, c(fits, list(paired = TRUE)))
#'
#' @examplesIf requireNamespace("psych", quietly = TRUE)
#' # Real missing data with groups: psych::bfi. The same 2800 respondents take
#' # all five Big Five scales; reverse-key the negatively-worded items first.
#' data(bfi, package = "psych")
#' neg <- c("A1", "C4", "C5", "E1", "E2", "O2", "O5")
#' bfi[neg] <- 7 - bfi[neg]
#'
#' # Dependent: are Neuroticism and Extraversion equally reliable? The same
#' # respondents answer both, so the estimates are paired.
#' N <- paste0("N", 1:5)
#' E <- paste0("E", 1:5)
#' alpha_test(N = alpha(bfi[, N], estimator = "nt_fiml"),
#'            E = alpha(bfi[, E], estimator = "nt_fiml"),
#'            paired = TRUE)
#'
#' # Independent: is Conscientiousness equally reliable across the two genders?
#' g <- split(seq_len(nrow(bfi)), bfi$gender)
#' C <- paste0("C", 1:5)
#' alpha_test(
#'   men   = alpha(bfi[g[["1"]], C], estimator = "nt_fiml"),
#'   women = alpha(bfi[g[["2"]], C], estimator = "nt_fiml"),
#'   paired = FALSE)
#'
#' @export
alpha_test <- function(..., paired = TRUE, value = 0) {
  fits <- list(...)
  if (length(fits) < 1L) stop("alpha_test() needs at least one fit.")
  .coef_equality_test(fits, coef = "alpha", paired = paired, value = value,
                      labels = .fit_labels(fits))
}
