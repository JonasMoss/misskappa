#' Joint asymptotic covariance across kappa fits on the same data
#'
#' @description
#' Given two or more `misskappa_estimate` objects fit on the **same
#' subjects in the same order**, returns the joint asymptotic covariance
#' of their stacked coefficient vector. The estimators can differ
#' (available-case, IPW, Gwet), use different weight schemes, or be
#' computed on different rater subsets of the same subjects.
#'
#' The joint matrix is assembled from per-subject influence functions:
#' `V = (1 / n^2) * crossprod(cbind(psi_1, ..., psi_K))`. This is the
#' standard nonparametric IF-based joint asymptotic covariance.
#'
#' All inputs must expose influence functions (currently the categorical
#' raw available-case / IPW / Gwet / FIML estimators, counts-format
#' available-case / FIML estimators, continuous MCAR estimators, and closed
#' rectangular g-wise estimators). Fits from estimators that do not yet expose
#' IFs (quadratic and quadratic-counts) are not supported.
#'
#' Row alignment is the caller's responsibility: the helper only checks
#' that the number of subjects matches across inputs and stacks `psi`
#' column-wise.
#'
#' @param ... Two or more `misskappa_estimate` objects. Optionally named;
#'   names are used as prefixes for the output `dimnames`.
#'
#' @return A square numeric matrix of dimension `K_total x K_total`, where
#'   `K_total` is the total number of coefficients across all inputs.
#'   Row / column names are `<prefix>.<coefficient>` for each fit.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 500; R <- 4
#' x <- matrix(sample.int(5, n * R, replace = TRUE), n, R)
#' x[matrix(runif(n * R) < 0.3, n, R)] <- NA
#' ac  <- kappa(x, method = "available")
#' ipw <- kappa(x, method = "ipw")
#' V   <- joint_vcov(ac = ac, ipw = ipw)
#'
#' # Hausman-style Wald test of kappa_AC = kappa_IPW (Conger).
#' delta <- coef(ac)[["Conger"]] - coef(ipw)[["Conger"]]
#' se <- sqrt(V["ac.Conger", "ac.Conger"]
#'          + V["ipw.Conger", "ipw.Conger"]
#'          - 2 * V["ac.Conger", "ipw.Conger"])
#' z <- delta / se
#' }
#'
#' @export
joint_vcov <- function(...) {
  fits <- list(...)
  if (length(fits) < 2L)
    stop("joint_vcov() requires at least two fits.")
  if (!all(vapply(fits, inherits, logical(1), "misskappa_estimate")))
    stop("All inputs must be 'misskappa_estimate' objects.")

  psis <- lapply(fits, stats::influence)
  null_idx <- vapply(psis, is.null, logical(1))
  if (any(null_idx))
    stop("Inputs ", paste(which(null_idx), collapse = ", "),
         " do not expose influence functions. ",
         "Only fits with non-null influence() currently support joint_vcov().")

  ns <- vapply(psis, nrow, integer(1))
  if (length(unique(ns)) != 1L)
    stop("All fits must have the same number of subjects; got ",
         paste(ns, collapse = ", "), ".")
  n <- ns[[1L]]

  psi_stack <- do.call(cbind, psis)
  V <- crossprod(psi_stack) / (as.numeric(n) * as.numeric(n))

  prefixes <- names(fits)
  if (is.null(prefixes) || any(!nzchar(prefixes)))
    prefixes <- paste0("fit", seq_along(fits))
  block_names <- unname(unlist(Map(function(p, mat) paste(p, colnames(mat),
                                                          sep = "."),
                                   prefixes, psis)))
  dimnames(V) <- list(block_names, block_names)
  V
}

#' Wald tests for one or more misskappa fits
#'
#' @description
#' Tests linear hypotheses of the form `contrast %*% beta = value`, where
#' `beta` is the coefficient vector from one `misskappa_estimate` or the
#' stacked coefficient vector from several fits on the same subjects.
#'
#' With one fit, `vcov(fit)` is used. With several fits, `joint_vcov()` is
#' used, so every fit must expose per-subject influence functions and have
#' rows aligned to the same subjects.
#'
#' @param ... One or more `misskappa_estimate` objects. Multiple fits may be
#'   named; names become prefixes such as `"ac.Conger"` in the stacked
#'   coefficient vector. Unnamed multiple fits are prefixed as `"fit1"`,
#'   `"fit2"`, and so on.
#' @param contrast Numeric vector or matrix defining the linear contrast. For
#'   multiple fits, named vector entries or matrix column names are matched
#'   against the stacked coefficient names.
#' @param value Null value for each contrast row. A scalar is recycled.
#'
#' @return An `htest` object with a chi-square Wald statistic.
#'
#' @export
wald_test <- function(..., contrast, value = 0) {
  fits <- list(...)
  if (length(fits) < 1L)
    stop("wald_test() requires at least one fit.")
  if (!all(vapply(fits, inherits, logical(1), "misskappa_estimate")))
    stop("All inputs must be 'misskappa_estimate' objects.")

  if (length(fits) == 1L) {
    beta <- stats::coef(fits[[1L]])
    V <- stats::vcov(fits[[1L]])
    data_name <- "misskappa_estimate"
  } else {
    V <- joint_vcov(...)
    beta_raw <- unname(unlist(lapply(fits, stats::coef)))
    beta <- stats::setNames(beta_raw, rownames(V))
    data_name <- paste(rownames(V), collapse = ", ")
  }

  L <- normalise_contrast(contrast, names(beta))
  if (!is.numeric(value) || !all(is.finite(value)))
    stop("'value' must be finite numeric.")
  if (length(value) == 1L) value <- rep(value, nrow(L))
  if (length(value) != nrow(L))
    stop("'value' must be scalar or have one entry per contrast row.")

  estimate <- as.numeric(L %*% beta)
  delta <- estimate - as.numeric(value)
  LVL <- L %*% V %*% t(L)
  stat <- as.numeric(t(delta) %*% qr.solve(LVL, delta))
  df <- nrow(L)
  names(estimate) <- rownames(L)
  names(value) <- rownames(L)

  structure(
    list(
      statistic = stats::setNames(stat, "X-squared"),
      parameter = stats::setNames(df, "df"),
      p.value = stats::pchisq(stat, df = df, lower.tail = FALSE),
      estimate = estimate,
      null.value = value,
      alternative = "two-sided",
      method = "Wald test for misskappa coefficients",
      data.name = data_name
    ),
    class = "htest"
  )
}

normalise_contrast <- function(contrast, beta_names) {
  if (is.character(contrast)) {
    missing <- setdiff(contrast, beta_names)
    if (length(missing) > 0L)
      stop("Unknown coefficient(s): ", paste(missing, collapse = ", "))
    L <- diag(length(beta_names))[match(contrast, beta_names), , drop = FALSE]
    rownames(L) <- contrast
    colnames(L) <- beta_names
    return(L)
  }

  if (!is.numeric(contrast))
    stop("'contrast' must be numeric or character.")

  if (is.null(dim(contrast))) {
    if (!is.null(names(contrast))) {
      missing <- setdiff(names(contrast), beta_names)
      if (length(missing) > 0L)
        stop("Unknown coefficient(s): ", paste(missing, collapse = ", "))
      L <- matrix(0, nrow = 1L, ncol = length(beta_names),
                  dimnames = list("contrast", beta_names))
      L[1L, names(contrast)] <- contrast
      return(L)
    }
    if (length(contrast) != length(beta_names))
      stop("Unnamed contrast vector must have length ", length(beta_names), ".")
    L <- matrix(contrast, nrow = 1L,
                dimnames = list("contrast", beta_names))
    return(L)
  }

  L <- as.matrix(contrast)
  if (!is.numeric(L))
    stop("'contrast' must be numeric.")
  if (!is.null(colnames(L))) {
    missing <- setdiff(colnames(L), beta_names)
    if (length(missing) > 0L)
      stop("Unknown coefficient(s): ", paste(missing, collapse = ", "))
    L_full <- matrix(0, nrow = nrow(L), ncol = length(beta_names),
                     dimnames = list(rownames(L), beta_names))
    L_full[, colnames(L)] <- L
    L <- L_full
  } else {
    if (ncol(L) != length(beta_names))
      stop("Unnamed contrast matrix must have ", length(beta_names), " columns.")
    colnames(L) <- beta_names
  }
  if (is.null(rownames(L))) rownames(L) <- paste0("contrast", seq_len(nrow(L)))
  L
}
