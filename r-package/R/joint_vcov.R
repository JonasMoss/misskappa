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
#' raw estimators: `kappa(method = "available")`, `kappa(method = "ipw")`,
#' `kappa(method = "gwet")`). Fits from estimators that do not yet expose
#' IFs (FIML, quadratic, continuous, counts) are not supported.
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
         "Only categorical raw estimators (available, ipw, gwet) ",
         "currently support joint_vcov().")

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
