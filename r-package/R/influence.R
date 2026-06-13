#' Per-subject influence functions for misskappa estimates
#'
#' @description
#' `stats::influence()` returns the same per-subject influence-function matrix
#' stored in `fit$psi`. Rows are subjects and columns are coefficients, with
#' scaling
#' \deqn{\mathrm{vcov}(\hat\theta) = n^{-2} \Psi^\top \Psi.}
#'
#' This is the public S3 adapter for R's `stats::influence()` generic; `fit$psi`
#' remains the stored component used internally by `joint_vcov()` and the
#' equality-test helpers.
#'
#' @param model A `misskappa_estimate` object.
#' @param ... Unused; present for S3 generic conformance.
#'
#' @return A numeric matrix of per-subject influence-function rows.
#'
#' @examples
#' fit <- kappa(dat.gwet2014, estimator = "ipw")
#' psi <- stats::influence(fit)
#' all.equal(vcov(fit), crossprod(psi) / nrow(psi)^2)
#'
#' @importFrom stats influence
#' @export
influence.misskappa_estimate <- function(model, ...) {
  psi <- model$psi
  if (is.null(psi)) {
    stop("This misskappa_estimate does not carry per-subject influence functions.")
  }
  psi
}
