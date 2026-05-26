# Internal constructor for misskappa result objects.
new_misskappa_estimate <- function(estimates,
                                  vcov,
                                  type,
                                  method,
                                  weight = NULL,
                                  call = NULL) {
  if (is.null(call)) call <- sys.call(-1)

  res <- list(estimates = estimates, vcov = vcov)
  attr(res, "type") <- type
  attr(res, "method") <- method
  attr(res, "weight") <- weight
  attr(res, "call") <- call
  class(res) <- "misskappa_estimate"
  res
}

#' @export
print.misskappa_estimate <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  type <- attr(x, "type", exact = TRUE)
  method <- attr(x, "method", exact = TRUE)
  weight <- attr(x, "weight", exact = TRUE)

  header <- sprintf(
    "misskappa (%s): method=%s%s",
    type,
    method,
    if (!is.null(weight)) paste0(", weight=", weight) else ""
  )
  cat(header, "\n", sep = "")

  est <- x$estimates
  if (is.null(names(est))) names(est) <- paste0("V", seq_along(est))

  se <- rep(NA_real_, length(est))
  if (is.matrix(x$vcov) && nrow(x$vcov) == length(est) && ncol(x$vcov) == length(est)) {
    se <- sqrt(diag(x$vcov))
  }

  tab <- data.frame(
    estimate = as.numeric(est),
    se = as.numeric(se),
    row.names = names(est),
    check.names = FALSE
  )

  print(signif(tab, digits = digits))
  invisible(x)
}

#' @export
as.data.frame.misskappa_estimate <- function(x, ...) {
  type <- attr(x, "type", exact = TRUE)
  method <- attr(x, "method", exact = TRUE)
  weight <- attr(x, "weight", exact = TRUE)

  est <- x$estimates
  if (is.null(names(est))) names(est) <- paste0("V", seq_along(est))

  se <- rep(NA_real_, length(est))
  if (is.matrix(x$vcov) && nrow(x$vcov) == length(est) && ncol(x$vcov) == length(est)) {
    se <- sqrt(diag(x$vcov))
  }

  data.frame(
    coefficient = names(est),
    estimate = as.numeric(est),
    se = as.numeric(se),
    type = type,
    method = method,
    weight = if (is.null(weight)) NA_character_ else weight,
    row.names = NULL,
    check.names = FALSE
  )
}

#' @export
coef.misskappa_estimate <- function(object, ...) {
  object$estimates
}

#' @export
vcov.misskappa_estimate <- function(object, ...) {
  object$vcov
}
