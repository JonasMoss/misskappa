#' Unified kappa interface
#'
#' @description
#' Convenience wrapper that dispatches to one of `kappa_raw()`, `kappa_continuous()`,
#' or `kappa_counts()` depending on `type`.
#'
#' @param x A matrix/data frame of ratings. Interpretation depends on `type`.
#' @param type Which data type to use. If `"auto"`, dispatch uses a simple heuristic:
#'   - If any non-integer values are present, uses `"continuous"`.
#'   - Otherwise uses `"raw"`, unless the matrix is non-missing, nonnegative integer
#'     and either `r` is supplied (via `...`) or all row sums are equal, in which case
#'     it uses `"counts"`.
#' @param ... Passed through to the underlying `kappa_*()` function.
#'
#' @return An object of class `misskappa_estimate`.
#' @export
kappa <- function(x, type = c("auto", "raw", "continuous", "counts"), ...) {
  type <- match.arg(type)

  if (!is.matrix(x) && !is.data.frame(x)) stop("'x' must be a matrix or data frame.")
  x_mat <- as.matrix(x)
  if (!is.numeric(x_mat)) stop("'x' must be numeric.")

  dots <- list(...)
  if (type == "auto") {
    type <- guess_kappa_type(x_mat, dots)
  }

  if (type == "raw") res <- kappa_raw(x_mat, ...)
  else if (type == "continuous") res <- kappa_continuous(x_mat, ...)
  else if (type == "counts") res <- kappa_counts(x_mat, ...)
  else stop("Unknown type: ", type)

  attr(res, "call") <- match.call()
  res
}

guess_kappa_type <- function(x_mat, dots, tol = 1e-8) {
  finite <- is.finite(x_mat)
  has_na <- any(is.na(x_mat))
  x_finite <- x_mat[finite]

  is_integer_like <- all(abs(x_finite - round(x_finite)) <= tol)
  if (!is_integer_like) return("continuous")

  if (has_na) return("raw")

  if (!is.null(dots$r)) return("counts")

  is_nonnegative <- all(x_finite >= 0)
  if (!is_nonnegative) return("raw")

  rs <- rowSums(x_mat)
  if (length(rs) > 0 && all(rs == rs[1])) return("counts")

  "raw"
}
