#' Calculate limits of a confidence interval.
#'
#' @param alternative Alternative choosen.
#' @param conf_level Confidence level.
#' @keywords internal
limits <- function(alternative, conf_level) {
  half <- (1 - conf_level) / 2
  if (alternative == "two.sided") {
    return(c(half, 1 - half))
  }
  if (alternative == "greater") {
    return(c(2 * half, 1))
  }
  if (alternative == "less") {
    return(c(0, conf_level))
  }
}

#' Half-vectorize matrix.
#'
#' @param x Matrix to vectorize.
#' @keywords internal
vech <- function(x) x[row(x) >= col(x)]

#' Get the indices of the diagonal for the variances.
#' @param r Number of raters.
#' @return A vector of `TRUE` where the variances are.
#' @keywords internal
get_diag_indices <- function(r) {
  indices <- rep(0, choose(r + 1, 2))
  indices[c(1, 1 + cumsum(r:2))] <- 1
  indices
}

#' Gamma matrix
#'
#' Calculate the gamma matrix from a matrix of observations.
#' @param x A numeric matrix of observations.
#' @param sigma Covariance matrix of the data.
#' @param type One of `adf`, `normal`, `elliptical` or `unbiased`.
#' @return The sample estimate of the gamma matrix.
#' @keywords internal
gamma_est <- function(x, sigma, type = "adf") {
  if (type == "adf") {
    gamma_est_adf(x, sigma)
  } else if (type == "unbiased") {
    gamma_est_unbiased(x, sigma)
  } else {
    gamma <- gamma_est_nt(sigma)
    gamma_est_nt(sigma) * kurtosis_correction(x, type = type)
  }
}

#' Asymptotically distribution free covariance matrix.
#' @param x Data.
#' @param sigma Covariance of the data.
#' @return Estimate of the ADF covariance matrix.
#' @keywords Internal.
gamma_est_adf <- function(x, sigma) {
  i_row <- \(n) unlist(lapply(seq_len(n), seq.int, n))
  i_col <- \(n) rep.int(seq_len(n), times = rev(seq_len(n)))
  rows <- i_row(ncol(x))
  cols <- i_col(ncol(x))
  y <- t(x) - colMeans(x, na.rm = TRUE)
  z <- y[cols, ] * y[rows, ]
  mat <- z - rowMeans(z, na.rm = TRUE)

  if (!anyNA(mat)) {
    div <- nrow(x)
  } else {
    nas <- is.na(mat)
    mat[nas] <- 0
    div <- tcrossprod(!nas)
  }
  tcrossprod(mat) / div
}

#' Normal theory gamma matrix
#'
#' Code obtained from `lavaan`:
#' https://github.com/yrosseel/lavaan/blob/6f047c800206d23f246d484b9522295257614222/R/lav_matrix.R
#'
#' Calculate the gamma matrix from a matrix of observations.
#' @param sigma Covariance matrix of the data.
#' @return Normal theory gamma matrix.
#' @keywords internal
gamma_est_nt <- function(sigma) {
  n <- ncol(sigma)

  lower <- lower_vec_indices(n)
  upper <- upper_vec_indices(n)

  y <- sigma %x% sigma
  out <- (y[lower, , drop = FALSE] + y[upper, , drop = FALSE]) / 2
  out[, lower, drop = FALSE] + out[, upper, drop = FALSE]
}

#' Unbiased asymptotic covariance matrix.
#'
#' @param x Data.
#' @param sigma Covariance matrix of the data.
#' @return Unbiased asymptotic covariance matrix.
#' @keywords internal
gamma_est_unbiased <- function(x, sigma) {
  gamma_adf <- gamma_est_adf(x, sigma)
  gamma_nt <- gamma_est_nt(sigma)
  gamma_rem <- tcrossprod(vech(sigma))
  n <- nrow(x)
  mult <- n / ((n - 2) * (n - 3))
  mult * ((n - 1) * gamma_adf - (gamma_nt - 2 / (n - 1) * gamma_rem))
}


#' Obtain indices of lower or upper triangular matrix in vec indices.
#'
#' Code obtained from `lavaan`:
#' https://github.com/yrosseel/lavaan/blob/6f047c800206d23f246d484b9522295257614222/R/lav_matrix.R
#'
#' @param n Dimension of square matrix.
#' @param diagonal If `TRUE`, includes the diagonal elements.
#' @returns Indices `x` so that `a[x] = c(a)[x]` returns the elements
#'    of the lower (upper) diagonal matrix in row-wise (column-wise)
#'    order.
#' @keywords internal
#'
lower_vec_indices <- function(n = 1L, diagonal = TRUE) {
  rows <- matrix(seq_len(n), n, n)
  cols <- matrix(seq_len(n), n, n, byrow = TRUE)
  if (diagonal) which(rows >= cols) else which(rows > cols)
}

upper_vec_indices <- function(n = 1L, diagonal = TRUE) {
  rows <- matrix(seq_len(n), n, n)
  cols <- matrix(seq_len(n), n, n, byrow = TRUE)
  tmp <- matrix(seq_len(n * n), n, n, byrow = TRUE)
  if (diagonal) tmp[rows >= cols] else tmp[rows > cols]
}

#' Calculate unbiased sample kurtosis.
#' @param x Matrix of valus.
#' @return Unbiased sample kurtosis.
#' @keywords internal
kurtosis <- function(x) {
  n <- nrow(x)
  g2 <- \(x) mean((x - mean(x, na.rm = TRUE))^4) / stats::var(x, na.rm = TRUE)^2
  kurtosis <- \(x) (n - 1) / ((n - 2) * (n - 3)) * ((n + 1) * g2(x) + 6)
  mean(apply(x, 2, kurtosis), na.rm = TRUE) - 3
}

#' Calculate kurtosis correction
#' @param x Matrix of values
#' @param type The type of correction, either "normal" or "elliptical".
#' @keywords internal
kurtosis_correction <- function(x, type) {
  kurt <- if (type == "normal") 0 else kurtosis(x)
  1 + kurt / 3
}

#' Calculates the empirical capital pi matrix.
#'
#' @param x Vector of probabilities for being missing.
#' @return The capital pi matrix.
#' @keywords internal
pi_mat_empirical <- \(x) {
  r <- ncol(x)
  if (!anyNA(x)) {
    return(matrix(1, choose(r + 1, 2), choose(r + 1, 2)))
  }

  ind2 <- arrangements::combinations(r, 2, replace = TRUE)
  ind4 <- arrangements::combinations(seq_len(nrow(ind2)), 2, replace = TRUE)

  nisna <- !is.na(x)
  combs <- apply(ind2, 1, \(i) nisna[, i[1]] & nisna[, i[2]])
  p2_hats <- colMeans(combs)

  hats <- apply(ind4, 1, \(i) {
    mean(combs[, i[1]] & combs[, i[2]]) / (p2_hats[i[1]] * p2_hats[i[2]])
  })

  new_mat <- matrix(NA, choose(r + 1, 2), choose(r + 1, 2))
  new_mat[ind4] <- hats
  as.matrix(Matrix::forceSymmetric(new_mat, uplo = "U"))
}

#' Get the required value of c1.
#' @param values Vector of values.
#' @param kind The kind of c1 requested.
#' @return The value of c1.
#' @keywords internal
bp_get_c1 <- \(values, kind) {
  if (kind == 1) {
    n_cat <- length(values)
    (2 * n_cat * sum(values^2) - 2 * sum(values)^2) / n_cat^2
  } else {
    0.5 * (max(values) - min(values))^2
  }
}

#' @export
print.quadagree <- function(x, digits = getOption("digits"), ...) {
  at <- \(y) attr(x, y)
  cat("Call: ", paste(deparse(at("call")),
    sep = "\n",
    collapse = "\n"
  ), "\n\n", sep = "")

  if (!is.null(x)) {
    cat(format(100 * at("conf_level")),
      "% confidence interval (n = ", at("n"), ").\n",
      sep = ""
    )
    print(x[1:2], digits = digits)
    cat("\n")
  }

  if (!is.null(at("estimate"))) {
    cat("Sample estimates.\n")
    print(
      c(
        kappa = at("estimate"),
        sd = at("sd")
      ),
      digits = digits
    )
  }
  invisible(x)
}

get_transformer <- function(transform) {
  transformers <- list(
    fisher = transformer_fisher,
    none = transformer_none,
    arcsin = transformer_arcsin,
    log = transformer_log
  )

  if (transform %in% names(transformers)) {
    transformers[[transform]]
  } else {
    stop(paste0("`transformer = ", transform, "` not supported."))
  }
}

transformer_fisher <- c(
  est = \(est) atanh(est),
  sd = \(est, sd) sd / (1 - est^2),
  inv = tanh
)

transformer_log <- c(
  est = \(est) log(1 - est),
  sd = \(est, sd) sd / abs(1 - est),
  inv = \(x) 1 - exp(x)
)

transformer_none <- c(
  est = \(est) est,
  sd = \(est, sd) sd,
  inv = \(x) x
)

transformer_arcsin <- c(
  est = asin,
  sd = \(est, sd) sd / sqrt(1 - est^2),
  inv = sin
)

tr <- \(x) {
  sum(x[1L + 0L:(dim(x)[1L] - 1L) * (dim(x)[1L] + 1L)])
}
