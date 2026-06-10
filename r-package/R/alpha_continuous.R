# Coefficient alpha for continuous items under ignorable missingness, via a
# saturated normal FIML / EM covariance and a sandwich delta-method standard
# error.
#
# This is the {continuous, missing} cell of the alpha-under-missingness grid:
# the categorical multinomial-EM route lives in alpha_cat_fiml() /
# rcpp_alpha_raw(). The fit here is pure R (a covariance-level routine, not an
# agreement-coefficient kernel), validated against magmaan's
# estimate_saturated_em_moments() and lavaan's saturated h1 estimator for the
# moments, and against
# coefficientalpha (varphi = 0) and a case bootstrap for the standard error.
#
# Conventions kept deliberately aligned with magmaan so the two implementations
# cross-check: the MLE covariance divides by n, and vech() is the column-major
# lower triangle (including the diagonal), the same order alpha_grad_s() weights.

# ---- vech bookkeeping -------------------------------------------------------

# Column-major lower-triangular half-vectorisation, matching magmaan's
# `for (j) for (i in j:p)` order.
.amc_vech <- function(S) S[lower.tri(S, diag = TRUE)]

# p x p integer matrix giving, for i >= j, the position of (i, j) in vech();
# symmetric for convenient lookup.
.amc_vech_pos <- function(p) {
  pos <- matrix(0L, p, p)
  idx <- which(lower.tri(pos, diag = TRUE), arr.ind = TRUE) # column-major
  for (k in seq_len(nrow(idx))) {
    i <- idx[k, 1L]; j <- idx[k, 2L]
    pos[i, j] <- k
    pos[j, i] <- k
  }
  pos
}

# Cronbach's alpha from a covariance matrix.
.amc_alpha_from_S <- function(S) {
  p <- ncol(S)
  (p / (p - 1)) * (1 - sum(diag(S)) / sum(S))
}

# Gradient of alpha wrt vech(S). Diagonal entries feed tr(S) (weight a),
# off-diagonals feed 1'S1 with weight 2 (weight b), since the off-diagonal
# vech entry is the single free parameter for the symmetric pair.
.amc_alpha_grad_s <- function(S) {
  p <- ncol(S)
  idx <- which(lower.tri(S, diag = TRUE), arr.ind = TRUE)
  diagonal <- idx[, 1L] == idx[, 2L]
  a <- ifelse(diagonal, 1, 0)
  b <- ifelse(diagonal, 1, 2)
  sig <- .amc_vech(S)
  as <- sum(a * sig)
  bs <- sum(b * sig)
  -(p / (p - 1)) * (a / bs - as * b / bs^2)
}

# ---- saturated normal EM ----------------------------------------------------

# Missing-data patterns: list of {rows, o} where o is the observed column set.
.amc_patterns <- function(R) {
  key <- apply(R, 1L, function(r) paste0(as.integer(r), collapse = ""))
  lapply(split(seq_len(nrow(R)), key), function(rows) {
    list(rows = rows, o = which(R[rows[1L], ]))
  })
}

# EM for the saturated multivariate-normal mean and covariance (MLE, /n).
.amc_em <- function(X, patterns, tol = 1e-8, max_iter = 10000L) {
  n <- nrow(X); p <- ncol(X)
  mu <- colMeans(X, na.rm = TRUE)
  v <- apply(X, 2L, stats::var, na.rm = TRUE)
  v[!is.finite(v) | v <= 0] <- 1
  Sigma <- diag(v, p)
  converged <- FALSE
  iter <- 0L
  for (iter in seq_len(max_iter)) {
    T1 <- numeric(p)
    T2 <- matrix(0, p, p)
    for (g in patterns) {
      o <- g$o; rows <- g$rows; ng <- length(rows)
      m <- setdiff(seq_len(p), o)
      Xo <- X[rows, o, drop = FALSE]
      Xc <- matrix(0, ng, p)
      Xc[, o] <- Xo
      if (length(m) > 0L) {
        Soo_inv <- solve(Sigma[o, o, drop = FALSE])
        B <- Sigma[m, o, drop = FALSE] %*% Soo_inv          # regress m on o
        Xo_c <- sweep(Xo, 2L, mu[o], "-")
        Xc[, m] <- sweep(Xo_c %*% t(B), 2L, mu[m], "+")     # E[x_m | x_o]
        Cmm <- Sigma[m, m, drop = FALSE] -
          B %*% Sigma[o, m, drop = FALSE]                   # residual cov
        T2[m, m] <- T2[m, m] + ng * Cmm
      }
      T1 <- T1 + colSums(Xc)
      T2 <- T2 + crossprod(Xc)
    }
    mu_new <- T1 / n
    Sigma_new <- T2 / n - tcrossprod(mu_new)
    Sigma_new <- (Sigma_new + t(Sigma_new)) / 2
    delta <- max(abs(mu_new - mu), abs(Sigma_new - Sigma))
    mu <- mu_new
    Sigma <- Sigma_new
    if (delta < tol) { converged <- TRUE; break }
  }
  list(mu = mu, Sigma = Sigma, iterations = iter, converged = converged)
}

# ---- casewise scores --------------------------------------------------------

# n x q matrix of casewise scores of the observed-data Gaussian log-likelihood
# wrt theta = (mu, vech(Sigma)). Analytic; the off-diagonal vech score carries
# the factor 2 of the symmetric pair, matching .amc_alpha_grad_s().
.amc_score_matrix <- function(mu, Sigma, X, patterns, vech_pos, p, pstar) {
  n <- nrow(X)
  S <- matrix(0, n, p + pstar)
  for (g in patterns) {
    o <- g$o; rows <- g$rows
    no <- length(o)
    if (no == 0L) next
    Soinv <- solve(Sigma[o, o, drop = FALSE])
    E <- sweep(X[rows, o, drop = FALSE], 2L, mu[o], "-")
    U <- E %*% Soinv                                        # rows are (Soinv e)'
    S[rows, o] <- S[rows, o] + U                            # d/dmu
    half_const <- 0.5 * Soinv
    for (a in seq_len(no)) {
      for (b in seq_len(a)) {
        wt <- if (a == b) 1 else 2
        col <- p + vech_pos[o[a], o[b]]
        S[rows, col] <- S[rows, col] +
          wt * (0.5 * U[, a] * U[, b] - half_const[a, b])
      }
    }
  }
  S
}

.amc_unpack <- function(theta, p, pstar) {
  Sig <- matrix(0, p, p)
  Sig[lower.tri(Sig, diag = TRUE)] <- theta[(p + 1L):(p + pstar)]
  Sig <- Sig + t(Sig) - diag(diag(Sig))
  list(mu = theta[seq_len(p)], Sigma = Sig)
}

# Observed information per case, H = -(1/n) d(sum_i score_i)/dtheta, by central
# finite differences on the analytic total score. Sidesteps duplication-matrix
# bookkeeping; H and the sandwich are then exact functions of the analytic score.
.amc_information <- function(theta, X, patterns, vech_pos, p, pstar, h = 1e-5) {
  n <- nrow(X); q <- p + pstar
  g_of <- function(th) {
    pr <- .amc_unpack(th, p, pstar)
    colSums(.amc_score_matrix(pr$mu, pr$Sigma, X, patterns, vech_pos, p, pstar))
  }
  H <- matrix(0, q, q)
  for (k in seq_len(q)) {
    tp <- theta; tp[k] <- tp[k] + h
    tm <- theta; tm[k] <- tm[k] - h
    H[, k] <- -(g_of(tp) - g_of(tm)) / (2 * h) / n
  }
  (H + t(H)) / 2
}

# ---- public estimator -------------------------------------------------------

#' Coefficient alpha for continuous items under missing data (normal FIML)
#'
#' @description
#' Backend for `alpha(estimator = "nt_fiml")`. Estimates Cronbach's
#' coefficient alpha for a battery of continuous items with missing entries,
#' via a saturated multivariate-normal covariance fitted by full-information
#' maximum likelihood (the EM algorithm), valid under ignorable (MCAR or MAR)
#' missingness. The standard error is the sandwich delta-method contraction of
#' the alpha gradient with the asymptotic covariance of the fitted moments.
#'
#' The fit and its sandwich covariance reproduce magmaan's
#' `estimate_saturated_em_moments()` and lavaan's saturated estimator; the
#' interval matches `coefficientalpha::alpha(varphi = 0)` and a nonparametric
#' case bootstrap.
#'
#' This backend is normally reached through [alpha()], which checks that every
#' item is observed and every item pair is jointly observed before dispatch.
#' Direct callers should enforce the same fixed-item condition; otherwise the
#' saturated covariance functional is not identified from the observed
#' missing-data pattern.
#'
#' @param x A subjects-by-items numeric matrix or data frame; `NA` marks
#'   missing entries. Rows that are entirely missing are dropped.
#' @param em_options Named list tuning the EM fit: `tol` (convergence
#'   tolerance on the moments, default `1e-8`), `max_iter` (default `10000`),
#'   and `fd_h` (finite-difference step for the information matrix, default
#'   `1e-5`). Pass any subset.
#'
#' @return An object of class `misskappa_estimate` carrying one coefficient
#'   named `alpha` and its asymptotic covariance. Additional fields: `moments`
#'   (the fitted `mu`, `Sigma`, EM `iterations`, and `converged` flag) and
#'   `psi` (per-subject influence-function rows). Methods: `print`, `coef`,
#'   `vcov`, `confint`, `as.data.frame`, and `stats::influence`.
#'
#' @references
#' Zhang, Z., & Yuan, K.-H. (2016). Robust coefficients alpha and omega and
#' confidence intervals with outlying observations and missing data.
#' *Educational and Psychological Measurement*, 76(3), 387-411.
#'
#' @keywords internal
alpha_continuous <- function(x, em_options = list()) {
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a matrix or data frame.")
  }
  X <- as.matrix(x)
  if (!is.numeric(X)) stop("'x' must be numeric.")
  storage.mode(X) <- "double"
  p <- ncol(X)
  if (p < 2L) stop("coefficient alpha requires at least two items.")

  R <- is.finite(X)
  keep <- rowSums(R) > 0L
  if (!all(keep)) {
    X <- X[keep, , drop = FALSE]
    R <- R[keep, , drop = FALSE]
  }
  n <- nrow(X)
  if (any(colSums(R) == 0L)) {
    stop("every item must be observed for at least one subject.")
  }

  opt <- utils::modifyList(
    list(tol = 1e-8, max_iter = 10000L, fd_h = 1e-5), em_options
  )
  pstar <- p * (p + 1L) / 2L
  vech_pos <- .amc_vech_pos(p)
  patterns <- .amc_patterns(R)

  em <- .amc_em(X, patterns, tol = opt$tol, max_iter = opt$max_iter)
  if (!em$converged) {
    warning("saturated EM did not converge in ", opt$max_iter, " iterations.")
  }
  Sigma <- em$Sigma
  alpha_hat <- .amc_alpha_from_S(Sigma)

  theta <- c(em$mu, .amc_vech(Sigma))
  scores <- .amc_score_matrix(em$mu, Sigma, X, patterns, vech_pos, p, pstar)
  H <- .amc_information(theta, X, patterns, vech_pos, p, pstar, h = opt$fd_h)
  Hinv <- solve(H)

  g <- c(numeric(p), .amc_alpha_grad_s(Sigma))            # alpha has no mu part

  estimates <- c(alpha = alpha_hat)
  # Per-subject influence psi_i = g' H^{-1} s_i; vcov = crossprod(psi) / n^2.
  psi <- matrix(scores %*% (Hinv %*% g), ncol = 1L,
                dimnames = list(NULL, "alpha"))
  var_alpha <- sum(psi^2) / n^2
  vcov_mat <- matrix(var_alpha, 1L, 1L, dimnames = list("alpha", "alpha"))

  structure(
    list(
      estimates = estimates,
      vcov = vcov_mat,
      psi = psi,
      method = "alpha-continuous-fiml",
      weight = "score",
      moments = list(mu = em$mu, Sigma = Sigma,
                     iterations = em$iterations, converged = em$converged)
    ),
    class = "misskappa_estimate"
  )
}
