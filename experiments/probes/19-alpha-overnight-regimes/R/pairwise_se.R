# Pairwise-available coefficient alpha: point estimate, gradient, and the
# overlapping-subsample asymptotic standard error.
#
# This is the PROTOTYPE of the paper's centerpiece SE (van Praag, Dijkstra &
# van Velzen 1985). The asymptotic covariance of the pairwise-deleted
# covariance entries is, under MCAR,
#
#   Cov(sigma_hat_{ab}, sigma_hat_{cd})
#       = ( n_{abcd} / (n_{ab} * n_{cd}) ) * ( mu4_{abcd} - sigma_{ab} sigma_{cd} ),
#
# where n_{ab} is the number observing both a,b; n_{abcd} the number observing
# all four; and mu4 the (centered) fourth cross-moment over that four-way
# overlap. Contracting this with the gradient of alpha(Sigma) gives the SE.
# When two entries share no respondents (n_{abcd} = 0) their estimates are
# asymptotically independent, so the contribution is 0 -- which also makes the
# routine degrade gracefully under zero-overlap planned designs.
#
# Validated against the case bootstrap inside this experiment; once it tracks,
# it ports to misskappa (joint_vcov()). Centering uses available-case means
# (consistent under MCAR, asymptotically equivalent to the pair-specific means
# R's cov(use = "pairwise") applies for the point estimate).

# Unique lower-triangle (vech, incl. diagonal) indices in column-major order;
# the SAME ordering feeds the gradient and the covariance so they line up.
amc_vech_index <- function(p) {
  which(lower.tri(matrix(0, p, p), diag = TRUE), arr.ind = TRUE)
}

alpha_point_from_cov <- function(S) {
  p <- ncol(S)
  (p / (p - 1)) * (1 - sum(diag(S)) / sum(S))
}

# Gradient of alpha wrt the UNIQUE covariance parameters vech(Sigma).
# alpha = c*(1 - t/s), c = p/(p-1), t = tr(Sigma), s = 1' Sigma 1.
# dalpha/dSigma = c*( -(1/s) I + (t/s^2) J ); off-diagonal unique params appear
# twice in the symmetric Sigma, so their derivative doubles.
alpha_grad_vech <- function(S) {
  p <- ncol(S)
  s <- sum(S); t <- sum(diag(S)); cc <- p / (p - 1)
  G <- cc * (-(1 / s) * diag(p) + (t / s^2) * matrix(1, p, p))
  idx <- amc_vech_index(p)
  vapply(seq_len(nrow(idx)), function(r) {
    a <- idx[r, 1]; b <- idx[r, 2]
    if (a == b) G[a, b] else G[a, b] + G[b, a]
  }, numeric(1))
}

# Optional self-check: analytic vech gradient vs numDeriv on alpha(theta).
alpha_grad_check <- function(S) {
  if (!requireNamespace("numDeriv", quietly = TRUE)) return(NA_real_)
  p <- ncol(S); idx <- amc_vech_index(p)
  f <- function(theta) {
    M <- matrix(0, p, p)
    M[idx] <- theta
    M[upper.tri(M)] <- t(M)[upper.tri(M)]
    alpha_point_from_cov(M)
  }
  g_num <- numDeriv::grad(f, S[idx])
  max(abs(g_num - alpha_grad_vech(S)))
}

# Asymptotic covariance of vech(Sigma_hat^PA) via the four-way overlap moments.
pairwise_vech_vcov <- function(X) {
  p <- ncol(X)
  R <- !is.na(X)
  idx <- amc_vech_index(p)
  m <- nrow(idx)
  mu <- colMeans(X, na.rm = TRUE)
  Xc <- sweep(X, 2L, mu, "-")
  Xc[is.na(Xc)] <- 0  # zeroed; masks below exclude unobserved rows from sums

  # pair counts and pairwise (centered) covariances, available-case means
  n_pair <- matrix(0L, m, 1L)
  s_pair <- numeric(m)
  for (r in seq_len(m)) {
    a <- idx[r, 1]; b <- idx[r, 2]
    obs <- R[, a] & R[, b]
    nab <- sum(obs)
    n_pair[r] <- nab
    s_pair[r] <- if (nab > 1L) sum(Xc[obs, a] * Xc[obs, b]) / (nab - 1L) else NA_real_
  }

  V <- matrix(0, m, m)
  for (r in seq_len(m)) {
    if (!is.finite(s_pair[r])) { V[r, ] <- NA; V[, r] <- NA; next }
    a <- idx[r, 1]; b <- idx[r, 2]
    for (q in r:m) {
      if (!is.finite(s_pair[q])) { V[r, q] <- V[q, r] <- NA; next }
      cc <- idx[q, 1]; dd <- idx[q, 2]
      both <- R[, a] & R[, b] & R[, cc] & R[, dd]
      nabcd <- sum(both)
      if (nabcd == 0L) { V[r, q] <- V[q, r] <- 0; next }
      mu4 <- sum(Xc[both, a] * Xc[both, b] * Xc[both, cc] * Xc[both, dd]) / nabcd
      cov_e <- (nabcd / (n_pair[r] * n_pair[q])) * (mu4 - s_pair[r] * s_pair[q])
      V[r, q] <- V[q, r] <- cov_e
    }
  }
  list(vcov = V, n_pair = drop(n_pair), s_pair = s_pair, idx = idx)
}

# Pairwise alpha point estimate + delta-method SE. Returns NA SE (and NA
# estimate) when the pairwise covariance is undefined (some pair never
# co-observed) or the variance is non-finite/negative.
pairwise_alpha <- function(X) {
  p <- ncol(X)
  S <- stats::cov(X, use = "pairwise.complete.obs")
  est <- alpha_point_from_cov(S)
  undefined <- !is.finite(est)
  se <- NA_real_
  if (!undefined) {
    vc <- tryCatch(pairwise_vech_vcov(X), error = function(e) NULL)
    if (!is.null(vc) && all(is.finite(vc$vcov))) {
      g <- alpha_grad_vech(S)
      v <- as.numeric(t(g) %*% vc$vcov %*% g)
      if (is.finite(v) && v > 0) se <- sqrt(v)
    }
  }
  npd <- if (undefined) NA else {
    ev <- tryCatch(min(eigen(S, only.values = TRUE)$values), error = function(e) NA_real_)
    is.finite(ev) && ev <= 0
  }
  list(estimate = if (undefined) NA_real_ else est, se = se,
       undefined = undefined, npd = npd)
}

# Feldt (1965; Feldt-Woodruff-Salih 1987) confidence interval with a single
# "average n" -- the de-facto software interval for pairwise alpha and the
# strawman this paper's SE is meant to beat. n defaults to the mean off-diagonal
# pairwise sample size.
feldt_alpha <- function(X, level = 0.95) {
  p <- ncol(X)
  S <- stats::cov(X, use = "pairwise.complete.obs")
  est <- alpha_point_from_cov(S)
  R <- !is.na(X)
  npair <- numeric(0)
  for (a in 1:(p - 1)) for (b in (a + 1):p) npair <- c(npair, sum(R[, a] & R[, b]))
  n_avg <- mean(npair[npair > 0])
  if (!is.finite(est) || !is.finite(n_avg) || n_avg < 3) {
    return(list(estimate = est, se = NA_real_, lwr = NA_real_, upr = NA_real_, n_avg = n_avg))
  }
  g <- 1 - level
  df1 <- n_avg - 1; df2 <- (n_avg - 1) * (p - 1)
  lwr <- 1 - (1 - est) * stats::qf(1 - g / 2, df1, df2)
  upr <- 1 - (1 - est) * stats::qf(g / 2, df1, df2)
  list(estimate = est, se = (upr - lwr) / (2 * stats::qnorm(level + (1 - level) / 2)),
       lwr = lwr, upr = upr, n_avg = n_avg)
}
