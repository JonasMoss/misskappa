# Free-mean compound-symmetry (CS) pairwise moment estimator: the "pairwise +
# variance_components" cell that completes the FIML-vs-moment grid (the moment
# foil of the free-rater-mean CS FIML, estimand = Conger's quadratic kappa).
#
# Model: Y_ij = mu_j + b_i + eps_ij, Var(b_i)=psi, Var(eps_ij)=e, so
# Sigma = psi J + e I (CS) with FREE rater means mu_j. Raw subjects-by-raters X
# with NA. The recipe (rater-center, then the counts pairwise machinery on the
# residuals, plus a mean-spread t3 with a bias correction):
#   mu_j   = available column mean of rater j            (the "vector of means")
#   r_ij   = Y_ij - mu_j                                  (rater-centred)
#   e      = sum_i sum_{j<k} (r_ij - r_ik)^2 / (2 sum_i P_i)    P_i = C(m_i,2)
#   psi    = sum_i sum_{j<k} r_ij r_ik / sum_i P_i
#   t3     = sum_j (mu_j - mubar)^2  -  bias              bias = tr(A V), see below
# Coefficients (CS-implied t1=R^2 psi + R e, t2=R(psi+e); .kqf_grad forms):
#   Conger = psi / (psi + e + t3/(R-1))
#   Fleiss = (R(R-1) psi - t3) / ((R-1)(R(psi+e) + t3))
#   BP     = 1 - (2e + 2 t3/(R-1)) / bp_chance
# Setting equal_mean = TRUE collapses mu_j to one grand mean and t3 = 0: this is
# the equal-mean (one-way / Fleiss-Cuzick-family) pairwise, demonstrating the
# nesting by mean-model restriction.

.cs_bp_chance_raw <- function(X, values) {
  if (is.null(values)) {
    o <- X[is.finite(X)]
    return(diff(range(o))^2 / 6)                 # continuous-uniform baseline
  }
  C <- length(values)
  (2 / C^2) * (C * sum(values^2) - sum(values)^2)
}

# t3 sampling bias = tr(A V), A = I - 11'/R, V = Cov(mu_hat):
#   V_jj = (psi+e)/n_j,   V_jk = psi * n_jk / (n_j n_k)   (n_jk = co-observations)
.cs_t3_bias <- function(psi, e, obs) {
  R <- ncol(obs)
  n_j <- colSums(obs)
  n_jk <- crossprod(obs)                          # R x R co-observation counts
  V <- psi * (n_jk / outer(n_j, n_j))
  diag(V) <- (psi + e) / n_j
  A <- diag(R) - matrix(1 / R, R, R)
  sum(diag(A %*% V))                              # tr(A V)
}

cs_free_pairwise <- function(X, values = NULL, equal_mean = FALSE) {
  X <- as.matrix(X); storage.mode(X) <- "double"
  obs <- is.finite(X)
  keep <- rowSums(obs) > 0
  X <- X[keep, , drop = FALSE]; obs <- obs[keep, , drop = FALSE]
  n <- nrow(X); R <- ncol(X)
  n_j <- colSums(obs)
  if (any(n_j == 0)) stop("every rater must be observed at least once.")

  # ---- mean vector + rater-centred residuals ----
  if (equal_mean) {
    mu <- rep(sum(X[obs]) / sum(obs), R)
  } else {
    mu <- colSums(ifelse(obs, X, 0)) / n_j
  }
  Xc <- X - matrix(mu, n, R, byrow = TRUE)

  # ---- pooled within-subject pairwise moments on residuals ----
  num_e <- 0; num_psi <- 0; Ptot <- 0
  for (i in seq_len(n)) {
    o <- which(obs[i, ]); m <- length(o)
    if (m < 2L) next
    ri <- Xc[i, o]; s1 <- sum(ri); s2 <- sum(ri^2)
    num_e   <- num_e   + (m * s2 - s1^2)          # sum_{j<k}(ri_j-ri_k)^2
    num_psi <- num_psi + (s1^2 - s2) / 2          # sum_{j<k} ri_j ri_k
    Ptot    <- Ptot + m * (m - 1) / 2
  }
  if (Ptot <= 0) stop("need >= 1 subject with two or more observed ratings.")
  e   <- num_e / (2 * Ptot)
  psi <- num_psi / Ptot

  # ---- mean-spread term ----
  if (equal_mean) {
    t3 <- 0; t3_naive <- 0; bias <- 0
  } else {
    mubar <- mean(mu)
    t3_naive <- sum((mu - mubar)^2)
    bias <- .cs_t3_bias(psi, e, obs)
    t3 <- max(t3_naive - bias, 0)
  }

  bp_chance <- .cs_bp_chance_raw(X, values)
  conger <- psi / (psi + e + t3 / (R - 1))
  fleiss <- (R * (R - 1) * psi - t3) / ((R - 1) * (R * (psi + e) + t3))
  d_obs  <- 2 * e + 2 * t3 / (R - 1)
  bp     <- 1 - d_obs / bp_chance

  list(estimates = c(Conger = conger, Fleiss = fleiss, `Brennan-Prediger` = bp),
       mu = mu, psi = psi, e = e, t3 = t3, t3_naive = t3_naive, bias = bias,
       n = n, R = R)
}

# subject-level nonparametric bootstrap SE (probe-grade; the package version
# would reuse the analytic raw-pairwise influence functions instead)
cs_free_pairwise_boot <- function(X, values = NULL, equal_mean = FALSE, B = 300) {
  X <- as.matrix(X); n <- nrow(X)
  est <- cs_free_pairwise(X, values, equal_mean)$estimates
  bs <- matrix(NA_real_, B, 3)
  for (b in seq_len(B)) {
    fit <- tryCatch(cs_free_pairwise(X[sample(n, n, TRUE), , drop = FALSE],
                                     values, equal_mean)$estimates,
                    error = function(e) rep(NA_real_, 3))
    bs[b, ] <- fit
  }
  list(estimates = est, se = apply(bs, 2, sd, na.rm = TRUE))
}
