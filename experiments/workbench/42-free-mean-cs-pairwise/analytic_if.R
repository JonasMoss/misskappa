# Analytic robust influence function for the free-mean CS pairwise moment
# estimator (the package version of the probe-42 estimator, replacing the
# bootstrap SE). Validated here against the bootstrap before porting to
# r-package/R/kappa_pairwise_cs.R.
#
# Stacked M-estimator, theta = (mu_1..mu_R, psi, e). Per-subject estimating
# functions:
#   g_mu_j = obs_ij (Y_ij - mu_j)
#   g_psi  = sum_{j<k obs} r_ij r_ik            - psi P_i
#   g_e    = sum_{j<k obs} (r_ij - r_ik)^2 / 2  - e   P_i
# At the truth E[r]=0 so E[d g_psi/d mu] = E[d g_e/d mu] = 0: A is block
# diagonal, A_sum = -diag(n_1..n_R, Ptot, Ptot), giving the closed-form IFs
#   IF_mu_j = g_mu_j / n_j,  IF_psi = g_psi / Ptot,  IF_e = g_e / Ptot.
# t3 = sum_j (mu_j - mubar)^2 enters the coefficients with d t3/d mu_j =
# 2(mu_j - mubar) (the O(1/n) bias correction does not affect the first-order
# IF). vcov(theta_hat) = sum_i IF_i IF_i'.

# Requires cs_free_pairwise() and .cs_bp_chance_raw() from free_cs_pairwise.R
# (source that first).

cs_free_pairwise_if <- function(X, values = NULL, equal_mean = FALSE) {
  X <- as.matrix(X); storage.mode(X) <- "double"
  obs <- is.finite(X)
  keep <- rowSums(obs) > 0
  X <- X[keep, , drop = FALSE]; obs <- obs[keep, , drop = FALSE]
  n <- nrow(X); R <- ncol(X)
  n_j <- colSums(obs)
  if (any(n_j == 0)) stop("every rater must be observed at least once.")

  fit <- cs_free_pairwise(X, values, equal_mean)
  mu <- fit$mu; psi <- fit$psi; e <- fit$e; t3 <- fit$t3
  Xc <- X - matrix(mu, n, R, byrow = TRUE)

  # per-subject IF rows for mu (column-mean equations)
  IF_mu <- matrix(0, n, R)
  for (j in seq_len(R)) IF_mu[, j] <- ifelse(obs[, j], Xc[, j], 0) / n_j[j]

  # per-subject estimating functions for psi, e
  g_psi <- numeric(n); g_e <- numeric(n); P <- numeric(n)
  for (i in seq_len(n)) {
    o <- which(obs[i, ]); m <- length(o)
    if (m < 2L) next
    ri <- Xc[i, o]; s1 <- sum(ri); s2 <- sum(ri^2); Pi <- m * (m - 1) / 2
    P[i] <- Pi
    g_psi[i] <- (s1^2 - s2) / 2 - psi * Pi
    g_e[i]   <- (m * s2 - s1^2) / 2 - e * Pi
  }
  Ptot <- sum(P)
  IF_psi <- g_psi / Ptot
  IF_e   <- g_e   / Ptot

  # coefficient gradients wrt (psi, e, t3)
  bp_chance <- .cs_bp_chance_raw(X, values)
  Rm1 <- R - 1
  D <- psi + e + t3 / Rm1
  dC <- c(psi = (e + t3 / Rm1) / D^2, e = -psi / D^2, t3 = -psi / (D^2 * Rm1))
  NF <- R * Rm1 * psi - t3; DF <- Rm1 * (R * (psi + e) + t3)
  dF <- c(psi = R * Rm1 * (DF - NF) / DF^2, e = -R * Rm1 * NF / DF^2,
          t3 = (-DF - NF * Rm1) / DF^2)
  dBP <- c(psi = 0, e = -2 / bp_chance, t3 = -2 / (Rm1 * bp_chance))

  mubar <- mean(mu)
  dt3_dmu <- if (equal_mean) rep(0, R) else 2 * (mu - mubar)
  IF_t3 <- as.numeric(IF_mu %*% dt3_dmu)

  mk <- function(g) g["psi"] * IF_psi + g["e"] * IF_e + g["t3"] * IF_t3
  PSI <- cbind(Conger = mk(dC), Fleiss = mk(dF), `Brennan-Prediger` = mk(dBP))
  vcov <- crossprod(PSI)                       # sum_i IF_i IF_i'
  list(estimates = fit$estimates, se = sqrt(diag(vcov)), vcov = vcov, psi = PSI)
}
