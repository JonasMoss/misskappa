# Free-mean compound-symmetry (CS) pairwise moment estimator for the quadratic
# kappas. This is the distribution-free MCAR moment companion to the
# compound-symmetry normal FIML (`estimator = "nt_fiml_cs"`): the
# "pairwise + compound-symmetry" cell of the moment-vs-FIML grid.
#
# Model Y_ij = mu_j + b_i + eps_ij, Var(b_i) = psi, Var(eps_ij) = e, so
# Sigma = psi J + e I (compound symmetry) with FREE rater means mu_j. Recipe:
# rater-center by the available column means, then the same pooled
# within-subject pairwise moments as the counts/exchangeable case, plus a
# mean-spread term t3 with a finite-sample bias correction.
#
#   mu_j = available column mean of rater j
#   r_ij = Y_ij - mu_j
#   e    = sum_i sum_{j<k obs} (r_ij - r_ik)^2 / (2 sum_i P_i),  P_i = C(m_i, 2)
#   psi  = sum_i sum_{j<k obs} r_ij r_ik / sum_i P_i
#   t3   = sum_j (mu_j - mubar)^2 - tr(A V),  A = I - 11'/R, V = Cov(mu_hat)
#
# Coefficients (CS-implied t1 = R^2 psi + R e, t2 = R(psi + e)):
#   Conger = psi / (psi + e + t3/(R-1))
#   Fleiss = (R(R-1) psi - t3) / ((R-1)(R(psi+e) + t3))
#   BP     = 1 - (2e + 2 t3/(R-1)) / bp_chance
#
# Inference is the analytic i.i.d.-subjects sandwich. The estimator is a stacked
# M-estimator in theta = (mu_1..mu_R, psi, e) whose expected cross-Jacobian
# d g_{psi,e} / d mu vanishes at the truth (E[r] = 0), so A is block diagonal and
# the per-subject influence is closed form:
#   IF_mu_j = obs_ij (Y_ij - mu_j) / n_j
#   IF_psi  = (sum_{j<k} r_ij r_ik         - psi P_i) / Ptot
#   IF_e    = (sum_{j<k} (r_ij - r_ik)^2/2 - e   P_i) / Ptot
# t3 enters the coefficients with d t3 / d mu_j = 2(mu_j - mubar); the O(1/n)
# bias correction does not affect the first-order influence. bp_chance is a fixed
# function of the rating support (treated as constant, as in the other quadratic
# paths). Validated against the subject-level bootstrap in
# experiments/workbench/42-free-mean-cs-pairwise/.

# t3 sampling bias = tr(A V), A = I - 11'/R, V = Cov(mu_hat):
#   V_jj = (psi + e)/n_j,   V_jk = psi * n_jk / (n_j n_k)   (n_jk = co-observations)
.cs_pairwise_t3_bias <- function(psi, e, obs) {
  R <- ncol(obs)
  n_j <- colSums(obs)
  n_jk <- crossprod(obs)                        # R x R co-observation counts
  V <- psi * (n_jk / outer(n_j, n_j))
  diag(V) <- (psi + e) / n_j
  A <- diag(R) - matrix(1 / R, R, R)
  sum(diag(A %*% V))
}

# Coefficient gradients d(Conger, Fleiss, Brennan-Prediger) / d(psi, e, t3).
.cs_pairwise_kappa_jac <- function(psi, e, t3, R, bp_chance) {
  Rm1 <- R - 1
  D <- psi + e + t3 / Rm1
  NF <- R * Rm1 * psi - t3
  DF <- Rm1 * (R * (psi + e) + t3)
  rbind(
    Conger = c(psi = (e + t3 / Rm1) / D^2, e = -psi / D^2, t3 = -psi / (D^2 * Rm1)),
    Fleiss = c(psi = R * Rm1 * (DF - NF) / DF^2, e = -R * Rm1 * NF / DF^2,
               t3 = (-DF - NF * Rm1) / DF^2),
    `Brennan-Prediger` = c(psi = 0, e = -2 / bp_chance, t3 = -2 / (Rm1 * bp_chance))
  )
}

# ---- entry point (called from kappa() for estimator = "pairwise_cs") ---------

.kappa_pairwise_cs <- function(x_mat, values = NULL) {
  X <- as.matrix(x_mat); storage.mode(X) <- "double"
  obs <- is.finite(X)
  keep <- rowSums(obs) > 0
  X <- X[keep, , drop = FALSE]; obs <- obs[keep, , drop = FALSE]
  n <- nrow(X); R <- ncol(X)
  if (R < 2L) stop("quadratic kappa requires at least two raters.", call. = FALSE)
  n_j <- colSums(obs)
  if (any(n_j == 0)) {
    stop("every rater must be observed by at least one subject.", call. = FALSE)
  }

  # mean vector and rater-centered residuals
  mu <- colSums(ifelse(obs, X, 0)) / n_j
  Xc <- X - matrix(mu, n, R, byrow = TRUE)

  # pooled within-subject pairwise moments, plus per-subject scores for the IF
  IF_mu <- matrix(0, n, R)
  for (j in seq_len(R)) IF_mu[, j] <- ifelse(obs[, j], Xc[, j], 0) / n_j[j]
  g_psi <- numeric(n); g_e <- numeric(n); P <- numeric(n)
  num_psi <- 0; num_e <- 0
  for (i in seq_len(n)) {
    o <- which(obs[i, ]); m <- length(o)
    if (m < 2L) next
    ri <- Xc[i, o]; s1 <- sum(ri); s2 <- sum(ri^2); Pi <- m * (m - 1) / 2
    P[i] <- Pi
    g_psi[i] <- (s1^2 - s2) / 2          # sum_{j<k} r_ij r_ik (psi P_i subtracted below)
    g_e[i]   <- (m * s2 - s1^2) / 2      # sum_{j<k} (r_ij - r_ik)^2 / 2
    num_psi <- num_psi + g_psi[i]
    num_e   <- num_e + g_e[i]
  }
  Ptot <- sum(P)
  if (Ptot <= 0) {
    stop("need at least one subject with two or more observed ratings.",
         call. = FALSE)
  }
  psi <- num_psi / Ptot
  e   <- num_e / Ptot

  # mean-spread term with finite-sample bias correction
  mubar <- mean(mu)
  t3 <- max(sum((mu - mubar)^2) - .cs_pairwise_t3_bias(psi, e, obs), 0)

  bp_chance <- .quadratic_bp_chance(X, values)
  estimates <- c(
    Conger = psi / (psi + e + t3 / (R - 1)),
    Fleiss = (R * (R - 1) * psi - t3) / ((R - 1) * (R * (psi + e) + t3)),
    `Brennan-Prediger` = 1 - (2 * e + 2 * t3 / (R - 1)) / bp_chance
  )

  # per-subject influence rows (theta_hat - theta ~ sum_i IF_i)
  IF_psi <- (g_psi - psi * P) / Ptot
  IF_e   <- (g_e - e * P) / Ptot
  IF_t3  <- as.numeric(IF_mu %*% (2 * (mu - mubar)))
  J <- .cs_pairwise_kappa_jac(psi, e, t3, R, bp_chance)
  IF <- cbind(IF_psi, IF_e, IF_t3) %*% t(J)         # n x 3, columns Conger/Fleiss/BP
  colnames(IF) <- names(estimates)

  psi_mat <- n * IF                                  # package convention: vcov = crossprod(psi)/n^2
  vcov_mat <- crossprod(psi_mat) / n^2
  dimnames(vcov_mat) <- list(names(estimates), names(estimates))

  structure(
    list(
      estimates = estimates,
      vcov = vcov_mat,
      psi = psi_mat,
      method = "pairwise_cs",
      weight = "quadratic",
      bp_chance = bp_chance,
      moments = list(mu = mu, psi = psi, e = e, t3 = t3, n = n, R = R)
    ),
    class = "misskappa_estimate"
  )
}
