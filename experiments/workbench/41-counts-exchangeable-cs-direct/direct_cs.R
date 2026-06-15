# Direct-on-counts exchangeable compound-symmetry (CS) quadratic kappa.
#
# Model: each subject i has a latent r-vector Y_i ~ N(mu 1, Sigma),
# Sigma = psi J + e I, with psi >= 0, e > 0 (common mean, exchangeable raters).
# A subject contributes only its OBSERVED sub-vector of size m_i; by
# exchangeability that marginal is again CS-normal with the SAME (mu, psi, e),
# so the per-subject likelihood depends on the data only through the sufficient
# statistics
#     m_i = number of ratings,
#     s_i = sum of observed scores      = sum_c n_ic v_c,
#     q_i = sum of squared scores       = sum_c n_ic v_c^2,
# plus, for the moment estimators, the within-subject pair statistics
#     p_i = C(m_i, 2)            (number of within-subject rater pairs),
#     c_i = sum_{j<k} y_ij y_ik = (s_i^2 - q_i) / 2.
# Nothing here needs a realizing subjects-by-raters table or r_total (except
# listwise, which needs r_total only to decide which rows are complete).
#
# Quadratic kappas as functions of (psi, e) (derivation in the dev note):
#     Fleiss            = psi / (psi + e)            ( = the ICC )
#     Brennan-Prediger  = 1 - 2 e / bp_chance
# where bp_chance = (2 / C^2) (C sum v^2 - (sum v)^2) is the package's uniform
# quadratic chance baseline. Both are r_total-invariant.
#
# Every estimator returns a misskappa-shaped list with O(1) per-subject
# influence rows `psi` obeying the package convention vcov = crossprod(psi)/n^2,
# so the SEs are robust (empirical sandwich) and joint_vcov()-compatible.

# ---- sufficient statistics --------------------------------------------------

.cs_suff <- function(counts, values) {
  counts <- as.matrix(counts)
  storage.mode(counts) <- "double"
  m <- rowSums(counts)
  s <- as.numeric(counts %*% values)
  q <- as.numeric(counts %*% (values^2))
  data.frame(m = m, s = s, q = q)
}

.cs_bp_chance <- function(values) {
  C <- length(values)
  (2 / C^2) * (C * sum(values^2) - sum(values)^2)
}

.cs_kappas <- function(psi, e, bp_chance) {
  c(Fleiss = psi / (psi + e),
    `Brennan-Prediger` = 1 - 2 * e / bp_chance)
}

# Jacobian d(Fleiss, BP) / d(mu, psi, e), a 2x3 matrix (mu column is 0).
.cs_kappa_jac <- function(psi, e, bp_chance) {
  tot <- psi + e
  J <- matrix(0, 2, 3, dimnames = list(c("Fleiss", "Brennan-Prediger"),
                                       c("mu", "psi", "e")))
  J["Fleiss", "psi"] <- e / tot^2
  J["Fleiss", "e"]   <- -psi / tot^2
  J["Brennan-Prediger", "e"] <- -2 / bp_chance
  J
}

.cs_estimate <- function(estimates, vcov, psi, method, moments) {
  structure(list(estimates = estimates, vcov = vcov, psi = psi,
                 method = method, weight = "quadratic", moments = moments),
            class = "misskappa_estimate")
}

# ---- per-subject CS log-likelihood (vectorised over subjects) ---------------

# Returns the n-vector of per-subject log densities for the observed sub-vector,
# from the sufficient statistics alone.
.cs_loglik_i <- function(mu, psi, e, m, s, q) {
  A <- q - 2 * mu * s + m * mu^2          # sum_c n_ic (v_c - mu)^2
  B <- s - m * mu                          # 1' (y - mu)
  denom <- e + m * psi
  -0.5 * (m * log(2 * pi) + log(denom) + (m - 1) * log(e)) -
    (1 / (2 * e)) * (A - (psi / denom) * B^2)
}

# Robust sandwich for an M-estimator theta = (mu, psi, e) defined by
# sum_i score_i(theta) = 0, where score_i is the gradient of .cs_loglik_i.
# Scores and the total Hessian are taken by central finite differences on the
# per-subject log-likelihood (mirrors the package's .amc_information). Produces
# O(1) influence rows for the two kappas and the matching 2x2 vcov.
.cs_mest_se <- function(mu, psi, e, m, s, q, bp_chance, method, n_full = NULL) {
  theta <- c(mu = mu, psi = psi, e = e)
  n <- length(m)
  ll <- function(th) .cs_loglik_i(th[1], th[2], th[3], m, s, q)
  h <- 1e-5 * pmax(abs(theta), 1)

  G <- matrix(0, n, 3)                      # per-subject scores
  for (j in 1:3) {
    tp <- theta; tp[j] <- tp[j] + h[j]
    tm <- theta; tm[j] <- tm[j] - h[j]
    G[, j] <- (ll(tp) - ll(tm)) / (2 * h[j])
  }
  tot_score <- function(th) {
    g <- numeric(3)
    for (j in 1:3) {
      tp <- th; tp[j] <- tp[j] + h[j]
      tm <- th; tm[j] <- tm[j] - h[j]
      g[j] <- sum(ll(tp) - ll(tm)) / (2 * h[j])
    }
    g
  }
  Htot <- matrix(0, 3, 3)                    # Hessian of total loglik
  for (j in 1:3) {
    tp <- theta; tp[j] <- tp[j] + h[j]
    tm <- theta; tm[j] <- tm[j] - h[j]
    Htot[, j] <- (tot_score(tp) - tot_score(tm)) / (2 * h[j])
  }
  Hbar <- -(Htot + t(Htot)) / 2 / n          # avg per-subject observed info
  Hinv <- solve(Hbar)

  psi_theta <- G %*% Hinv                     # O(1) IF rows for (mu, psi, e)
  Jk <- .cs_kappa_jac(psi, e, bp_chance)
  psi_kappa <- psi_theta %*% t(Jk)            # n x 2
  colnames(psi_kappa) <- c("Fleiss", "Brennan-Prediger")
  vcov <- crossprod(psi_kappa) / n^2

  est <- .cs_kappas(psi, e, bp_chance)
  .cs_estimate(est, vcov, psi_kappa, method,
               list(mu = mu, psi = psi, e = e, n = n,
                    n_full = if (is.null(n_full)) n else n_full))
}

# ---- FIML -------------------------------------------------------------------

cs_fiml <- function(counts, values = NULL, r_total = NULL, method = "fiml") {
  C <- ncol(counts)
  if (is.null(values)) values <- seq_len(C)
  ss <- .cs_suff(counts, values)
  ss <- ss[ss$m > 0, , drop = FALSE]
  m <- ss$m; s <- ss$s; q <- ss$q
  bp <- .cs_bp_chance(values)

  mu0 <- sum(s) / sum(m)
  v0  <- max(sum(q) / sum(m) - mu0^2, 1e-6)
  negll <- function(par) {
    -sum(.cs_loglik_i(par[1], exp(par[3]), exp(par[2]), m, s, q))
  }
  opt <- stats::optim(c(mu0, log(0.7 * v0), log(0.3 * v0)), negll,
                      method = "BFGS",
                      control = list(reltol = 1e-13, maxit = 2000))
  mu <- opt$par[1]; e <- exp(opt$par[2]); psi <- exp(opt$par[3])
  .cs_mest_se(mu, psi, e, m, s, q, bp, method)
}

# ---- listwise (FIML on complete rows; trivial) ------------------------------

cs_listwise <- function(counts, values = NULL, r_total = NULL) {
  if (is.null(r_total)) r_total <- max(rowSums(counts))
  complete <- rowSums(counts) == r_total
  if (sum(complete) < 2L) stop("listwise needs >= 2 complete rows.")
  fit <- cs_fiml(counts[complete, , drop = FALSE], values = values,
                 method = "listwise")
  fit$moments$n_full <- nrow(counts)
  fit
}

# ---- pairwise-available covariance ------------------------------------------

# Estimate Sigma's entries from jointly-observed data, pooled by exchangeability
# (the ANOVA / pooled pairwise-available decomposition):
#   mu       = (sum s_i) / (sum m_i)                        (all tokens)
#   e        = (sum w_i) / (2 sum p_i),  w_i = m_i q_i - s_i^2   (within-pair
#              squared differences: E[(Y_j - Y_k)^2] = 2e; mu-free)
#   psi + e  = (sum q_i) / (sum m_i) - mu^2                 (all tokens: variance)
#   psi      = (psi + e) - e                                 (pairwise covariance)
# Singletons (m_i = 1) inform the variance but not the covariance, exactly as
# pairwise-available covariance does. Robust IF by the delta method on the
# per-subject sum vector z_i = (m_i, s_i, q_i, p_i, w_i).
cs_pairwise <- function(counts, values = NULL, r_total = NULL) {
  C <- ncol(counts)
  if (is.null(values)) values <- seq_len(C)
  ss <- .cs_suff(counts, values)
  ss <- ss[ss$m > 0, , drop = FALSE]
  m <- ss$m; s <- ss$s; q <- ss$q
  p <- m * (m - 1) / 2
  w <- m * q - s^2                          # sum_{j<k} (y_ij - y_ik)^2
  bp <- .cs_bp_chance(values)
  n <- length(m)
  Z <- cbind(m = m, s = s, q = q, p = p, w = w)
  zbar <- colMeans(Z)

  # kappas as a function of the mean sufficient-stat vector
  kap_of <- function(z) {
    mu  <- z["s"] / z["m"]
    e   <- z["w"] / (2 * z["p"])
    tot <- z["q"] / z["m"] - mu^2          # psi + e
    psi <- tot - e
    unname(.cs_kappas(psi, e, bp))
  }
  # central-difference Jacobian d kappa / d zbar  (2 x 5)
  hh <- 1e-6 * pmax(abs(zbar), 1)
  Jz <- matrix(0, 2, 5)
  for (j in 1:5) {
    zp <- zbar; zp[j] <- zp[j] + hh[j]
    zm <- zbar; zm[j] <- zm[j] - hh[j]
    Jz[, j] <- (kap_of(zp) - kap_of(zm)) / (2 * hh[j])
  }
  psi_kappa <- sweep(Z, 2, zbar) %*% t(Jz)   # IF_i = Jz (z_i - zbar)
  colnames(psi_kappa) <- c("Fleiss", "Brennan-Prediger")
  vcov <- crossprod(psi_kappa) / n^2

  mu  <- zbar["s"] / zbar["m"]
  e   <- zbar["w"] / (2 * zbar["p"])
  tot <- zbar["q"] / zbar["m"] - mu^2
  psi <- tot - e
  est <- .cs_kappas(psi, e, bp)
  .cs_estimate(est, vcov, psi_kappa, "pairwise",
               list(mu = unname(mu), psi = unname(psi), e = unname(e), n = n))
}

# convenience dispatcher
cs_counts <- function(counts, estimator = c("fiml", "pairwise", "listwise"),
                      values = NULL, r_total = NULL) {
  estimator <- match.arg(estimator)
  switch(estimator,
         fiml     = cs_fiml(counts, values, r_total),
         pairwise = cs_pairwise(counts, values, r_total),
         listwise = cs_listwise(counts, values, r_total))
}
