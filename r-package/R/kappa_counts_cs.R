# Direct-on-counts exchangeable compound-symmetry (CS) quadratic kappa.
#
# Fits the exchangeable common-mean CS normal model Sigma = psi J + e I straight
# from the counts, with no realizing subjects-by-raters table, no r_total, and no
# EM. Backs agreement_counts(estimator = "nt_fiml_cs") for weight = "quadratic" (the
# moment / available-case slot is filled by the Fleiss-Cuzick estimator).
# A subject contributes only its observed sub-vector of
# size m_i; under exchangeability that marginal is again CS-normal with the same
# (mu, psi, e), so the per-subject likelihood depends on the data only through
#   m_i = sum_c n_ic,  s_i = sum_c n_ic v_c,  q_i = sum_c n_ic v_c^2.
# The kappas are r_total-invariant smooth maps of (psi, e):
#   Fleiss = psi / (psi + e),  Brennan-Prediger = 1 - 2 e / bp_chance,
# matching .kqf_grad() at the equal-mean CS point. All routes return O(1)
# per-subject influence rows obeying vcov = crossprod(psi) / n^2, so the SEs are
# robust (empirical sandwich) and research_joint_vcov()-compatible. Derivation and
# validation: dev/notes/counts-exchangeable-cs-direct.md, probe 41.

# ---- sufficient statistics and coefficient maps -----------------------------

.cs_counts_suff <- function(counts, values) {
  counts <- as.matrix(counts)
  storage.mode(counts) <- "double"
  list(m = rowSums(counts),
       s = as.numeric(counts %*% values),
       q = as.numeric(counts %*% (values^2)))
}

.cs_counts_kappas <- function(psi, e, bp_chance) {
  c(Fleiss = psi / (psi + e),
    `Brennan-Prediger` = 1 - 2 * e / bp_chance)
}

# d(Fleiss, BP) / d(mu, psi, e); the mu column is zero.
.cs_counts_kappa_jac <- function(psi, e, bp_chance) {
  tot <- psi + e
  J <- matrix(0, 2, 3, dimnames = list(c("Fleiss", "Brennan-Prediger"),
                                       c("mu", "psi", "e")))
  J["Fleiss", "psi"] <- e / tot^2
  J["Fleiss", "e"]   <- -psi / tot^2
  J["Brennan-Prediger", "e"] <- -2 / bp_chance
  J
}

# Per-subject CS log density of the observed sub-vector, from (m, s, q).
.cs_counts_loglik_i <- function(mu, psi, e, m, s, q) {
  A <- q - 2 * mu * s + m * mu^2
  B <- s - m * mu
  denom <- e + m * psi
  -0.5 * (m * log(2 * pi) + log(denom) + (m - 1) * log(e)) -
    (1 / (2 * e)) * (A - (psi / denom) * B^2)
}

# ---- robust sandwich for the M-estimator (mu, psi, e) -----------------------

# Casewise scores and the total observed information by central finite
# differences on the per-subject log-likelihood (mirrors .amc_information).
# Returns O(1) kappa influence rows and the matching 2x2 vcov.
.cs_counts_mest_se <- function(mu, psi, e, m, s, q, bp_chance) {
  theta <- c(mu, psi, e)
  n <- length(m)
  ll <- function(th) .cs_counts_loglik_i(th[1], th[2], th[3], m, s, q)
  h <- 1e-5 * pmax(abs(theta), 1)

  G <- matrix(0, n, 3)
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
  Htot <- matrix(0, 3, 3)
  for (j in 1:3) {
    tp <- theta; tp[j] <- tp[j] + h[j]
    tm <- theta; tm[j] <- tm[j] - h[j]
    Htot[, j] <- (tot_score(tp) - tot_score(tm)) / (2 * h[j])
  }
  Hbar <- -(Htot + t(Htot)) / 2 / n
  Hinv <- solve(Hbar)

  Jk <- .cs_counts_kappa_jac(psi, e, bp_chance)
  psi_kappa <- (G %*% Hinv) %*% t(Jk)
  colnames(psi_kappa) <- c("Fleiss", "Brennan-Prediger")
  list(estimates = .cs_counts_kappas(psi, e, bp_chance),
       vcov = crossprod(psi_kappa) / n^2,
       psi = psi_kappa)
}

# ---- the three estimators ---------------------------------------------------

.cs_counts_fiml <- function(counts, values, bp_chance, em_options = list()) {
  ss <- .cs_counts_suff(counts, values)
  keep <- ss$m > 0
  m <- ss$m[keep]; s <- ss$s[keep]; q <- ss$q[keep]
  if (length(m) < 2L) stop("counts normal-FIML needs >= 2 rated subjects.",
                           call. = FALSE)

  mu0 <- sum(s) / sum(m)
  v0  <- max(sum(q) / sum(m) - mu0^2, 1e-6)
  reltol <- if (!is.null(em_options$tol)) em_options$tol else 1e-13
  maxit  <- if (!is.null(em_options$max_iter)) em_options$max_iter else 2000L
  negll <- function(par) {
    -sum(.cs_counts_loglik_i(par[1], exp(par[3]), exp(par[2]), m, s, q))
  }
  opt <- stats::optim(c(mu0, log(0.7 * v0), log(0.3 * v0)), negll,
                      method = "BFGS",
                      control = list(reltol = reltol, maxit = maxit))
  mu <- opt$par[1]; e <- exp(opt$par[2]); psi <- exp(opt$par[3])
  se <- .cs_counts_mest_se(mu, psi, e, m, s, q, bp_chance)
  list(se = se, moments = list(mu = mu, psi = psi, e = e, n = length(m),
                               converged = opt$convergence == 0L))
}

# ---- entry point (agreement_counts(estimator = "nt_fiml_cs")) --------------------

.kappa_counts_cs <- function(x_mat, values, em_options = list()) {
  C <- ncol(x_mat)
  if (is.null(values)) values <- seq_len(C)
  bp_chance <- .quadratic_bp_chance(x_mat, values)
  res <- .cs_counts_fiml(x_mat, values, bp_chance, em_options)
  .new_estimate(res$se, c("Fleiss", "Brennan-Prediger"),
                method = "nt_fiml_cs", weight = "quadratic",
                moments = res$moments, bp_chance = bp_chance)
}
