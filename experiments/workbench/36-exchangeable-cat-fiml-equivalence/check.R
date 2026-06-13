# Question: if we impose rater exchangeability on the saturated cat-FIML
# (joint multinomial over {1..C}^R with FIML over observed rater subsets),
# do we get the count-FIML of kappa_counts(estimator = "cat_fiml"), or a
# different estimator?
#
# Claim: exactly the count-FIML. Exchangeability makes the joint uniform on
# each permutation orbit, so the free parameter is theta on the composition
# simplex; the marginal likelihood of any specific observed rater subset's
# ordered ratings is proportional (in theta) to the multivariate-
# hypergeometric counts likelihood the C++ backend maximises.
#
# Check 1: raw rater-identified exchangeable loglik(theta) minus counts
#          loglik(theta) is constant in theta (additive constant only).
# Check 2: direct ML on the raw rater-identified data reproduces
#          kappa_counts(estimator = "cat_fiml") estimates and vcov.

library(misskappa)
set.seed(20260612)

R_tot <- 4L
C <- 3L
n <- 300L

# --- compositions of R_tot into C parts --------------------------------------
enumerate_compositions <- function(n_total, c_parts) {
  if (c_parts == 1L) return(matrix(n_total, nrow = 1L))
  out <- NULL
  for (i in 0:n_total) {
    rest <- enumerate_compositions(n_total - i, c_parts - 1L)
    out <- rbind(out, cbind(i, rest))
  }
  unname(out)
}
comps <- enumerate_compositions(R_tot, C)  # J x C
J <- nrow(comps)
log_multinom <- function(z) lgamma(sum(z) + 1) - sum(lgamma(z + 1))

# --- simulate exchangeable data with rater-specific missingness --------------
theta0 <- as.numeric(rdirichlet <- {g <- rgamma(J, 1); g / sum(g)})
Z <- comps[sample.int(J, n, replace = TRUE, prob = theta0), , drop = FALSE]
X <- t(apply(Z, 1L, function(z) sample(rep.int(1:C, z))))   # n x R, exchangeable
miss_prob <- c(0.0, 0.15, 0.35, 0.55)                       # rater-specific
for (j in 1:R_tot) X[runif(n) < miss_prob[j], j] <- NA

counts <- t(apply(X, 1L, function(r) tabulate(r[!is.na(r)], nbins = C)))

# --- the two log-likelihoods as functions of theta ----------------------------
# Raw rater-identified, exchangeable joint: P(specific ordered observed tuple)
#   = sum_z theta_z * multinom(R - m; z - c) / multinom(R; z).
loglik_raw <- function(theta) {
  s <- 0
  for (i in 1:n) {
    c_i <- counts[i, ]
    m <- sum(c_i)
    terms <- vapply(1:J, function(j) {
      z <- comps[j, ]
      if (any(z < c_i)) return(0)
      theta[j] * exp(log_multinom(z - c_i) - log_multinom(z))
    }, 0)
    s <- s + log(sum(terms))
  }
  s
}

# Counts likelihood (what the C++ backend maximises):
#   P(c | m) = sum_z theta_z * prod_k choose(z_k, c_k) / choose(R, m).
loglik_counts <- function(theta) {
  s <- 0
  for (i in 1:n) {
    c_i <- counts[i, ]
    m <- sum(c_i)
    terms <- vapply(1:J, function(j) {
      z <- comps[j, ]
      if (any(z < c_i)) return(0)
      theta[j] * prod(choose(z, c_i)) / choose(R_tot, m)
    }, 0)
    s <- s + log(sum(terms))
  }
  s
}

# --- Check 1: difference constant in theta ------------------------------------
diffs <- replicate(25, {
  g <- rgamma(J, 1); th <- g / sum(g)
  loglik_raw(th) - loglik_counts(th)
})
cat(sprintf("Check 1: loglik_raw - loglik_counts over 25 random theta:\n  range = [%.12f, %.12f], spread = %.3e\n",
            min(diffs), max(diffs), diff(range(diffs))))
stopifnot(diff(range(diffs)) < 1e-8)

# --- Check 2: direct raw-data ML vs kappa_counts(cat_fiml) --------------------
softmax <- function(eta) {p <- exp(c(eta, 0) - max(c(eta, 0))); p / sum(p)}
fit <- optim(rep(0, J - 1L), function(eta) -loglik_raw(softmax(eta)),
             method = "BFGS", control = list(maxit = 1000, reltol = 1e-14))
theta_hat <- softmax(fit$par)

# Map theta -> (Fleiss, BP), nominal loss L = 1 - I (mirrors map_to_kappa).
L <- matrix(1, C, C) - diag(C)
d_z <- apply(comps, 1L, function(z) {
  (drop(z %*% L %*% z) - sum(diag(L) * z)) / (R_tot * (R_tot - 1))
})
pd <- sum(d_z * theta_hat)
p_hat <- colSums(comps * theta_hat) / R_tot
pe <- drop(p_hat %*% L %*% p_hat)
fleiss_raw <- 1 - pd / pe
bp_raw <- 1 - pd / (sum(L) / C^2)

pkg <- kappa_counts(counts, estimator = "cat_fiml", r_total = R_tot,
                    em_options = list(tol = 1e-12, max_iter = 100000))
cat(sprintf("\nCheck 2: estimates\n  raw-ML  Fleiss = %.10f, BP = %.10f\n  package Fleiss = %.10f, BP = %.10f\n",
            fleiss_raw, bp_raw, pkg$estimates[1], pkg$estimates[2]))
stopifnot(abs(fleiss_raw - pkg$estimates[1]) < 1e-6,
          abs(bp_raw - pkg$estimates[2]) < 1e-6)

# SE: observed-information sandwich-free ML vcov on the reference-removed
# parametrisation, delta method to (Fleiss, BP); compare with the package's
# Louis-information vcov.
# Prune near-zero compositions (the package does the same before its Louis
# information pass); FD steps would otherwise leave the simplex.
active <- which(theta_hat > 1e-6)
ref <- active[length(active)]
free_idx <- active[-length(active)]
loglik_free <- function(th_free) {
  th <- numeric(J)
  th[free_idx] <- th_free
  th[ref] <- 1 - sum(th_free)
  loglik_raw(th)
}
H <- optimHess(theta_hat[free_idx], loglik_free,
               control = list(ndeps = rep(1e-6, length(free_idx))))
V_theta_free <- solve(-H)
grad_kappa <- function(th) {  # 2 x J jacobian of (Fleiss, BP) wrt full theta
  pd_ <- sum(d_z * th)
  p_ <- colSums(comps * th) / R_tot
  pe_ <- drop(p_ %*% L %*% p_)
  g_pe <- 2 * drop(comps %*% (L %*% p_)) / R_tot
  rbind(-(d_z * pe_ - pd_ * g_pe) / pe_^2,
        -d_z / (sum(L) / C^2))
}
Jk <- grad_kappa(theta_hat)
Jk_free <- Jk[, free_idx, drop = FALSE] -
  Jk[, ref] %o% rep(1, length(free_idx))  # theta_ref = 1 - sum(free)
V_kappa_raw <- Jk_free %*% V_theta_free %*% t(Jk_free)
cat(sprintf("\nCheck 2: SEs\n  raw-ML  se(Fleiss) = %.8f, se(BP) = %.8f\n  package se(Fleiss) = %.8f, se(BP) = %.8f\n",
            sqrt(V_kappa_raw[1, 1]), sqrt(V_kappa_raw[2, 2]),
            sqrt(pkg$vcov[1, 1]), sqrt(pkg$vcov[2, 2])))
stopifnot(abs(sqrt(V_kappa_raw[1, 1]) - sqrt(pkg$vcov[1, 1])) < 1e-4,
          abs(sqrt(V_kappa_raw[2, 2]) - sqrt(pkg$vcov[2, 2])) < 1e-4)

cat("\nAll checks passed: exchangeability-constrained cat-FIML == count-FIML (estimates + SEs).\n")
