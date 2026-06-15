# Validate the analytic robust IF for the free-mean CS pairwise estimator
# against the subject-level bootstrap, and check MCAR coverage. Run:
#   Rscript validate_if.R
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
if (length(here) == 0 || !nzchar(here)) here <- "."
source(file.path(here, "free_cs_pairwise.R"))
source(file.path(here, "analytic_if.R"))

set.seed(20260615)

# CS DGP with free rater means: Y_ij = mu_j + b_i + eps_ij
sim_cs <- function(n, R, mu, psi, e, p_obs = 1) {
  b <- rnorm(n, 0, sqrt(psi))
  eps <- matrix(rnorm(n * R, 0, sqrt(e)), n, R)
  Y <- matrix(mu, n, R, byrow = TRUE) + b + eps
  if (p_obs < 1) {
    drop <- matrix(runif(n * R) > p_obs, n, R)   # MCAR
    Y[drop] <- NA
    Y[rowSums(!is.na(Y)) < 2, ] <- NA            # keep >=2 obs/subject
    Y <- Y[rowSums(!is.na(Y)) >= 2, , drop = FALSE]
  }
  Y
}

cat("=== analytic SE vs bootstrap SE (single large dataset per design) ===\n")
designs <- list(
  list(tag = "complete  R=4", n = 4000, R = 4, mu = c(0, .3, .6, 1), p = 1),
  list(tag = "MCAR .7   R=4", n = 6000, R = 4, mu = c(0, .3, .6, 1), p = .7),
  list(tag = "MCAR .5   R=6", n = 8000, R = 6, mu = seq(0, 1, length = 6), p = .5)
)
for (d in designs) {
  X <- sim_cs(d$n, d$R, d$mu, psi = 1, e = .5, p_obs = d$p)
  an <- cs_free_pairwise_if(X)
  bo <- cs_free_pairwise_boot(X, B = 400)
  cat(sprintf("\n%s  (n_used=%d)\n", d$tag, nrow(X)))
  print(round(rbind(estimate = an$estimates, se_analytic = an$se,
                    se_boot = bo$se, ratio = an$se / bo$se), 4))
}

cat("\n=== MCAR coverage (B reps, nominal 95%) ===\n")
cover_check <- function(reps, n, R, mu, psi, e, p) {
  truth <- {                       # population coefficients at the CS truth
    t3 <- sum((mu - mean(mu))^2); Rm1 <- R - 1
    conger <- psi / (psi + e + t3 / Rm1)
    fleiss <- (R * Rm1 * psi - t3) / (Rm1 * (R * (psi + e) + t3))
    c(Conger = conger, Fleiss = fleiss)
  }
  hit <- matrix(0, reps, 2); est <- matrix(NA_real_, reps, 2)
  for (r in seq_len(reps)) {
    X <- sim_cs(n, R, mu, psi, e, p_obs = p)
    f <- cs_free_pairwise_if(X)
    z <- qnorm(.975)
    for (k in 1:2) {
      lo <- f$estimates[k] - z * f$se[k]; hi <- f$estimates[k] + z * f$se[k]
      hit[r, k] <- truth[k] >= lo && truth[k] <= hi
      est[r, k] <- f$estimates[k]
    }
  }
  cat(sprintf("truth: Conger=%.3f Fleiss=%.3f\n", truth[1], truth[2]))
  cat(sprintf("coverage: Conger=%.3f Fleiss=%.3f | mean est: %.3f %.3f | emp SD: %.4f %.4f\n",
              mean(hit[, 1]), mean(hit[, 2]), mean(est[, 1]), mean(est[, 2]),
              sd(est[, 1]), sd(est[, 2])))
}
cover_check(reps = 600, n = 200, R = 4, mu = c(0, .3, .6, 1), psi = 1, e = .5, p = .7)
cover_check(reps = 600, n = 120, R = 5, mu = seq(0, 1, length = 5), psi = 1, e = .6, p = 1)
