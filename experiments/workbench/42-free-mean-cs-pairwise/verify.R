# Probe 42 — free-mean CS pairwise moment estimator ("pairwise +
# variance_components"). Is it a principled generalisation of the pairwise?
#   A NESTING.  Does it reduce to the equal-mean (one-way) pairwise when means
#               are equal (and t3_hat -> 0)?
#   B FIML.     Does it match the free-mean CS FIML (raw nt_fiml,
#               variance_components, mean=free) in the limit (complete data)?
#   C t3 FIX.   Is the t3 mean-spread bias correction necessary and correct
#               (without it, Conger is spuriously biased at small n)?
#   D CALIB.    MCAR-consistent / MAR-biased (moment estimator); FIML consistent
#               under both; bootstrap SE calibrated under MCAR.

suppressMessages(library(misskappa))
source("experiments/workbench/42-free-mean-cs-pairwise/free_cs_pairwise.R")
set.seed(42)

R <- 5L; PSI <- 0.8; E <- 1.2
MU_EQ <- rep(3, R)
MU_FR <- c(2, 2.5, 3, 3.5, 4)          # t3_pop = 2.5

pop <- function(mu) {
  t3 <- sum((mu - mean(mu))^2)
  c(Conger = PSI / (PSI + E + t3 / (R - 1)),
    Fleiss = (R * (R - 1) * PSI - t3) / ((R - 1) * (R * (PSI + E) + t3)))
}
sim_raw <- function(n, mu) {
  b <- rnorm(n, 0, sqrt(PSI))
  outer(b, rep(1, R)) + matrix(rep(mu, each = n), n, R) +
    matrix(rnorm(n * R, 0, sqrt(E)), n, R)
}
amputate <- function(Y, mech, rate = 0.3) {
  n <- nrow(Y)
  if (mech == "complete") return(Y)
  if (mech == "mcar") { Y[matrix(runif(n * R) < rate, n, R)] <- NA; return(Y) }
  if (mech == "mar") {
    pr <- plogis(-0.5 + 0.8 * scale(Y[, 1])[, 1])     # depends on observed r1
    for (j in 2:R) Y[runif(n) < pr, j] <- NA
    return(Y)
  }
}
fiml <- function(Y) {
  f <- kappa(Y, estimator = "nt_fiml", weight = "quadratic",
             em_options = list(covariance = "variance_components"))
  list(est = coef(f)[c("Conger", "Fleiss", "Brennan-Prediger")],
       se  = sqrt(diag(vcov(f)))[c("Conger", "Fleiss", "Brennan-Prediger")])
}

cat("pop: equal-mean", round(pop(MU_EQ), 4), "| free-mean", round(pop(MU_FR), 4), "\n\n")

# ============================================================ A NESTING
cat("== A. NESTING: free reduces to equal-mean pairwise (equal-mean DGP) ==\n")
for (n in c(60, 250, 2000)) {
  Y <- sim_raw(n, MU_EQ)
  fr <- cs_free_pairwise(Y, equal_mean = FALSE)
  eq <- cs_free_pairwise(Y, equal_mean = TRUE)
  cat(sprintf("n=%4d  free:  C=%.4f F=%.4f t3hat=%.4f | equal: C=%.4f F=%.4f | |free-equal|_C=%.2e\n",
              n, fr$estimates[1], fr$estimates[2], fr$t3,
              eq$estimates[1], eq$estimates[2], abs(fr$estimates[1] - eq$estimates[1])))
}

# ============================================================ B vs FIML
cat("\n== B. vs free-mean CS FIML (complete data, mean over 200 reps) ==\n")
for (lab in c("equal", "free")) {
  mu <- if (lab == "equal") MU_EQ else MU_FR
  D <- t(replicate(200, {
    Y <- sim_raw(200, mu)
    mo <- cs_free_pairwise(Y)$estimates
    fi <- fiml(Y)$est
    c(mo[1:2], fi[1:2])
  }))
  cat(sprintf("%-5s  moment C=%.4f F=%.4f | FIML C=%.4f F=%.4f | mean|mo-fiml|_C=%.4f\n",
      lab, mean(D[,1]), mean(D[,2]), mean(D[,3]), mean(D[,4]), mean(abs(D[,1]-D[,3]))))
}

# ============================================================ C t3 correction
cat("\n== C. t3 bias correction (equal-mean DGP: true t3=0, Conger should = Fleiss) ==\n")
cat("   naive uses t3=sum(mu_j-bar)^2 (no correction); corrected subtracts tr(AV).\n")
conger_from <- function(psi, e, t3) psi / (psi + e + t3 / (R - 1))
for (n in c(30, 60, 120, 500)) {
  D <- t(replicate(400, {
    Y <- sim_raw(n, MU_EQ)
    f <- cs_free_pairwise(Y)
    c(corrected = f$estimates[1],
      naive = conger_from(f$psi, f$e, f$t3_naive),
      fleiss = f$estimates[2])
  }))
  cat(sprintf("n=%4d  Conger naive=%.4f  corrected=%.4f  Fleiss=%.4f  (truth 0.4000)\n",
              n, mean(D[,1+1]), mean(D[,1]), mean(D[,3])))
}

# ============================================================ D calibration
cat("\n== D. CALIBRATION (free-mean DGP, n=200, B=300) ==\n")
truth <- fiml(sim_raw(200000, MU_FR))$est        # complete-data estimand
cat(sprintf("estimand (huge-n FIML): Conger=%.4f Fleiss=%.4f BP=%.4f\n\n",
            truth[1], truth[2], truth[3]))
B <- 300; n <- 200
rows <- list()
for (mech in c("complete", "mcar", "mar")) {
  M <- matrix(NA_real_, B, 6); Fl <- matrix(NA_real_, B, 6)   # moment / fiml: est1-3, se1-3
  for (b in seq_len(B)) {
    Y <- amputate(sim_raw(n, MU_FR), mech)
    mb <- tryCatch(cs_free_pairwise_boot(Y, B = 150), error = function(e) NULL)
    fb <- tryCatch(fiml(Y), error = function(e) NULL)
    if (!is.null(mb)) M[b, ] <- c(mb$estimates, mb$se)
    if (!is.null(fb)) Fl[b, ] <- c(fb$est, fb$se)
  }
  for (ci in 1:3) {
    nm <- c("Conger","Fleiss","BP")[ci]
    for (meth in c("moment","fiml")) {
      A <- if (meth=="moment") M else Fl
      est <- A[, ci]; se <- A[, ci+3]; ok <- is.finite(est)
      cov <- mean(abs(est[ok]-truth[ci]) <= qnorm(0.975)*se[ok])
      rows[[length(rows)+1]] <- data.frame(mech=mech, coef=nm, method=meth,
        bias=mean(est[ok])-truth[ci], sd_emp=sd(est[ok]),
        mean_se=mean(se[ok]), se_ratio=mean(se[ok])/sd(est[ok]), cov=cov)
    }
  }
}
out <- do.call(rbind, rows)
op <- out; num <- sapply(op, is.numeric); op[num] <- lapply(op[num], round, 4)
print(op[op$coef %in% c("Conger","Fleiss"), ], row.names = FALSE)
saveRDS(out, "experiments/workbench/42-free-mean-cs-pairwise/calibration.rds")
cat("\ndone.\n")
