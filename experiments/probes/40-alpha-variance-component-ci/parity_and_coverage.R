# Probe 40 — variance-component ("ten Hove") alpha: missing-data oracle parity
# and CI coverage.
#
# Two questions left open after commit ff96f5c "Add variance-component normal
# FIML":
#   (1) PARITY. The C++ path fits the compound-symmetry (CS) covariance by a
#       saturated EM with a CS projection in the M-step. The probe-39 oracle
#       (lme4 MLE-RE) was only ever checked on a *different* implementation
#       (Woodbury + optim) and only the point estimate. Does the package's
#       EM-projection estimate match lme4 on MISSING data?
#   (2) COVERAGE. No CI calibration of the package's analytic sandwich exists
#       for alpha at all. How does the variance-component alpha CI compare to
#       saturated nt_fiml and pairwise-available, under complete/MCAR/MAR, with
#       a true-CS DGP (model true) and a congeneric DGP (CS misspecified)?
#
# Note on the estimand: alpha depends on Sigma only through tr(Sigma) and
# sum(Sigma), and the CS projection PRESERVES both (algebraic fact, see
# NOTES.md). So on the *moments* the CS route and saturated route target the
# same alpha even when CS is false. The open question is whether the constrained
# EM fixed point keeps that invariance under MAR, and how the SEs calibrate.

suppressMessages({
  library(misskappa)
  library(MASS)
  library(lme4)
})

set.seed(40)

# ---- DGP covariance builders ------------------------------------------------

cs_sigma <- function(p, psi = 1.0, e = 1.0) psi * matrix(1, p, p) + e * diag(p)

congeneric_sigma <- function(p) {
  lambda <- seq(0.6, 1.4, length.out = p)   # unequal loadings -> not CS
  theta  <- seq(0.5, 1.3, length.out = p)   # unequal uniquenesses
  tcrossprod(lambda) + diag(theta)
}

pop_alpha <- function(S) {
  p <- ncol(S)
  (p / (p - 1)) * (1 - sum(diag(S)) / sum(S))
}

# ---- missingness mechanisms (return X with NAs) -----------------------------

amputate <- function(X, mech, rate = 0.30) {
  n <- nrow(X); p <- ncol(X)
  if (mech == "complete") return(X)
  if (mech == "mcar") {
    X[matrix(runif(n * p) < rate, n, p)] <- NA
    return(X)
  }
  if (mech == "mar") {
    # observed-anchor MAR: item 1 always observed; items 2..p go missing with
    # probability increasing in the (observed) value of item 1.
    a <- scale(X[, 1])
    pr <- plogis(-0.4 + 0.9 * a)            # ~30-40% on the missing-eligible items
    for (j in 2:p) X[runif(n) < pr, j] <- NA
    return(X)
  }
  stop("unknown mech")
}

# ---- lme4 oracle for the CS-constrained alpha (point estimate) --------------

lme4_alpha <- function(X) {
  n <- nrow(X); p <- ncol(X)
  long <- data.frame(
    value   = as.vector(X),
    subject = factor(rep(seq_len(n), times = p)),
    item    = factor(rep(seq_len(p), each = n))
  )
  long <- long[is.finite(long$value), ]
  # fixed item means, random subject intercept, homoscedastic residual =
  # Sigma = psi J + e I with a free mean vector. FIML (REML = FALSE).
  fit <- lmer(value ~ 0 + item + (1 | subject), data = long, REML = FALSE,
              control = lmerControl(calc.derivs = FALSE))
  vc  <- as.data.frame(VarCorr(fit))
  psi <- vc$vcov[vc$grp == "subject"]
  e   <- vc$vcov[vc$grp == "Residual"]
  (p / (p - 1)) * (1 - (e + psi) / (e + p * psi))
}

# =============================================================================
# PART 1 — missing-data parity: package variance_components vs lme4
# =============================================================================

cat("== PART 1: package variance_components alpha vs lme4 (missing data) ==\n")
p <- 6
parity_rows <- list()
for (dgp in c("cs", "congeneric")) {
  S0 <- if (dgp == "cs") cs_sigma(p) else congeneric_sigma(p)
  for (mech in c("mcar", "mar")) {
    X <- mvrnorm(600, rep(0, p), S0)
    Xm <- amputate(X, mech, rate = 0.30)
    pkg <- alpha(Xm, estimator = "nt_fiml",
                 em_options = list(covariance = "variance_components"))
    a_pkg  <- as.numeric(coef(pkg))
    a_lme4 <- lme4_alpha(Xm)
    parity_rows[[length(parity_rows) + 1]] <- data.frame(
      dgp = dgp, mech = mech, pkg = a_pkg, lme4 = a_lme4,
      abs_diff = abs(a_pkg - a_lme4))
  }
}
parity <- do.call(rbind, parity_rows)
print(parity, row.names = FALSE, digits = 8)
cat(sprintf("\nworst |pkg - lme4| = %.2e\n\n", max(parity$abs_diff)))

# =============================================================================
# PART 2 — CI coverage simulation
# =============================================================================

cat("== PART 2: alpha CI coverage (B reps) ==\n")
B <- 800
n <- 250
p <- 6

fit_one <- function(Xm, estimator, covariance = NULL) {
  res <- tryCatch({
    if (estimator == "pairwise") {
      f <- alpha(Xm, estimator = "pairwise")
    } else {
      f <- alpha(Xm, estimator = "nt_fiml",
                 em_options = if (is.null(covariance)) list()
                              else list(covariance = covariance))
    }
    est <- as.numeric(coef(f))
    se  <- sqrt(diag(vcov(f)))
    ci_w <- confint(f, transform = "none")
    ci_f <- confint(f, transform = "fisher")
    c(est = est, se = se,
      lo_w = ci_w[1], hi_w = ci_w[2], lo_f = ci_f[1], hi_f = ci_f[2])
  }, error = function(e) rep(NA_real_, 6))
  res
}

methods <- list(
  pairwise  = list(estimator = "pairwise"),
  saturated = list(estimator = "nt_fiml", covariance = "saturated"),
  varcomp   = list(estimator = "nt_fiml", covariance = "variance_components")
)

scenarios <- expand.grid(
  dgp  = c("cs", "congeneric"),
  mech = c("complete", "mcar", "mar"),
  stringsAsFactors = FALSE
)

summ <- list()
for (s in seq_len(nrow(scenarios))) {
  dgp <- scenarios$dgp[s]; mech <- scenarios$mech[s]
  S0 <- if (dgp == "cs") cs_sigma(p) else congeneric_sigma(p)
  truth <- pop_alpha(S0)

  acc <- lapply(names(methods), function(m) {
    matrix(NA_real_, B, 6,
           dimnames = list(NULL, c("est","se","lo_w","hi_w","lo_f","hi_f")))
  })
  names(acc) <- names(methods)

  for (b in seq_len(B)) {
    X  <- mvrnorm(n, rep(0, p), S0)
    Xm <- amputate(X, mech, rate = 0.30)
    for (m in names(methods)) {
      acc[[m]][b, ] <- do.call(fit_one, c(list(Xm = Xm), methods[[m]]))
    }
  }

  for (m in names(methods)) {
    A <- acc[[m]]
    ok <- is.finite(A[, "est"])
    est <- A[ok, "est"]; se <- A[ok, "se"]
    cov_w <- mean(A[ok, "lo_w"] <= truth & truth <= A[ok, "hi_w"])
    cov_f <- mean(A[ok, "lo_f"] <= truth & truth <= A[ok, "hi_f"])
    summ[[length(summ) + 1]] <- data.frame(
      dgp = dgp, mech = mech, method = m, truth = truth,
      conv = mean(ok),
      bias = mean(est) - truth,
      sd_emp = sd(est),
      mean_se = mean(se, na.rm = TRUE),
      se_ratio = mean(se, na.rm = TRUE) / sd(est),
      cov_wald = cov_w,
      cov_fisher = cov_f
    )
  }
}

out <- do.call(rbind, summ)
out_print <- out
num <- sapply(out_print, is.numeric)
out_print[num] <- lapply(out_print[num], function(z) round(z, 4))
print(out_print, row.names = FALSE)

saveRDS(out, "experiments/probes/40-alpha-variance-component-ci/coverage_summary.rds")
cat("\ndone.\n")
