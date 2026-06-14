# Probe 39 — comparison sim: ten Hove MLE-RE (structured ICC) vs misskappa
# nt_fiml (saturated Conger quadratic kappa) under incomplete raw data.
#
# Estimand: population quadratically-weighted Conger kappa for the FIXED rater
# panel (what nt_fiml targets; = ICC(A,1) under the additive Gaussian model).
# Four DGPs: model-true, then three misspecifications of MLE-RE's assumptions.
#   true        : additive Gaussian, homoscedastic, unit loadings  (MLE-RE correct)
#   hetero      : rater-specific error variances                   (breaks equal residual var)
#   congeneric  : unequal subject-effect loadings per rater        (breaks parallel/unit loadings)
#   categorical : additive Gaussian latent thresholded to C codes  (breaks Gaussian + additivity)
#
# Reports per DGP, per estimator: %bias (= inconsistency), SD & RMSE (= efficiency),
# mean SE and 95% coverage. MLE-RE SE in two flavors: model-based and two-way
# (subject x rater) cluster-robust sandwich (CGM 2011) via merDeriv scores.

suppressMessages({ library(misskappa); library(lme4); library(merDeriv) })

S <- 80; R <- 6; per <- 3; mu <- 3; s2s <- 2; s2sr <- 1; C <- 4
B <- if (exists("B_OVERRIDE")) B_OVERRIDE else 150
set.seed(101)
b      <- rnorm(R, 0, sqrt(0.8))                 # FIXED rater means (fixed-rater estimand)
vr     <- c(0.3, 0.5, 1.0, 1.5, 2.0, 3.0)        # hetero: rater error variances (mean ~1.4)
lambda <- c(0.5, 0.7, 1.0, 1.0, 1.3, 1.5)        # congeneric: rater loadings
# fixed category thresholds from the pooled latent marginal
zmarg_sd <- sqrt(s2s + s2sr)
thr <- mu + mean(b) + qnorm(seq_len(C-1)/C) * zmarg_sd

gen_complete <- function(dist, n) {
  a  <- rnorm(n, 0, sqrt(s2s))
  Bm <- matrix(b, n, R, byrow = TRUE)
  E  <- matrix(rnorm(n*R), n, R)
  if (dist == "true")        mu + outer(a, rep(1,R)) + Bm + E*sqrt(s2sr)
  else if (dist == "hetero") mu + outer(a, rep(1,R)) + Bm + sweep(E, 2, sqrt(vr), "*")
  else if (dist == "congeneric") mu + outer(a, lambda) + Bm + E*sqrt(s2sr)
  else if (dist == "categorical") {
    z <- mu + outer(a, rep(1,R)) + Bm + E*sqrt(s2sr)
    matrix(findInterval(z, thr) + 1, n, R)
  } else stop(dist)
}

conger_q <- function(m, Sig) {                   # population quadratic Conger kappa
  num <- den <- 0
  for (i in 1:(R-1)) for (j in (i+1):R) {
    md2 <- (m[i]-m[j])^2
    num <- num + Sig[i,i]+Sig[j,j]-2*Sig[i,j]+md2
    den <- den + Sig[i,i]+Sig[j,j]+md2
  }
  1 - num/den
}
truth_of <- function(dist) {
  set.seed(9999); Y <- gen_complete(dist, 2e5); conger_q(colMeans(Y), cov(Y))
}

ig <- function(th) { s<-th[1]; r<-th[2]; e<-th[3]; Tt<-sum(th); c((r+e)/Tt^2, -s/Tt^2, -s/Tt^2) }
fit_mlere <- function(d) {
  m <- suppressWarnings(lmer(y ~ 1 + (1|subj) + (1|rat), d, REML = FALSE,
                             control = lmerControl(calc.derivs = FALSE)))
  vc <- as.data.frame(VarCorr(m))
  th <- c(vc$vcov[vc$grp=="subj"], vc$vcov[vc$grp=="rat"], vc$vcov[vc$grp=="Residual"])
  A  <- vcov(m, full = TRUE, ranpar = "var", information = "expected")
  idx <- c(grep("subj",colnames(A)), grep("rat",colnames(A)), grep("residual",colnames(A)))
  g4 <- numeric(ncol(A)); g4[idx] <- ig(th)
  se_mod <- sqrt(as.numeric(t(g4) %*% A %*% g4))
  ef <- estfun(m, level = 1, ranpar = "var")
  M  <- crossprod(rowsum(ef, d$subj)) + crossprod(rowsum(ef, d$rat)) - crossprod(ef)
  se_2w <- sqrt(max(as.numeric(t(g4) %*% (A %*% M %*% A) %*% g4), 0))
  c(icc = th[1]/sum(th), se_mod = se_mod, se_2w = se_2w)
}

one_rep <- function(dist, seed) {
  set.seed(seed)
  Y <- gen_complete(dist, S)
  X <- Y; for (s in 1:S) X[s, sample(R, R-per)] <- NA
  nt <- tryCatch(kappa(X, estimator="nt_fiml", weight="quadratic"), error=function(e) NULL)
  d  <- data.frame(subj = factor(as.integer(row(X))), rat = factor(as.integer(col(X))),
                   y = as.vector(X))
  d  <- d[!is.na(d$y), ]
  re <- tryCatch(fit_mlere(d), error = function(e) c(icc=NA, se_mod=NA, se_2w=NA))
  c(est_nt = if (is.null(nt)) NA else unname(coef(nt)["Conger"]),
    se_nt  = if (is.null(nt)) NA else unname(sqrt(diag(vcov(nt))["Conger"])),
    est_re = unname(re["icc"]), se_mod = unname(re["se_mod"]), se_2w = unname(re["se_2w"]))
}

dgps <- c("true", "hetero", "congeneric", "categorical")
rows <- list()
for (dist in dgps) {
  tru <- truth_of(dist)
  M <- t(vapply(seq_len(B), function(k) one_rep(dist, 7000 + k), numeric(5)))
  M <- M[stats::complete.cases(M[, c("est_nt","est_re")]), , drop = FALSE]
  message(sprintf("  %-12s truth=%.3f  reps=%d", dist, tru, nrow(M)))
  summ <- function(est, se, tru) {
    bias <- mean(est) - tru; sdv <- sd(est)
    c(pct_bias = 100*bias/tru, sd = sdv, rmse = sqrt(mean((est-tru)^2)),
      mean_se = mean(se), cov = mean(abs(est - tru) <= 1.96*se))
  }
  nt <- summ(M[,"est_nt"], M[,"se_nt"], tru)
  rm <- summ(M[,"est_re"], M[,"se_mod"], tru)
  rr <- summ(M[,"est_re"], M[,"se_2w"], tru)
  rows[[dist]] <- data.frame(
    DGP = dist, truth = round(tru,3), reps = nrow(M),
    estimator = c("nt_fiml(sat)", "MLE-RE(mod)", "MLE-RE(2way)"),
    pct_bias = round(c(nt["pct_bias"], rm["pct_bias"], rr["pct_bias"]),1),
    SD = round(c(nt["sd"], rm["sd"], rr["sd"]),4),
    RMSE = round(c(nt["rmse"], rm["rmse"], rr["rmse"]),4),
    mean_SE = round(c(nt["mean_se"], rm["mean_se"], rr["mean_se"]),4),
    cov95 = round(c(nt["cov"], rm["cov"], rr["cov"]),3))
}
out <- do.call(rbind, rows)
cat(sprintf("\nComparison: MLE-RE (structured ICC) vs nt_fiml (saturated Conger)\nS=%d R=%d %d%% missing, B=%d, fixed raters\n\n",
            S, R, round(100*(1-per/R)), B))
print(out, row.names = FALSE)
