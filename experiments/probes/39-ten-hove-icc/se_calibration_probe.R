# Probe 39 — settle the SE question. Under non-normal (finite-kurtosis) data with
# CORRECT second moments, is the model-based delta SE for ICC(A,1) calibrated, or
# do we need the two-way (subject x rater) cluster-robust sandwich?
# Calibration = mean(SE) / sd(estimate) across reps; also coverage.

suppressMessages({ library(lme4); library(merDeriv) })
set.seed(11)
S <- 80; R <- 6; per <- 3; mu <- 3; s2s <- 2; s2r <- 0.8; s2sr <- 1
b <- rnorm(R, 0, sqrt(s2r))                          # raters held FIXED
icc_true <- s2s / (s2s + var(b)*(R-1)/R + s2sr)      # estimand for THESE raters

# variance-preserving leptokurtic mixture (finite 4th moment, ~10x normal kurtosis)
rcontam <- function(n, v, eps = 0.1, k = 5) {
  cc <- v / ((1-eps) + eps*k^2)
  ifelse(runif(n) < eps, rnorm(n, 0, sqrt(cc*k^2)), rnorm(n, 0, sqrt(cc)))
}
gen <- function(dist) {
  rf <- if (dist == "normal") function(n, v) rnorm(n, 0, sqrt(v)) else rcontam
  a <- rf(S, s2s)
  do.call(rbind, lapply(seq_len(S), function(s) {
    rs <- sort(sample.int(R, per))
    data.frame(subj = s, rat = rs, y = mu + a[s] + b[rs] + rf(per, s2sr))
  }))
}
ig <- function(th) { s<-th[1]; r<-th[2]; e<-th[3]; Tt<-sum(th); c((r+e)/Tt^2, -s/Tt^2, -s/Tt^2) }

fit_one <- function(d) {
  m <- suppressWarnings(lmer(y ~ 1 + (1|subj) + (1|rat), d, REML = FALSE,
                             control = lmerControl(calc.derivs = FALSE)))
  vc <- as.data.frame(VarCorr(m))
  th <- c(vc$vcov[vc$grp=="subj"], vc$vcov[vc$grp=="rat"], vc$vcov[vc$grp=="Residual"])
  A  <- vcov(m, full = TRUE, ranpar = "var", information = "expected")  # 4x4 inv-info (bread)
  idx <- c(grep("subj", colnames(A)), grep("rat", colnames(A)), grep("residual", colnames(A)))
  g4  <- numeric(ncol(A)); g4[idx] <- ig(th)
  icc <- th[1]/sum(th)
  se_mod <- sqrt(as.numeric(t(g4) %*% A %*% g4))
  # two-way cluster-robust meat: M_subj + M_rat - M_obs  (CGM 2011)
  ef <- estfun(m, level = 1, ranpar = "var")          # n x 4 obs-level scores
  M  <- crossprod(rowsum(ef, d$subj)) + crossprod(rowsum(ef, d$rat)) - crossprod(ef)
  Vr <- A %*% M %*% A
  se_2w <- sqrt(max(as.numeric(t(g4) %*% Vr %*% g4), 0))
  c(icc = icc, se_mod = se_mod, se_2w = se_2w)
}

run <- function(dist, K = 200) {
  out <- t(vapply(seq_len(K), function(k) {
    set.seed(1000*(dist=="contam") + k)
    tryCatch(fit_one(gen(dist)), error = function(e) c(NA,NA,NA))
  }, numeric(3)))
  out <- out[stats::complete.cases(out), , drop = FALSE]
  sdic <- sd(out[,"icc"])
  cov <- function(se) mean(abs(out[,"icc"] - icc_true) <= 1.96*out[,se])
  data.frame(
    dist = dist, reps = nrow(out), icc_true = round(icc_true,3),
    sd_emp = round(sdic,4),
    se_mod = round(mean(out[,"se_mod"]),4), cal_mod = round(mean(out[,"se_mod"])/sdic,2),
    cov_mod = round(cov("se_mod"),3),
    se_2w  = round(mean(out[,"se_2w"]),4),  cal_2w  = round(mean(out[,"se_2w"])/sdic,2),
    cov_2w  = round(cov("se_2w"),3))
}

res <- rbind(run("normal"), run("contam"))
cat("\nCalibration of ICC(A,1) SEs  [cal = mean(SE)/sd(est); want ~1.0, cov ~0.95]\n\n")
print(res, row.names = FALSE)
