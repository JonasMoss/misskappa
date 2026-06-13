# Probe 39 — verify the C++ MLE-RE ICC estimator against lme4 (the oracle the
# paper itself used). Simulate an incomplete two-way design, fit both, compare
# variance components and the ICCs.

suppressMessages({
  library(Rcpp)
  library(lme4)
})

here <- "experiments/probes/39-ten-hove-icc"
sourceCpp(file.path(here, "icc_fiml.cpp"))

## --- estimator wrapper: optim over log-variances, mu profiled in C++ --------
fit_icc <- function(subj, rat, y, S, R, start = NULL) {
  nll <- function(par) {
    v <- exp(par)
    -icc_loglik(v[1], v[2], v[3], subj, rat, y, S, R)$loglik
  }
  if (is.null(start)) {
    v0 <- stats::var(y) / 3
    start <- log(c(v0, v0, v0))
  }
  opt <- optim(start, nll, method = "BFGS",
               control = list(reltol = 1e-12, maxit = 500))
  v <- exp(opt$par)
  names(v) <- c("s2s", "s2r", "s2sr")
  list(varcomp = v,
       mu      = icc_loglik(v[1], v[2], v[3], subj, rat, y, S, R)$mu,
       loglik  = -opt$value,
       icc_A1  = v[["s2s"]] / sum(v),
       icc_C1  = v[["s2s"]] / (v[["s2s"]] + v[["s2sr"]]))
}

## --- simulate a planned-incomplete two-way design ---------------------------
set.seed(2025)
S <- 60; R <- 12; per_subj <- 6          # each subject rated by 6 of 12 raters
mu <- 3; s2s <- 2; s2r <- 1; s2sr <- 1.5 # truth (verification is vs lme4, not truth)

a <- rnorm(S, 0, sqrt(s2s))
b <- rnorm(R, 0, sqrt(s2r))
rows <- do.call(rbind, lapply(seq_len(S), function(s) {
  rs <- sort(sample.int(R, per_subj))
  data.frame(subj = s, rat = rs,
             y = mu + a[s] + b[rs] + rnorm(per_subj, 0, sqrt(s2sr)))
}))
cat(sprintf("design: S=%d R=%d  n=%d cells  (%.0f%% of cells missing)\n",
            S, R, nrow(rows), 100 * (1 - nrow(rows) / (S * R))))

## --- fit both ---------------------------------------------------------------
mine <- fit_icc(rows$subj, rows$rat, rows$y, S, R)

m  <- lmer(y ~ 1 + (1 | subj) + (1 | rat), data = rows, REML = FALSE)
vc <- as.data.frame(VarCorr(m))
lme <- c(s2s  = vc$vcov[vc$grp == "subj"],
         s2r  = vc$vcov[vc$grp == "rat"],
         s2sr = vc$vcov[vc$grp == "Residual"])
lme_A1 <- lme[["s2s"]] / sum(lme)
lme_C1 <- lme[["s2s"]] / (lme[["s2s"]] + lme[["s2sr"]])

## --- report -----------------------------------------------------------------
cmp <- data.frame(
  quantity = c("s2s", "s2r", "s2sr", "mu", "logLik", "ICC(A,1)", "ICC(C,1)"),
  cpp  = c(mine$varcomp, mine$mu, mine$loglik, mine$icc_A1, mine$icc_C1),
  lme4 = c(lme, fixef(m)[[1]], as.numeric(logLik(m)), lme_A1, lme_C1)
)
cmp$abs_diff <- abs(cmp$cpp - cmp$lme4)
cmp$rel_diff <- cmp$abs_diff / pmax(abs(cmp$lme4), 1e-8)
print(cmp, digits = 8, row.names = FALSE)

stopifnot(max(cmp$rel_diff[cmp$quantity %in%
        c("s2s","s2r","s2sr","ICC(A,1)","ICC(C,1)")]) < 1e-4)
cat("\nPASS: C++ MLE-RE matches lme4 to <1e-4 relative on all variance components and ICCs.\n")
