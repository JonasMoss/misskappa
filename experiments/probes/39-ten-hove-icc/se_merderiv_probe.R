# Probe 39 — does merDeriv give a usable robust (sandwich) SE for ICC(A,1) on a
# CROSSED lmer fit, and does it actually fix coverage under non-normality where
# the model-based SE fails?  De-risk before building the sim around it.

if (!requireNamespace("merDeriv", quietly = TRUE)) {
  install.packages("merDeriv", repos = "https://packagemanager.posit.co/cran/__linux__/noble/latest")
}
suppressMessages({ library(lme4); library(merDeriv); library(numDeriv) })

set.seed(7)
S <- 80; R <- 6; per <- 3
mu <- 3; s2s <- 2; s2r <- 0.8; s2sr <- 1
b <- rnorm(R, 0, sqrt(s2r))                 # raters held FIXED across reps
icc_true <- s2s / (s2s + var(b) * (R-1)/R + s2sr)  # placeholder; recomputed in sim

## one heavy-tailed (non-normal) data set: t_3 subject + error effects ---------
gen <- function(seed) {
  set.seed(seed)
  a <- rt(S, df = 3) * sqrt(s2s / 3)         # scaled so Var = s2s, but kurtotic
  do.call(rbind, lapply(seq_len(S), function(s) {
    rs <- sort(sample.int(R, per))
    e  <- rt(per, df = 3) * sqrt(s2sr / 3)
    data.frame(subj = s, rat = rs, y = mu + a[s] + b[rs] + e)
  }))
}
d <- gen(1)

fit_lmer <- function(d) suppressWarnings(
  lmer(y ~ 1 + (1 | subj) + (1 | rat), data = d, REML = FALSE,
       control = lmerControl(calc.derivs = FALSE)))
m <- fit_lmer(d)

## pull variance components + their vcov from merDeriv -------------------------
vc <- as.data.frame(VarCorr(m))
theta <- c(s2s = vc$vcov[vc$grp=="subj"], s2r = vc$vcov[vc$grp=="rat"],
           s2sr = vc$vcov[vc$grp=="Residual"])
icc <- function(th) th[1] / sum(th)
g <- local({ s<-theta[1]; r<-theta[2]; e<-theta[3]; Tt<-sum(theta)
             c((r+e)/Tt^2, -s/Tt^2, -s/Tt^2) })

cat("--- merDeriv objects ---\n")
Vmod <- vcov(m, full = TRUE, ranpar = "var", information = "expected")
cat("vcov(full,ranpar=var) dimnames:\n"); print(colnames(Vmod))

ef <- tryCatch(estfun(m, level = 2, ranpar = "var"), error = function(e) {
  cat("estfun level=2 failed:", conditionMessage(e), "\n")
  tryCatch(estfun(m, ranpar = "var"), error = function(e2) {
    cat("estfun default failed:", conditionMessage(e2), "\n"); NULL })
})
cat("estfun dim:", if (is.null(ef)) "NULL" else paste(dim(ef), collapse="x"), "\n")
if (!is.null(ef)) cat("estfun colnames:\n  ", paste(colnames(ef), collapse=" | "), "\n")

## locate the 3 variance-component columns by name ----------------------------
vc_idx <- function(nm) {
  hit <- function(pat) grep(pat, nm, ignore.case = TRUE)
  c(s2s = hit("subj")[1], s2r = hit("rat")[1],
    s2sr = hit("residual|sigma")[1])
}
idx <- vc_idx(colnames(Vmod))
cat("var-comp column idx in Vmod:", paste(idx, collapse=","), "\n")

se_model <- sqrt(as.numeric(t(g) %*% Vmod[idx, idx] %*% g))

se_robust <- NA
if (!is.null(ef)) {
  efi <- vc_idx(colnames(ef))
  A <- Vmod[idx, idx]                       # bread (model-based vcov of varcomps)
  meat <- crossprod(ef[, efi])              # sum of casewise score outer products
  Vrob <- A %*% meat %*% A
  se_robust <- sqrt(as.numeric(t(g) %*% Vrob %*% g))
}

## ground truth: nonparametric subject bootstrap of ICC(A,1) -------------------
B <- 300
boot <- replicate(B, {
  ss <- sample.int(S, S, replace = TRUE)
  db <- do.call(rbind, lapply(seq_along(ss), function(k) {
    r <- d[d$subj == ss[k], ]; r$subj <- k; r }))
  mb <- tryCatch(fit_lmer(db), error = function(e) NULL)
  if (is.null(mb)) return(NA)
  v <- as.data.frame(VarCorr(mb))
  th <- c(v$vcov[v$grp=="subj"], v$vcov[v$grp=="rat"], v$vcov[v$grp=="Residual"])
  th[1]/sum(th)
})
sd_boot <- sd(boot, na.rm = TRUE)

cat(sprintf("\nICC(A,1) = %.4f\n", icc(theta)))
cat(sprintf("SE model-based : %.4f\n", se_model))
cat(sprintf("SE robust(sand): %.4f\n", se_robust))
cat(sprintf("SD subj-bootstrap (truth, non-normal): %.4f\n", sd_boot))
cat(sprintf("\nmodel/boot = %.2f   robust/boot = %.2f\n",
            se_model/sd_boot, se_robust/sd_boot))
