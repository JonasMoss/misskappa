# Probe 39 — delta-method SEs for the MLE-RE ICC estimator.
#
# ICCs are smooth functions g(theta) of the variance components
# theta = (s2s, s2r, s2sr).  Delta method: SE = sqrt(grad' V_theta grad), with
#   V_theta = [ -d2/dtheta2  profile-loglik ]^{-1}   (mu already profiled out,
#             so estimation of the mean is correctly accounted for).
# Analytic ICC gradients; numeric Hessian (numDeriv) of the C++ profile loglik.

suppressMessages({ library(Rcpp); library(numDeriv) })
here <- "experiments/probes/39-ten-hove-icc"
sourceCpp(file.path(here, "icc_fiml.cpp"))

## profile negative log-lik on the natural (variance) scale --------------------
nll_nat <- function(theta, subj, rat, y, S, R) {
  if (any(theta <= 0)) return(1e10)
  -icc_loglik(theta[1], theta[2], theta[3], subj, rat, y, S, R)$loglik
}

## point estimate (log-scale BFGS) --------------------------------------------
fit_pt <- function(subj, rat, y, S, R) {
  v0 <- stats::var(y) / 3
  p  <- optim(log(rep(v0, 3)),
              function(p) nll_nat(exp(p), subj, rat, y, S, R),
              method = "BFGS", control = list(reltol = 1e-12, maxit = 500))$par
  v <- exp(p); names(v) <- c("s2s", "s2r", "s2sr")
  iccs <- c(A1 = v[["s2s"]] / sum(v),
            C1 = v[["s2s"]] / (v[["s2s"]] + v[["s2sr"]]))
  list(v = v, A1 = iccs[["A1"]], C1 = iccs[["C1"]],
       mu = icc_loglik(v[1], v[2], v[3], subj, rat, y, S, R)$mu)
}

## point estimate + delta-method SEs ------------------------------------------
fit_se <- function(subj, rat, y, S, R) {
  f <- fit_pt(subj, rat, y, S, R)
  s <- f$v[["s2s"]]; r <- f$v[["s2r"]]; e <- f$v[["s2sr"]]
  Tt <- s + r + e; D <- s + e

  H   <- numDeriv::hessian(nll_nat, f$v, subj = subj, rat = rat, y = y, S = S, R = R)
  Vth <- tryCatch(solve(H), error = function(...) matrix(NA, 3, 3))

  gA <- c((r + e) / Tt^2, -s / Tt^2, -s / Tt^2)   # d ICC(A,1) / d theta
  gC <- c(e / D^2,          0,        -s / D^2)    # d ICC(C,1) / d theta
  seA <- sqrt(as.numeric(t(gA) %*% Vth %*% gA))
  seC <- sqrt(as.numeric(t(gC) %*% Vth %*% gC))

  ci_logit <- function(g, se) {                    # [0,1]-respecting interval
    h <- qlogis(g); seh <- se / (g * (1 - g))
    plogis(h + c(-1, 1) * 1.96 * seh)
  }
  list(v = f$v, se_v = sqrt(diag(Vth)), mu = f$mu,
       A1 = f$A1, seA = seA, ciA = f$A1 + c(-1, 1) * 1.96 * seA, ciA_logit = ci_logit(f$A1, seA),
       C1 = f$C1, seC = seC, ciC = f$C1 + c(-1, 1) * 1.96 * seC, ciC_logit = ci_logit(f$C1, seC))
}

## parametric bootstrap ground-truth SD (honest: refit each replicate) --------
boot_sd <- function(d, S, R, f, B = 500, seed = 99) {
  set.seed(seed)
  v <- f$v; est <- matrix(NA_real_, B, 2)
  for (b in seq_len(B)) {
    a  <- rnorm(S, 0, sqrt(v[["s2s"]])); bb <- rnorm(R, 0, sqrt(v[["s2r"]]))
    yb <- f$mu + a[d$subj] + bb[d$rat] + rnorm(nrow(d), 0, sqrt(v[["s2sr"]]))
    fb <- tryCatch(fit_pt(d$subj, d$rat, yb, S, R), error = function(...) NULL)
    if (!is.null(fb)) est[b, ] <- c(fb$A1, fb$C1)
  }
  c(sdA = sd(est[, 1], na.rm = TRUE), sdC = sd(est[, 2], na.rm = TRUE))
}

## designs: interior (easy) vs boundary (small S, heavy missing) ---------------
sim <- function(seed, S, R, per, s2s, s2r, s2sr, mu = 3) {
  set.seed(seed)
  a <- rnorm(S, 0, sqrt(s2s)); b <- rnorm(R, 0, sqrt(s2r))
  do.call(rbind, lapply(seq_len(S), function(s) {
    rs <- sort(sample.int(R, per))
    data.frame(subj = s, rat = rs, y = mu + a[s] + b[rs] + rnorm(per, 0, sqrt(s2sr)))
  }))
}

designs <- list(
  interior = list(S = 60, R = 12, per = 6, s2s = 2,   s2r = 1,   s2sr = 1.5),
  boundary = list(S = 30, R = 12, per = 2, s2s = 4,   s2r = 0.1, s2sr = 0.5)  # high ICC, tiny rater var
)

cat("Delta-method SE vs honest parametric-bootstrap SD (B=500):\n\n")
for (nm in names(designs)) {
  g <- designs[[nm]]
  d <- sim(2025, g$S, g$R, g$per, g$s2s, g$s2r, g$s2sr)
  f <- fit_se(d$subj, d$rat, d$y, g$S, g$R)
  bs <- boot_sd(d, g$S, g$R, f)
  cat(sprintf("[%s]  S=%d R=%d  %.0f%% missing  (ICC(A,1)=%.3f  ICC(C,1)=%.3f)\n",
              nm, g$S, g$R, 100*(1 - g$per/g$R), f$A1, f$C1))
  tab <- data.frame(
    ICC   = c("A,1", "C,1"),
    est   = c(f$A1, f$C1),
    se_delta = c(f$seA, f$seC),
    sd_boot  = c(bs[["sdA"]], bs[["sdC"]]),
    ratio    = c(f$seA / bs[["sdA"]], f$seC / bs[["sdC"]]),
    CI_Wald  = c(sprintf("[%.3f, %.3f]", f$ciA[1], f$ciA[2]),
                 sprintf("[%.3f, %.3f]", f$ciC[1], f$ciC[2])),
    CI_logit = c(sprintf("[%.3f, %.3f]", f$ciA_logit[1], f$ciA_logit[2]),
                 sprintf("[%.3f, %.3f]", f$ciC_logit[1], f$ciC_logit[2]))
  )
  print(tab, row.names = FALSE, digits = 4)
  cat("\n")
}
