# Probe 44 — finite-sample performance of the vector-valued quadratic kappa
# estimators (subjects x raters x features arrays, squared vector loss
# (x - y)' W (x - y)). Probe 43 locked the MATH on complete data at large n
# (block-averaging projection == Gaussian MLE; point consistency at n=8k; exact
# d=1 reduction; vcov == crossprod(psi)/n^2). This probe is the finite-sample
# PERFORMANCE story, the vector analogue of probe 40 (alpha variance-component
# CI): bias / RMSE / coverage / SE calibration under complete / MCAR / MAR, for
#
#   pairwise    -- pairwise-available covariance moment (MCAR; Conger/Fleiss),
#   nt_fiml     -- saturated normal-theory FIML (+ Brennan-Prediger),
#   nt_fiml_cs  -- block compound symmetry  Sigma = J_R (x) Psi + I_R (x) E.
#
# Headline question. Under MCAR the block-CS projection PRESERVES the contractions
# Conger/Fleiss read off:  tB = sum((J_R (x) W) . Sigma)  depends only on the sum
# of all blocks, tT = sum((I_R (x) W) . Sigma) on the sum of diagonal blocks, and
# block-averaging preserves both; tG is mean-only. So nt_fiml_cs should target the
# SAME Conger/Fleiss as nt_fiml even when block-CS is FALSE -- consistent under
# MCAR-misspecification, exactly as scalar alpha is. The genuinely open cell is
# MAR x misspecified: the constrained EM fixed point need not reproduce that
# invariance under MAR reweighting, so nt_fiml_cs may carry a small bias there
# while saturated nt_fiml stays consistent. Confirm or refute.
#
# Usage:  Rscript parity_and_coverage.R            # screen (default, minutes)
#         Rscript parity_and_coverage.R --targeted # heavier grid

suppressMessages(library(misskappa))

args <- commandArgs(trailingOnly = TRUE)
MODE <- if ("--targeted" %in% args) "targeted" else "screen"
HERE <- "experiments/probes/44-vector-quadratic-coverage"
set.seed(44)

# ---- grid by mode -----------------------------------------------------------
if (MODE == "screen") {
  DGPS <- c("blockcs", "camps")
  N_GRID <- c(120L, 300L)
  D_GRID <- 2L
  W_GRID <- "identity"
  B <- 300L
} else {
  DGPS <- c("blockcs", "camps", "mixt")
  N_GRID <- c(80L, 150L, 300L, 600L)
  D_GRID <- c(1L, 2L, 3L)
  W_GRID <- c("identity", "coupling")
  B <- 1000L
}
R <- 4L
MECHS <- c("complete", "mcar", "mar")
METHODS <- c("pairwise", "nt_fiml", "nt_fiml_cs")
COEFS3 <- c("Conger", "Fleiss", "Brennan-Prediger")
coefs_of <- function(method)
  if (method == "pairwise") c("Conger", "Fleiss") else COEFS3
cat(sprintf("probe 44 [%s]  R=%d  B=%d  dgps={%s}  n={%s}  d={%s}  W={%s}\n\n",
            MODE, R, B, paste(DGPS, collapse = ","),
            paste(N_GRID, collapse = ","), paste(D_GRID, collapse = ","),
            paste(W_GRID, collapse = ",")))

# ---- block-CS / moment building blocks (reuse probe 43 conventions) ---------
make_sigma <- function(R, Psi, E)                       # J_R (x) Psi + I_R (x) E
  kronecker(matrix(1, R, R), Psi) + kronecker(diag(R), E)
blk <- function(S, r, s, d) S[((r - 1) * d + 1):(r * d), ((s - 1) * d + 1):(s * d)]

project_block_cs <- function(S, R, d) {                 # block averaging
  Dbar <- matrix(0, d, d); Obar <- matrix(0, d, d)
  for (r in 1:R) Dbar <- Dbar + blk(S, r, r, d)
  Dbar <- Dbar / R
  for (r in 1:R) for (s in 1:R) if (r != s) Obar <- Obar + blk(S, r, s, d)
  Obar <- Obar / (R * (R - 1))
  Psi <- (Obar + t(Obar)) / 2
  E <- (Dbar - Psi); E <- (E + t(E)) / 2
  list(Psi = Psi, E = E)
}

# population Conger / Fleiss from (mu stacked Rd, Sigma Rd x Rd, W d x d)
vec_coef <- function(mu, Sigma, R, d, W) {
  Bd <- kronecker(matrix(1, R, R), W)
  Td <- kronecker(diag(R), W)
  Gd <- kronecker(diag(R) - matrix(1 / R, R, R), W)
  tB <- sum(Bd * Sigma); tT <- sum(Td * Sigma)
  tG <- as.numeric(t(mu) %*% Gd %*% mu)
  c(Conger = (tB - tT) / ((R - 1) * tT + R * tG),
    Fleiss = (tB - tT - tG) / ((R - 1) * (tT + tG)))
}

psd_floor <- function(S, eps = 1e-6) {                  # symmetric + PD repair
  S <- (S + t(S)) / 2
  ev <- eigen(S, symmetric = TRUE)
  if (min(ev$values) < eps)
    S <- ev$vectors %*% diag(pmax(ev$values, eps)) %*% t(ev$vectors)
  (S + t(S)) / 2
}

# ---- d-general DGP / W parameters (d=2 reproduces probe 43 exactly) ---------
make_mu <- function(R, d) {
  if (d == 2L && R == 4L) return(c(0, 0.2, 0.5, 0.1, 1, -0.3, 0.7, 0.4))
  m <- numeric(R * d)
  for (r in 1:R) for (a in 1:d)
    m[(r - 1) * d + a] <- 0.3 * (r - 1) + 0.18 * (a - 1) - 0.1 * ((r + a) %% 2)
  m
}
make_psi <- function(d) { M <- matrix(0.4, d, d); diag(M) <- seq(1.0, 0.8, length.out = d); M }
make_E   <- function(d) { M <- matrix(0.1, d, d); diag(M) <- seq(0.6, 0.5, length.out = d); M }
make_W <- function(id, d) {
  if (id == "identity") return(diag(d))
  if (d == 1L) return(matrix(1))
  M <- matrix(0.25, d, d); diag(M) <- seq(1.0, 1.2, length.out = d)
  M[1, 2] <- M[2, 1] <- 0.3
  psd_floor(M, eps = 1e-8)
}

make_dgp <- function(id, d) {
  Psi <- make_psi(d); E <- make_E(d); mu <- make_mu(R, d)
  if (id == "blockcs")
    return(list(mu = mu, Sigma = make_sigma(R, Psi, E), kind = "gauss"))
  if (id == "mixt")
    return(list(mu = mu, Sigma = make_sigma(R, Psi, E), kind = "t6"))
  if (id == "camps") {
    # saturated-but-not-block-CS: unequal rater diagonal blocks (s_r * E) +
    # an extra correlated rater pair (1,2). Block-CS misspecified, saturated OK.
    s <- seq(0.7, 1.6, length.out = R)
    Delta <- diag(seq(0.35, 0.2, length.out = d), d, d)
    S <- kronecker(matrix(1, R, R), Psi)                # every block = Psi
    for (r in 1:R) {
      ix <- ((r - 1) * d + 1):(r * d)
      S[ix, ix] <- S[ix, ix] + s[r] * E
    }
    i1 <- 1:d; i2 <- (d + 1):(2 * d)
    S[i1, i2] <- S[i1, i2] + Delta
    S[i2, i1] <- S[i2, i1] + t(Delta)
    return(list(mu = mu, Sigma = psd_floor(S), kind = "gauss"))
  }
  stop("unknown dgp ", id)
}

# ---- simulator: n x R x d array, rater-major / feature-minor stacking --------
# column (r-1)*d + a of the stacked draw must land in arr[, r, a] to match the
# backend (.kvq_flatten / C++ vector_col). Do NOT rely on array()'s fill order.
sim_arr <- function(n, dgp, R, d) {
  Sig <- if (dgp$kind == "t6") (6 - 2) / 6 * dgp$Sigma else dgp$Sigma
  U <- chol(Sig)                                        # U'U = Sig
  Z <- matrix(rnorm(n * R * d), n, R * d) %*% U
  if (dgp$kind == "t6") Z <- Z / sqrt(rgamma(n, 6 / 2, 6 / 2))  # t_6 scale mixture
  Y <- sweep(Z, 2L, dgp$mu, "+")
  arr <- array(NA_real_, c(n, R, d))
  for (r in 1:R) for (a in 1:d) arr[, r, a] <- Y[, (r - 1) * d + a]
  arr
}

# ---- rater-level amputation (whole d-vector observed or missing) ------------
amputate_arr <- function(arr, mech, rate = 0.30) {
  n <- dim(arr)[1]; R <- dim(arr)[2]
  if (mech == "complete") return(arr)
  if (mech == "mcar") {
    drop <- matrix(runif(n * R) < rate, n, R)
    for (r in 1:R) arr[drop[, r], r, ] <- NA
    return(arr)
  }
  if (mech == "mar") {
    # rater 1 is the always-observed anchor -> co-observation graph stays
    # connected; raters 2..R drop with probability rising in rater-1 feature 1.
    a <- as.numeric(scale(arr[, 1, 1]))
    pr <- plogis(-0.4 + 0.9 * a)
    for (r in 2:R) arr[runif(n) < pr, r, ] <- NA
    return(arr)
  }
  stop("unknown mech ", mech)
}
drop_empty <- function(arr) {
  keep <- apply(arr, 1L, function(s) any(is.finite(s)))
  arr[keep, , , drop = FALSE]
}

# ---- one fit ----------------------------------------------------------------
fit_one <- function(arr, method, W) {
  tryCatch({
    f <- kappa(arr, estimator = method, W = W)
    est <- coef(f); se <- sqrt(diag(vcov(f)))
    cw <- suppressWarnings(confint(f, transform = "none"))
    cf <- suppressWarnings(confint(f, transform = "fisher"))
    cv <- f$moments$converged
    list(ok = TRUE, est = est, se = se, cw = cw, cf = cf,
         conv = if (is.null(cv)) TRUE else isTRUE(cv))
  }, error = function(e) list(ok = FALSE, msg = conditionMessage(e)))
}

# ---- truth -----------------------------------------------------------------
# Conger / Fleiss: exact analytic, from the TRUE (mu, Sigma). For 'camps' this is
# the saturated estimand; by the contraction-preservation argument it is ALSO the
# block-CS estimand, so all three methods share one Conger/Fleiss target.
# Brennan-Prediger: no fixed analytic value (its chance baseline uses the observed
# per-feature range, which drifts with n), so use a large-n Monte-Carlo BP from
# the estimator's own convention -- coverage of BP is therefore indicative only.
truth_key <- function(dgp, d, Wn, coef) paste(dgp, d, Wn, coef, sep = "|")
TRUTH <- new.env(parent = emptyenv())
for (dn in DGPS) for (d in D_GRID) for (Wn in W_GRID) {
  dg <- make_dgp(dn, d); W <- make_W(Wn, d)
  cf2 <- vec_coef(dg$mu, dg$Sigma, R, d, W)
  assign(truth_key(dn, d, Wn, "Conger"), unname(cf2["Conger"]), TRUTH)
  assign(truth_key(dn, d, Wn, "Fleiss"), unname(cf2["Fleiss"]), TRUTH)
  big <- sim_arr(100000L, dg, R, d)                     # MC BP truth (complete)
  bpf <- tryCatch(coef(kappa(big, estimator = "nt_fiml", W = W))[["Brennan-Prediger"]],
                  error = function(e) NA_real_)
  assign(truth_key(dn, d, Wn, "Brennan-Prediger"), bpf, TRUTH)
}

# ---- self-test: stacking order + analytic-vs-fit parity ---------------------
{
  dg <- make_dgp("blockcs", 2L); W <- make_W("identity", 2L)
  set.seed(999)
  big <- sim_arr(60000L, dg, R, 2L)
  fit <- kappa(big, estimator = "nt_fiml", W = W)
  ana <- vec_coef(dg$mu, dg$Sigma, R, 2L, W)
  err <- max(abs(coef(fit)[c("Conger", "Fleiss")] - ana))
  cat(sprintf("self-test  max|fit - analytic| (Conger/Fleiss, n=60k) = %.4f\n", err))
  if (err > 0.02) stop("self-test FAILED: stacking order or contraction mismatch.")
  cat("self-test OK\n\n")
}

# ---- scenario grid ----------------------------------------------------------
scen <- expand.grid(dgp = DGPS, mech = MECHS, n = N_GRID, d = D_GRID,
                    Wn = W_GRID, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

# preallocate replicate-row columns
maxrows <- nrow(scen) * B * length(COEFS3)
R_dgp <- character(maxrows); R_mech <- character(maxrows); R_n <- integer(maxrows)
R_d <- integer(maxrows); R_Wn <- character(maxrows); R_method <- character(maxrows)
R_coef <- character(maxrows); R_rep <- integer(maxrows)
R_est <- numeric(maxrows); R_se <- numeric(maxrows)
R_lw <- numeric(maxrows); R_hw <- numeric(maxrows)
R_lf <- numeric(maxrows); R_hf <- numeric(maxrows)
R_conv <- logical(maxrows)
k <- 0L
push <- function(dgp, mech, n, d, Wn, method, coef, rep, est, se, lw, hw, lf, hf, conv) {
  k <<- k + 1L
  R_dgp[k] <<- dgp; R_mech[k] <<- mech; R_n[k] <<- n; R_d[k] <<- d
  R_Wn[k] <<- Wn; R_method[k] <<- method; R_coef[k] <<- coef; R_rep[k] <<- rep
  R_est[k] <<- est; R_se[k] <<- se
  R_lw[k] <<- lw; R_hw[k] <<- hw; R_lf[k] <<- lf; R_hf[k] <<- hf; R_conv[k] <<- conv
}

# ---- main loop --------------------------------------------------------------
t0 <- proc.time()[["elapsed"]]
for (s in seq_len(nrow(scen))) {
  dn <- scen$dgp[s]; mech <- scen$mech[s]; n <- scen$n[s]
  d <- scen$d[s]; Wn <- scen$Wn[s]
  dg <- make_dgp(dn, d); W <- make_W(Wn, d)
  for (b in seq_len(B)) {
    set.seed(44000L + s * 100000L + b)
    arr <- drop_empty(amputate_arr(sim_arr(n, dg, R, d), mech))
    for (m in METHODS) {
      mc <- coefs_of(m)
      f <- fit_one(arr, m, W)
      if (isTRUE(f$ok)) {
        for (cn in mc) {
          e <- f$est[[cn]]; se <- f$se[[cn]]
          push(dn, mech, n, d, Wn, m, cn, b, e, se,
               f$cw[cn, 1], f$cw[cn, 2], f$cf[cn, 1], f$cf[cn, 2], f$conv)
        }
      } else {
        for (cn in mc)
          push(dn, mech, n, d, Wn, m, cn, b, NA, NA, NA, NA, NA, NA, FALSE)
      }
    }
  }
  cat(sprintf("  cell %2d/%2d  %-8s %-8s n=%-3d d=%d W=%-8s  (%.0fs)\n",
              s, nrow(scen), dn, mech, n, d, Wn, proc.time()[["elapsed"]] - t0))
}
replicates <- data.frame(
  dgp = R_dgp[1:k], mech = R_mech[1:k], n = R_n[1:k], d = R_d[1:k], W = R_Wn[1:k],
  method = R_method[1:k], coefficient = R_coef[1:k], rep = R_rep[1:k],
  est = R_est[1:k], se = R_se[1:k], lo_w = R_lw[1:k], hi_w = R_hw[1:k],
  lo_f = R_lf[1:k], hi_f = R_hf[1:k], conv = R_conv[1:k],
  stringsAsFactors = FALSE)

# ---- aggregate --------------------------------------------------------------
gkey <- with(replicates, paste(dgp, mech, n, d, W, method, coefficient, sep = "|"))
summ <- do.call(rbind, lapply(split(replicates, gkey), function(g) {
  ok <- is.finite(g$est)
  est <- g$est[ok]; se <- g$se[ok]
  t0 <- get(truth_key(g$dgp[1], g$d[1], g$W[1], g$coefficient[1]), TRUTH)
  fin_w <- ok & is.finite(g$lo_w) & is.finite(g$hi_w)
  fin_f <- ok & is.finite(g$lo_f) & is.finite(g$hi_f)
  data.frame(
    dgp = g$dgp[1], mech = g$mech[1], n = g$n[1], d = g$d[1], W = g$W[1],
    method = g$method[1], coefficient = g$coefficient[1], truth = t0,
    conv = mean(is.finite(g$est)),
    bias = mean(est) - t0, sd_emp = sd(est),
    rmse = sqrt((mean(est) - t0)^2 + (length(est) - 1) / length(est) * var(est)),
    mean_se = mean(se),
    se_ratio = mean(se) / sd(est),
    cov_wald = mean(g$lo_w[fin_w] <= t0 & t0 <= g$hi_w[fin_w]),
    cov_fisher = mean(g$lo_f[fin_f] <= t0 & t0 <= g$hi_f[fin_f]),
    stringsAsFactors = FALSE)
}))
ord <- with(summ, order(dgp, factor(mech, MECHS), n, d, W,
                        factor(method, METHODS),
                        factor(coefficient, COEFS3)))
summ <- summ[ord, ]
rownames(summ) <- NULL

# ---- report -----------------------------------------------------------------
pr <- summ
num <- sapply(pr, is.numeric); pr[num] <- lapply(pr[num], function(z) round(z, 4))
cat("\n==== probe 44 summary (Conger/Fleiss truth exact; BP truth MC large-n) ====\n")
print(pr, row.names = FALSE)
saveRDS(summ, file.path(HERE, "coverage_summary.rds"))
write.csv(summ, file.path(HERE, "coverage_summary.csv"), row.names = FALSE)
write.csv(replicates, file.path(HERE, "replicates.csv"), row.names = FALSE)
cat(sprintf("\ndone (%.0fs). wrote coverage_summary.{rds,csv}, replicates.csv\n",
            proc.time()[["elapsed"]] - t0))
