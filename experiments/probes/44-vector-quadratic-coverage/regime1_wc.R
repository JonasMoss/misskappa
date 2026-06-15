# Probe 44 — REGIME 1 REBUILD with a WELL-CONDITIONED adversarial knob.
#
# The committed stress.R regime 1 used a two-camp `camps_sigma` whose diagonal
# rater-variance blocks go NEGATIVE for s > 1.67 (1 + s*spread[r], spread in
# [-0.6,0.6]), so it gets psd_floor()ed into a near-singular truth at s >= 2.
# That confounded two distinct stories: (i) nt_fiml_cs's MAR bias from block-CS
# misspecification, and (ii) saturated nt_fiml "failing to converge" — which at
# s >= 2 was really "the TRUTH is degenerate", not "missingness broke the EM".
#
# This rebuild uses a latent CAMP-FACTOR covariance that stays positive definite
# at every s, with the smallest eigenvalue PINNED at lambda_min(E):
#
#   Y_ij = mu_j + b_i + c_i^{camp(j)} + eps_ij
#     b_i        ~ N(0, Psi)     shared subject factor      -> J_R (x) Psi
#     c_i^{camp} ~ N(0, s*G)     per-camp factor            -> (1_camp 1_camp') (x) sG
#     eps_ij     ~ N(0, E)       rater residual             -> I_R (x) E
#
#   Sigma = J_R(x)Psi + I_R(x)E + sum_camp (1_camp 1_camp')(x)(s*G)
#     diagonal block      = Psi + s*G + E   (CONSTANT across raters, PD)
#     within-camp off-diag = Psi + s*G
#     cross-camp off-diag  = Psi
#
# s = 0 is exact block-CS; s > 0 makes block-CS FALSE in a graded way (within-camp
# pairs share an extra factor the CS model averages away) while Sigma_oo >= I(x)E
# for every observed sub-pattern, so the saturated EM's per-pattern inversions are
# never population-singular. This isolates the CS bias curve from DGP degeneracy
# AND tests whether saturated's non-convergence was an artifact of the old knob.
#
# Methods: pairwise (MCAR floor), nt_fiml (saturated, consistent yardstick),
#          nt_fiml_cs (block-CS).  Metric headline: Conger bias/cov vs s, PLUS
#          the conv rate of saturated nt_fiml at high s (the disentangling test).
#
# Usage:  Rscript experiments/probes/44-vector-quadratic-coverage/regime1_wc.R

suppressMessages(library(misskappa))
library(parallel)

HERE  <- "experiments/probes/44-vector-quadratic-coverage"
CORES <- max(1L, detectCores() - 1L)
RNGkind("L'Ecuyer-CMRG")
R <- 4L; d <- 2L

# ---- building blocks (probe 43 / stress.R conventions) ----------------------
psd_floor <- function(S, eps = 1e-8) {
  S <- (S + t(S)) / 2
  ev <- eigen(S, symmetric = TRUE)
  if (min(ev$values) < eps)
    S <- ev$vectors %*% diag(pmax(ev$values, eps)) %*% t(ev$vectors)
  (S + t(S)) / 2
}
make_psi <- function(d) { M <- matrix(0.4, d, d); diag(M) <- seq(1.0, 0.8, length.out = d); M }
make_E   <- function(d) { M <- matrix(0.1, d, d); diag(M) <- seq(0.6, 0.5, length.out = d); M }
make_mu  <- function(R, d) {
  if (d == 2L && R == 4L) return(c(0, 0.2, 0.5, 0.1, 1, -0.3, 0.7, 0.4))
  m <- numeric(R * d)
  for (r in 1:R) for (a in 1:d) m[(r - 1) * d + a] <- 0.3 * (r - 1) + 0.18 * (a - 1) - 0.1 * ((r + a) %% 2)
  m
}
make_W <- function(id, d) {
  if (id == "identity") return(diag(d))
  if (d == 1L) return(matrix(1))
  M <- matrix(0.25, d, d); diag(M) <- seq(1.0, 1.2, length.out = d)
  M[1, 2] <- M[2, 1] <- 0.3
  psd_floor(M, eps = 1e-8)
}

vec_coef <- function(mu, Sigma, R, d, W, bp_chance = NA_real_) {
  Bd <- kronecker(matrix(1, R, R), W); Td <- kronecker(diag(R), W)
  Gd <- kronecker(diag(R) - matrix(1 / R, R, R), W)
  tB <- sum(Bd * Sigma); tT <- sum(Td * Sigma)
  tG <- as.numeric(t(mu) %*% Gd %*% mu)
  conger <- (tB - tT) / ((R - 1) * tT + R * tG)
  fleiss <- (tB - tT - tG) / ((R - 1) * (tT + tG))
  c(Conger = conger, Fleiss = fleiss)
}

# ---- WELL-CONDITIONED two-camp covariance: s=0 is exact block-CS -------------
# camp membership: {1,2} and {3,4} for R=4. G is the per-camp factor covariance.
camp_of <- function(R) if (R == 4L) c(1, 1, 2, 2) else rep(1:2, length.out = R)
camps_wc <- function(R, d, s) {
  Psi <- make_psi(d); E <- make_E(d); G <- 0.6 * make_psi(d)
  camp <- camp_of(R)
  S <- kronecker(matrix(1, R, R), Psi) + kronecker(diag(R), E)   # block-CS base
  for (cmp in unique(camp)) {
    m  <- which(camp == cmp)
    Mc <- matrix(0, R, R); Mc[m, m] <- 1                          # 1_camp 1_camp'
    S  <- S + kronecker(Mc, s * G)
  }
  (S + t(S)) / 2                                                  # PD by construction
}

# ---- simulate / amputate ----------------------------------------------------
sim_arr <- function(n, mu, Sigma, R, d) {
  U <- chol(Sigma)
  Y <- sweep(matrix(rnorm(n * R * d), n, R * d) %*% U, 2L, mu, "+")
  arr <- array(NA_real_, c(n, R, d))
  for (r in 1:R) for (a in 1:d) arr[, r, a] <- Y[, (r - 1) * d + a]
  arr
}
amp_mcar <- function(arr, rate) {
  n <- dim(arr)[1]; R <- dim(arr)[2]
  drop <- matrix(runif(n * R) < rate, n, R)
  for (r in 1:R) arr[drop[, r], r, ] <- NA
  arr
}
# camp-aimed MAR: keep camp A (raters 1,2) always; drop camp B (raters 3..R)
# with probability rising in rater-1 feature-1 -> the cross-camp blocks the CS
# model mis-averages are under-observed in a value-dependent way.
amp_mar_camps <- function(arr) {
  n <- dim(arr)[1]; R <- dim(arr)[2]
  pr <- plogis(-0.4 + 0.9 * as.numeric(scale(arr[, 1, 1])))
  for (r in 3:R) arr[runif(n) < pr, r, ] <- NA
  arr
}
drop_empty <- function(arr) {
  keep <- apply(arr, 1L, function(s) any(is.finite(s)))
  arr[keep, , , drop = FALSE]
}

# ---- one fit / one replicate / one cell -------------------------------------
coefs_of <- function(m) c("Conger", "Fleiss")
fit_one <- function(arr, method, W) tryCatch({
  f <- kappa(arr, estimator = method, W = W)
  list(ok = TRUE, est = coef(f), se = sqrt(diag(vcov(f))),
       cw = suppressWarnings(confint(f, transform = "none")))
}, error = function(e) list(ok = FALSE))

one_rep <- function(sim_fun, amp_fun, methods, W) {
  arr <- drop_empty(amp_fun(sim_fun()))
  rows <- list()
  for (m in methods) {
    f <- fit_one(arr, m, W)
    for (cn in coefs_of(m)) {
      ok <- isTRUE(f$ok) && cn %in% names(f$est)
      rows[[length(rows) + 1L]] <- data.frame(
        method = m, coefficient = cn,
        est = if (ok) f$est[[cn]] else NA_real_,
        se  = if (ok) f$se[[cn]]  else NA_real_,
        lo_w = if (ok) f$cw[cn, 1] else NA_real_, hi_w = if (ok) f$cw[cn, 2] else NA_real_,
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

run_cell <- function(seed, B, sim_fun, amp_fun, methods, W, truth) {
  set.seed(seed)
  reps <- mclapply(seq_len(B), function(b)
    tryCatch(one_rep(sim_fun, amp_fun, methods, W), error = function(e) NULL),
    mc.cores = CORES, mc.set.seed = TRUE)
  df <- do.call(rbind, reps[!vapply(reps, is.null, logical(1))])
  split(df, paste(df$method, df$coefficient, sep = "|")) |>
    lapply(function(g) {
      t0 <- truth[[g$coefficient[1]]]
      ok <- is.finite(g$est); est <- g$est[ok]; se <- g$se[ok]
      fw <- ok & is.finite(g$lo_w)
      data.frame(
        method = g$method[1], coefficient = g$coefficient[1], truth = t0,
        conv = mean(is.finite(g$est)),
        bias = if (length(est)) mean(est) - t0 else NA_real_,
        sd_emp = if (length(est) > 1) sd(est) else NA_real_,
        mean_se = mean(se),
        se_ratio = if (length(est) > 1) mean(se) / sd(est) else NA_real_,
        cov_wald = if (any(fw)) mean(g$lo_w[fw] <= t0 & t0 <= g$hi_w[fw]) else NA_real_,
        stringsAsFactors = FALSE)
    }) |> (\(x) do.call(rbind, x))()
}
rnd <- function(df) { num <- sapply(df, is.numeric); df[num] <- lapply(df[num], round, 4); df }

# ---- prove the knob is well-conditioned (vs old psd-floored knob) -----------
cat(sprintf("probe 44 regime-1 REBUILD (well-conditioned knob)  cores=%d\n\n", CORES))
cat("-- conditioning of the camp-factor truth (min eig PINNED at lambda_min(E)) --\n")
cnd <- do.call(rbind, lapply(c(0, 0.5, 1, 2, 4), function(s) {
  S <- camps_wc(R, d, s); ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  data.frame(s = s, min_eig = min(ev), max_eig = max(ev), rcond = min(ev) / max(ev))
}))
print(rnd(cnd), row.names = FALSE)

# ---- run grid ---------------------------------------------------------------
methods <- c("pairwise", "nt_fiml", "nt_fiml_cs")
grid <- expand.grid(s = c(0, 0.5, 1, 2, 4), n = c(400L, 1500L),
                    mech = c("mcar", "mar"), stringsAsFactors = FALSE)
B <- 300L; mu <- make_mu(R, d); W <- make_W("identity", d)
out <- list(); t0 <- proc.time()[["elapsed"]]
for (i in seq_len(nrow(grid))) {
  s <- grid$s[i]; n <- grid$n[i]; mech <- grid$mech[i]
  Sig <- camps_wc(R, d, s)
  truth <- vec_coef(mu, Sig, R, d, W)
  amp <- if (mech == "mcar") function(a) amp_mcar(a, 0.30) else amp_mar_camps
  res <- run_cell(40000L + i, B, function() sim_arr(n, mu, Sig, R, d), amp, methods, W, truth)
  out[[i]] <- cbind(s = s, n = n, mech = mech, res)
  cat(sprintf("  s=%-3g n=%-4d %-4s  (%.0fs)\n", s, n, mech, proc.time()[["elapsed"]] - t0))
}
r1 <- do.call(rbind, out)
saveRDS(r1, file.path(HERE, "stress_regime1_wc.rds"))

# ---- report -----------------------------------------------------------------
cat("\n== MAR: Conger bias / coverage by strength (the CS bias curve, now clean) ==\n")
v <- r1[r1$coefficient == "Conger" & r1$mech == "mar", ]
print(rnd(v[order(v$n, v$s, v$method), c("s", "n", "method", "conv", "bias", "cov_wald")]),
      row.names = FALSE)

cat("\n== DISENTANGLE: saturated nt_fiml conv rate by s (was 0.00-0.40 at s>=2 on old knob) ==\n")
vc <- r1[r1$coefficient == "Conger" & r1$method == "nt_fiml", ]
print(rnd(vc[order(vc$mech, vc$n, vc$s), c("s", "n", "mech", "conv", "bias", "cov_wald")]),
      row.names = FALSE)

cat("\n== MCAR control: nt_fiml_cs should stay ~unbiased at ALL s (contraction-invariance) ==\n")
vm <- r1[r1$coefficient == "Conger" & r1$mech == "mcar" & r1$method == "nt_fiml_cs", ]
print(rnd(vm[order(vm$n, vm$s), c("s", "n", "conv", "bias", "cov_wald")]), row.names = FALSE)

cat("\ndone.\n")
