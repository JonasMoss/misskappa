# Probe 44 — STRESS study for the vector-valued quadratic kappa estimators.
# The screen run (parity_and_coverage.R) settled the headline map; this targets
# the three things that can still surprise us:
#
#   Regime 1  ADVERSARIAL MAR.  A two-camp `camps` covariance with a tunable
#             misspecification strength s (s=0 is exact block-CS) plus a MAR
#             mechanism aimed at the cross-camp blocks the CS model averages
#             away. Does nt_fiml_cs ever go inconsistent under MAR (like the
#             variance-component alpha did in probe 40), and does the bias grow
#             with s and persist at large n? Saturated nt_fiml is the consistent
#             yardstick; pairwise is the MAR-broken floor.
#
#   Regime 2  ROBUSTNESS / conv-rate.  Small R + sparse/designed missingness
#             where the saturated pairwise/nt_fiml paths FAIL (missing pair
#             overlap / singular per-pattern covariance) but the connected-only
#             nt_fiml_cs survives. The screen run never triggered this (all
#             conv=1); here conv rate IS the headline metric.
#
#   Regime 3  FIXED-RANGE BP.  With the new `ranges=` hook the Brennan-Prediger
#             chance baseline is deterministic, so BP has a genuine population
#             value (analytic) and its SE should calibrate. Confirms the fix.
#
# Parallel over replicates via mclapply. Usage:
#   Rscript stress.R                 # all three regimes
#   Rscript stress.R 1 3             # only regimes 1 and 3

suppressMessages(library(misskappa))
library(parallel)

HERE <- "experiments/probes/44-vector-quadratic-coverage"
CORES <- max(1L, detectCores() - 1L)
sel <- as.integer(commandArgs(trailingOnly = TRUE))
if (!length(sel)) sel <- 1:3
RNGkind("L'Ecuyer-CMRG")
cat(sprintf("probe 44 stress  cores=%d  regimes={%s}\n\n", CORES,
            paste(sel, collapse = ",")))

# ---- shared building blocks (probe 43 conventions) --------------------------
make_sigma <- function(R, Psi, E)
  kronecker(matrix(1, R, R), Psi) + kronecker(diag(R), E)
psd_floor <- function(S, eps = 1e-6) {
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
  for (r in 1:R) for (a in 1:d)
    m[(r - 1) * d + a] <- 0.3 * (r - 1) + 0.18 * (a - 1) - 0.1 * ((r + a) %% 2)
  m
}
make_W <- function(id, d) {
  if (id == "identity") return(diag(d))
  if (d == 1L) return(matrix(1))
  M <- matrix(0.25, d, d); diag(M) <- seq(1.0, 1.2, length.out = d)
  M[1, 2] <- M[2, 1] <- 0.3
  psd_floor(M, eps = 1e-8)
}

# Conger / Fleiss / Brennan-Prediger from (mu, Sigma, W); BP needs a chance.
vec_coef <- function(mu, Sigma, R, d, W, bp_chance = NA_real_) {
  Bd <- kronecker(matrix(1, R, R), W)
  Td <- kronecker(diag(R), W)
  Gd <- kronecker(diag(R) - matrix(1 / R, R, R), W)
  tB <- sum(Bd * Sigma); tT <- sum(Td * Sigma)
  tG <- as.numeric(t(mu) %*% Gd %*% mu)
  conger <- (tB - tT) / ((R - 1) * tT + R * tG)
  fleiss <- (tB - tT - tG) / ((R - 1) * (tT + tG))
  d_obs <- (2 / (R - 1)) * (tT + tG - tB / R)
  bp <- if (is.na(bp_chance)) NA_real_ else 1 - d_obs / bp_chance
  c(Conger = conger, Fleiss = fleiss, `Brennan-Prediger` = bp)
}

# ---- two-camp adversarial covariance: s=0 is exact block-CS -----------------
# Camps {1,2} and {3,4}. Strength s scales (a) unequal rater diagonal blocks and
# (b) extra WITHIN-camp coupling, both departures the block-CS model cannot see.
camps_sigma <- function(R, d, s) {
  Psi <- make_psi(d); E <- make_E(d)
  spread <- if (R == 4L) c(-0.6, -0.2, 0.2, 0.6) else
            seq(-0.6, 0.6, length.out = R)
  S <- kronecker(matrix(1, R, R), Psi)                 # every block = Psi
  for (r in 1:R) {
    ix <- ((r - 1) * d + 1):(r * d)
    S[ix, ix] <- S[ix, ix] + (1 + s * spread[r]) * E   # s=0 -> +E  (block-CS)
  }
  addb <- function(S, r, c, M) {
    ir <- ((r - 1) * d + 1):(r * d); ic <- ((c - 1) * d + 1):(c * d)
    S[ir, ic] <- S[ir, ic] + M; S
  }
  G <- 0.3 * E
  camps <- if (R >= 4L) list(c(1, 2), c(3, 4)) else list(c(1, 2))
  for (cc in camps) { S <- addb(S, cc[1], cc[2], s * G); S <- addb(S, cc[2], cc[1], s * t(G)) }
  psd_floor(S)
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
# camp-aimed MAR: anchor camp A (raters 1,2 always observed), camp B (raters
# 3..R) drops with probability rising in rater-1 feature-1 -> the cross-camp
# blocks the CS model mis-averages are under-observed in a value-dependent way.
amp_mar_camps <- function(arr) {
  n <- dim(arr)[1]; R <- dim(arr)[2]
  pr <- plogis(-0.4 + 0.9 * as.numeric(scale(arr[, 1, 1])))
  for (r in 3:R) arr[runif(n) < pr, r, ] <- NA
  arr
}
amp_random_k <- function(arr, k) {                     # each subject keeps k raters
  n <- dim(arr)[1]; R <- dim(arr)[2]
  for (i in 1:n) arr[i, -sample.int(R, k), ] <- NA
  arr
}
drop_empty <- function(arr) {
  keep <- apply(arr, 1L, function(s) any(is.finite(s)))
  arr[keep, , , drop = FALSE]
}

# ---- one fit ----------------------------------------------------------------
coefs_of <- function(m) if (m == "pairwise") c("Conger", "Fleiss") else
  c("Conger", "Fleiss", "Brennan-Prediger")
fit_one <- function(arr, method, W, ranges = NULL) {
  tryCatch({
    f <- if (is.null(ranges)) kappa(arr, estimator = method, W = W)
         else kappa(arr, estimator = method, W = W, ranges = ranges)
    list(ok = TRUE, est = coef(f), se = sqrt(diag(vcov(f))),
         cw = suppressWarnings(confint(f, transform = "none")),
         cf = suppressWarnings(confint(f, transform = "fisher")))
  }, error = function(e) list(ok = FALSE))
}

# one replicate -> long data.frame of (method, coefficient, est, se, ci) rows
one_rep <- function(sim_fun, amp_fun, methods, W, score, ranges) {
  arr <- drop_empty(amp_fun(sim_fun()))
  rows <- list()
  for (m in methods) {
    mc <- intersect(coefs_of(m), score)
    if (!length(mc)) next
    f <- fit_one(arr, m, W, ranges)
    for (cn in mc) {
      ok <- isTRUE(f$ok) && cn %in% names(f$est)
      rows[[length(rows) + 1L]] <- data.frame(
        method = m, coefficient = cn,
        est = if (ok) f$est[[cn]] else NA_real_,
        se  = if (ok) f$se[[cn]]  else NA_real_,
        lo_w = if (ok) f$cw[cn, 1] else NA_real_, hi_w = if (ok) f$cw[cn, 2] else NA_real_,
        lo_f = if (ok) f$cf[cn, 1] else NA_real_, hi_f = if (ok) f$cf[cn, 2] else NA_real_,
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

# run B replicates of a cell in parallel; attach truth; aggregate
run_cell <- function(seed, B, sim_fun, amp_fun, methods, W, truth, ranges = NULL) {
  set.seed(seed)
  reps <- mclapply(seq_len(B), function(b)
    tryCatch(one_rep(sim_fun, amp_fun, methods, W, names(truth), ranges),
             error = function(e) NULL),
    mc.cores = CORES, mc.set.seed = TRUE)
  df <- do.call(rbind, reps[!vapply(reps, is.null, logical(1))])
  split(df, paste(df$method, df$coefficient, sep = "|")) |>
    lapply(function(g) {
      t0 <- truth[[g$coefficient[1]]]
      ok <- is.finite(g$est); est <- g$est[ok]; se <- g$se[ok]
      fw <- ok & is.finite(g$lo_w); ff <- ok & is.finite(g$lo_f)
      data.frame(
        method = g$method[1], coefficient = g$coefficient[1], truth = t0,
        conv = mean(is.finite(g$est)),
        bias = if (length(est)) mean(est) - t0 else NA_real_,
        sd_emp = if (length(est) > 1) sd(est) else NA_real_,
        mean_se = mean(se),
        se_ratio = if (length(est) > 1) mean(se) / sd(est) else NA_real_,
        cov_wald = if (any(fw)) mean(g$lo_w[fw] <= t0 & t0 <= g$hi_w[fw]) else NA_real_,
        cov_fisher = if (any(ff)) mean(g$lo_f[ff] <= t0 & t0 <= g$hi_f[ff]) else NA_real_,
        stringsAsFactors = FALSE)
    }) |> (\(x) do.call(rbind, x))()
}
rnd <- function(df) { num <- sapply(df, is.numeric); df[num] <- lapply(df[num], round, 4); df }
R <- 4L; d <- 2L

# =============================================================================
# Regime 1 — adversarial MAR stress
# =============================================================================
if (1 %in% sel) {
  cat("=== Regime 1: adversarial MAR stress (two-camp camps, strength s) ===\n")
  methods <- c("pairwise", "nt_fiml", "nt_fiml_cs")
  grid <- expand.grid(s = c(0, 0.5, 1, 2, 4), n = c(400L, 1500L),
                      mech = c("mcar", "mar"), stringsAsFactors = FALSE)
  B <- 300L; mu <- make_mu(R, d); W <- make_W("identity", d)
  out <- list(); t0 <- proc.time()[["elapsed"]]
  for (i in seq_len(nrow(grid))) {
    s <- grid$s[i]; n <- grid$n[i]; mech <- grid$mech[i]
    Sig <- camps_sigma(R, d, s)
    truth <- vec_coef(mu, Sig, R, d, W)[c("Conger", "Fleiss")]
    amp <- if (mech == "mcar") function(a) amp_mcar(a, 0.30) else amp_mar_camps
    res <- run_cell(10000L + i, B, function() sim_arr(n, mu, Sig, R, d),
                    amp, methods, W, truth)
    res <- cbind(s = s, n = n, mech = mech, res)
    out[[i]] <- res
    cat(sprintf("  s=%-3g n=%-4d %-4s  (%.0fs)\n", s, n, mech,
                proc.time()[["elapsed"]] - t0))
  }
  r1 <- do.call(rbind, out)
  saveRDS(r1, file.path(HERE, "stress_regime1.rds"))
  cat("\n-- nt_fiml_cs vs nt_fiml: Conger bias by strength (MAR) --\n")
  v <- r1[r1$coefficient == "Conger" & r1$mech == "mar" &
            r1$method %in% c("nt_fiml", "nt_fiml_cs"),
          c("s", "n", "method", "bias", "cov_wald")]
  print(rnd(v[order(v$n, v$s, v$method), ]), row.names = FALSE)
  cat("\n(MCAR control: nt_fiml_cs should stay ~unbiased at ALL s)\n")
  vc <- r1[r1$coefficient == "Conger" & r1$mech == "mcar" &
             r1$method == "nt_fiml_cs" & r1$n == 1500, c("s", "bias", "cov_wald")]
  print(rnd(vc), row.names = FALSE)
  cat("\n")
}

# =============================================================================
# Regime 2 — robustness / conv-rate (clean block-CS DGP)
# =============================================================================
if (2 %in% sel) {
  cat("=== Regime 2: robustness / conv-rate (sparse designs, small R) ===\n")
  methods <- c("pairwise", "nt_fiml", "nt_fiml_cs")
  grid <- expand.grid(R = c(3L, 4L, 5L), design = c("random2", "mcar70"),
                      n = c(40L, 80L, 150L), stringsAsFactors = FALSE)
  B <- 400L; W <- make_W("identity", d)
  out <- list(); t0 <- proc.time()[["elapsed"]]
  for (i in seq_len(nrow(grid))) {
    Rg <- grid$R[i]; design <- grid$design[i]; n <- grid$n[i]
    mu <- make_mu(Rg, d); Sig <- make_sigma(Rg, make_psi(d), make_E(d))
    truth <- vec_coef(mu, Sig, Rg, d, W)[c("Conger", "Fleiss")]
    amp <- if (design == "random2") function(a) amp_random_k(a, 2L) else
                                    function(a) amp_mcar(a, 0.70)
    res <- run_cell(20000L + i, B, function() sim_arr(n, mu, Sig, Rg, d),
                    amp, methods, W, truth)
    res <- cbind(R = Rg, design = design, n = n, res)
    out[[i]] <- res
    cat(sprintf("  R=%d %-8s n=%-3d  (%.0fs)\n", Rg, design, n,
                proc.time()[["elapsed"]] - t0))
  }
  r2 <- do.call(rbind, out)
  saveRDS(r2, file.path(HERE, "stress_regime2.rds"))
  cat("\n-- conv rate by method (Conger row) --\n")
  v <- r2[r2$coefficient == "Conger", c("R", "design", "n", "method", "conv")]
  w <- reshape(v, idvar = c("R", "design", "n"), timevar = "method",
               direction = "wide")
  names(w) <- sub("conv.", "", names(w), fixed = TRUE)
  print(rnd(w[order(w$design, w$R, w$n), ]), row.names = FALSE)
  cat("\n")
}

# =============================================================================
# Regime 3 — fixed-range Brennan-Prediger calibration
# =============================================================================
if (3 %in% sel) {
  cat("=== Regime 3: fixed-range BP calibration (ranges=) ===\n")
  methods <- c("nt_fiml", "nt_fiml_cs")
  ranges <- rbind(c(-6, 6), c(-5, 5)); widths <- ranges[, 2] - ranges[, 1]
  grid <- expand.grid(Wn = c("identity", "coupling"),
                      mech = c("complete", "mcar", "mar"),
                      n = c(120L, 300L, 1000L), stringsAsFactors = FALSE)
  B <- 400L; mu <- make_mu(R, d); Sig <- make_sigma(R, make_psi(d), make_E(d))
  out <- list(); t0 <- proc.time()[["elapsed"]]
  for (i in seq_len(nrow(grid))) {
    Wn <- grid$Wn[i]; mech <- grid$mech[i]; n <- grid$n[i]
    W <- make_W(Wn, d)
    bp_chance <- sum(diag(W) * widths^2 / 6)            # deterministic baseline
    truth <- vec_coef(mu, Sig, R, d, W, bp_chance)      # all three, BP analytic
    amp <- switch(mech, complete = identity,
                  mcar = function(a) amp_mcar(a, 0.30),
                  mar = function(a) {
                    pr <- plogis(-0.4 + 0.9 * as.numeric(scale(a[, 1, 1])))
                    for (r in 2:R) a[runif(dim(a)[1]) < pr, r, ] <- NA; a
                  })
    res <- run_cell(30000L + i, B, function() sim_arr(n, mu, Sig, R, d),
                    amp, methods, W, truth, ranges = ranges)
    res <- cbind(W = Wn, mech = mech, n = n, res)
    out[[i]] <- res
    cat(sprintf("  W=%-8s %-8s n=%-4d  (%.0fs)\n", Wn, mech, n,
                proc.time()[["elapsed"]] - t0))
  }
  r3 <- do.call(rbind, out)
  saveRDS(r3, file.path(HERE, "stress_regime3.rds"))
  cat("\n-- Brennan-Prediger calibration with fixed range (was se_ratio~0.3, cov 0) --\n")
  v <- r3[r3$coefficient == "Brennan-Prediger",
          c("W", "mech", "n", "method", "bias", "se_ratio", "cov_wald")]
  print(rnd(v[order(v$W, v$mech, v$n, v$method), ]), row.names = FALSE)
  cat("\n")
}

cat("done.\n")
