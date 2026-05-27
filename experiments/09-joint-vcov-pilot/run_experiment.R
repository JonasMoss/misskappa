#!/usr/bin/env Rscript
#
# 09-joint-vcov-pilot
#
# Pilot for joint inference across non-independent kappa estimates on the
# same data ("Case 3" of the equality-of-two-kappas question — see plan
# file mutable-rolling-wombat). Pure R; no C++ changes. The joint vcov
# is assembled from a nonparametric bootstrap over subjects.
#
# Four sub-questions:
#   UC8 - Fleiss vs Brennan-Prediger contrast on the same fit
#         (analytic; uses existing 3x3 vcov directly).
#   UC1 - All pairwise Cohen kappas + joint vcov; Wald test of pairwise
#         homogeneity. Size under exchangeable raters, power against
#         one miscalibrated rater.
#   UC4 - Joint vcov of (available-case, IPW) + Hausman-style Wald test
#         of equality. Size under MCAR + exchangeable, power under MCAR
#         + non-exchangeable (the experiment-03 DGP).
#   UC2 - Joint vcov of (identity, linear, quadratic) weighted kappas;
#         report joint CIs and bootstrap correlations.
#
# Outputs (under results/):
#   uc8_replicates.csv   per-replicate Fleiss/BP estimates + contrast z
#   uc8_summary.csv      cell-level rejection rate at 0.05
#   uc1_replicates.csv   per-replicate pairwise kappas + Wald stat
#   uc1_summary.csv      cell-level size / power at 0.05
#   uc4_replicates.csv   per-replicate AC + IPW + Wald z for contrast
#   uc4_summary.csv      cell-level size / power at 0.05
#   uc2_replicates.csv   per-replicate (identity, linear, quadratic) kappas
#                        with joint vcov entries
#   uc2_summary.csv      cell-level: mean/SD of each estimate, mean cor matrix
#   metadata.csv         run metadata

suppressPackageStartupMessages({
  library(misskappa)
})

# ---- CLI ----------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
has_flag <- function(name) name %in% args
get_val  <- function(name, default, parser = as.character) {
  i <- which(args == name)
  if (length(i) == 0L) return(default)
  parser(args[i + 1L])
}

if (has_flag("--help") || has_flag("-h")) {
  cat("Usage: Rscript run_experiment.R [options]\n",
      " --smoke         Fast smoke run (reps=20, boots=50).\n",
      " --reps N        MC replicates per cell (default 200).\n",
      " --boots B       Bootstrap reps for joint vcov (default 200).\n",
      " --n N           Subjects per replicate (default 500).\n",
      " --seed-base K   Seed base (default 1).\n",
      " --only UC       Run only one UC: uc8|uc1|uc4|uc2 (default all).\n",
      " --help, -h      This help.\n", sep = "")
  quit("no", status = 0)
}

smoke      <- has_flag("--smoke")
seed_base  <- get_val("--seed-base", 1L, as.integer)
reps_user  <- get_val("--reps", NA_integer_, as.integer)
boots_user <- get_val("--boots", NA_integer_, as.integer)
n_user     <- get_val("--n", NA_integer_, as.integer)
only_user  <- get_val("--only", NA_character_, as.character)

reps  <- if (smoke) 20L else 200L
boots <- if (smoke) 50L else 200L
n     <- 500L
if (!is.na(reps_user))  reps  <- reps_user
if (!is.na(boots_user)) boots <- boots_user
if (!is.na(n_user))     n     <- n_user

ucs <- if (is.na(only_user)) c("uc8", "uc1", "uc4", "uc2") else only_user

results_dir <- "results"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

t0 <- Sys.time()

# ---- DGP helpers --------------------------------------------------------
C <- 5L

simulate_truth_guess <- function(n, R, rho_vec, p_truth, guess_mat) {
  truth <- sample.int(C, n, replace = TRUE, prob = p_truth)
  x_star <- matrix(NA_integer_, nrow = n, ncol = R)
  for (j in seq_len(R)) {
    correct <- stats::rbinom(n, 1L, prob = rho_vec[j]) == 1L
    guessed <- sample.int(C, n, replace = TRUE, prob = guess_mat[j, ])
    x_star[, j] <- ifelse(correct, truth, guessed)
  }
  x_star
}

apply_mcar_independent <- function(x, pi_vec) {
  R <- ncol(x)
  M <- matrix(stats::runif(nrow(x) * R), nrow = nrow(x), ncol = R) <
       matrix(pi_vec, nrow = nrow(x), ncol = R, byrow = TRUE)
  x[!M] <- NA_integer_
  x
}

# Robust quadratic form: handles a possibly-rank-deficient covariance via
# Moore-Penrose pseudoinverse. Returns (stat, df) for the chi-sq test.
wald_quadform <- function(theta, V, A = diag(length(theta)), rcond_tol = 1e-8) {
  d <- A %*% theta
  M <- A %*% V %*% t(A)
  sv <- svd(M)
  keep <- sv$d > sv$d[1] * rcond_tol
  M_pinv <- sv$v[, keep, drop = FALSE] %*%
            diag(1 / sv$d[keep], nrow = sum(keep)) %*%
            t(sv$u[, keep, drop = FALSE])
  stat <- as.numeric(t(d) %*% M_pinv %*% d)
  list(stat = stat, df = sum(keep))
}

# Bootstrap helper. fit_fn(x) returns a numeric vector (one element per
# scalar moment we want a joint vcov for). Returns the boot * K matrix.
bootstrap_fits <- function(x, fit_fn, B, seed_offset) {
  n <- nrow(x)
  ref <- fit_fn(x)
  K <- length(ref)
  out <- matrix(NA_real_, nrow = B, ncol = K)
  for (b in seq_len(B)) {
    set.seed(seed_base + 7919L * seed_offset + b)
    idx <- sample.int(n, n, replace = TRUE)
    val <- try(fit_fn(x[idx, , drop = FALSE]), silent = TRUE)
    if (!inherits(val, "try-error") && length(val) == K) out[b, ] <- val
  }
  attr(out, "point") <- ref
  out
}

# Tag a chunk of timing output.
tic <- function(tag) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), tag))

# ===================== UC8: Fleiss vs BP free contrast ===================
#
# Same fit, two coefficients with known covariance. Population contrast
# delta = kappa_F - BP is zero iff the chance baselines coincide, which
# happens iff rater marginals are exactly uniform (under identity loss).
# Under skewed marginals delta != 0 and the contrast has power.

run_uc8 <- function() {
  tic("UC8 start")
  R <- 6L
  rho_vec <- rep(0.85, R)
  marg_levels <- list(
    uniform = rep(1 / C, C),
    skewed  = c(0.05, 0.10, 0.20, 0.30, 0.35)
  )
  z975 <- qnorm(0.975)

  per_rep <- list()
  for (m_idx in seq_along(marg_levels)) {
    marg_name <- names(marg_levels)[m_idx]
    p_truth   <- marg_levels[[m_idx]]
    guess_mat <- matrix(rep(p_truth, each = R), nrow = R, byrow = TRUE)
    for (b in seq_len(reps)) {
      set.seed(seed_base + 101000L * m_idx + b)
      x <- simulate_truth_guess(n, R, rho_vec, p_truth, guess_mat)
      k <- misskappa::kappa(x, method = "available", weight = "identity")
      est <- coef(k)
      V <- vcov(k)
      c_vec <- c(0, 1, -1)
      delta <- as.numeric(crossprod(c_vec, est))
      v <- as.numeric(crossprod(c_vec, V %*% c_vec))
      se <- sqrt(max(v, 0))
      z <- if (se > 0) delta / se else NA_real_
      per_rep[[length(per_rep) + 1L]] <- data.frame(
        marg = marg_name, b = b,
        kappa_F = est[["Fleiss"]], kappa_BP = est[["Brennan-Prediger"]],
        delta = delta, se = se, z = z,
        rej_05 = if (is.na(z)) NA else abs(z) > z975
      )
    }
  }
  rep_df <- do.call(rbind, per_rep)
  write.csv(rep_df, file.path(results_dir, "uc8_replicates.csv"), row.names = FALSE)

  summ <- do.call(rbind, by(rep_df, rep_df$marg, function(d) {
    data.frame(
      marg       = unique(d$marg),
      reps       = nrow(d),
      mean_kappa_F = mean(d$kappa_F),
      mean_kappa_BP = mean(d$kappa_BP),
      mean_delta = mean(d$delta),
      sd_delta   = sd(d$delta),
      mean_se    = mean(d$se),
      reject_05  = mean(d$rej_05, na.rm = TRUE)
    )
  }))
  rownames(summ) <- NULL
  write.csv(summ, file.path(results_dir, "uc8_summary.csv"), row.names = FALSE)
  tic("UC8 done")
}

# ===================== UC1: pairwise Cohen joint vcov ====================
#
# All choose(R,2) pairwise Cohen kappas + joint vcov via bootstrap.
# Wald test of pairwise homogeneity uses contrasts kappa_p - kappa_1 for
# p = 2..K.

pairwise_kappas <- function(x) {
  R <- ncol(x)
  pairs <- combn(R, 2)
  out <- numeric(ncol(pairs))
  for (p in seq_len(ncol(pairs))) {
    r <- pairs[1, p]; s <- pairs[2, p]
    xrs <- x[, c(r, s), drop = FALSE]
    keep <- rowSums(!is.na(xrs)) == 2L
    if (sum(keep) < 5L) { out[p] <- NA_real_; next }
    fit <- try(
      misskappa::kappa(xrs[keep, , drop = FALSE],
                       method = "available", weight = "identity"),
      silent = TRUE
    )
    if (inherits(fit, "try-error")) { out[p] <- NA_real_; next }
    out[p] <- as.numeric(coef(fit)[["Conger"]])
  }
  out
}

run_uc1 <- function() {
  tic("UC1 start")
  R <- 6L
  p_truth <- c(0.05, 0.10, 0.20, 0.30, 0.35)
  guess_mat <- matrix(rep(p_truth, each = R), nrow = R, byrow = TRUE)
  cells <- list(
    null    = rep(0.92, R),                          # all raters equal
    one_bad = c(0.92, 0.92, 0.92, 0.92, 0.92, 0.60)  # rater 6 miscalibrated
  )
  n_pairs <- choose(R, 2)
  pair_names <- apply(combn(R, 2), 2,
                      function(rs) sprintf("k_%d_%d", rs[1], rs[2]))

  # Contrast matrix: zero out kappa_1; rows give kappa_p - kappa_1.
  A <- cbind(rep(-1, n_pairs - 1L), diag(n_pairs - 1L))

  per_rep <- list()
  for (cell_idx in seq_along(cells)) {
    cell_name <- names(cells)[cell_idx]
    rho_vec <- cells[[cell_idx]]
    for (b in seq_len(reps)) {
      set.seed(seed_base + 202000L * cell_idx + b)
      x <- simulate_truth_guess(n, R, rho_vec, p_truth, guess_mat)
      boot <- bootstrap_fits(x, pairwise_kappas, boots, seed_offset = b)
      pt <- attr(boot, "point")
      ok <- apply(is.finite(boot), 2, all) & is.finite(pt)
      if (!all(ok)) {
        # If any pair lost all bootstrap fits or fails at point, skip.
        next
      }
      V <- cov(boot)
      w <- wald_quadform(pt, V, A)
      p_chi <- pchisq(w$stat, df = w$df, lower.tail = FALSE)
      row <- as.data.frame(t(setNames(pt, pair_names)))
      row$cell <- cell_name
      row$b <- b
      row$wald_stat <- w$stat
      row$wald_df <- w$df
      row$p_chi <- p_chi
      row$reject_05 <- p_chi < 0.05
      per_rep[[length(per_rep) + 1L]] <- row
    }
  }
  rep_df <- do.call(rbind, per_rep)
  write.csv(rep_df, file.path(results_dir, "uc1_replicates.csv"), row.names = FALSE)

  summ <- do.call(rbind, by(rep_df, rep_df$cell, function(d) {
    pair_cols <- grep("^k_", colnames(d), value = TRUE)
    means <- sapply(pair_cols, function(cn) mean(d[[cn]], na.rm = TRUE))
    data.frame(
      cell = unique(d$cell),
      reps = nrow(d),
      mean_min_kappa = min(means),
      mean_max_kappa = max(means),
      mean_range     = max(means) - min(means),
      mean_wald      = mean(d$wald_stat, na.rm = TRUE),
      median_wald    = median(d$wald_stat, na.rm = TRUE),
      reject_05      = mean(d$reject_05, na.rm = TRUE)
    )
  }))
  rownames(summ) <- NULL
  write.csv(summ, file.path(results_dir, "uc1_summary.csv"), row.names = FALSE)
  tic("UC1 done")
}

# ===================== UC4: Hausman AC vs IPW ============================
#
# Joint vcov of (kappa_AC, kappa_IPW) on the same incomplete data via
# bootstrap. Single contrast c = (1, -1) gives the Hausman-style Wald.

ac_ipw_kappas <- function(x) {
  ac <- try(misskappa::kappa(x, method = "available", weight = "identity"),
            silent = TRUE)
  ipw <- try(misskappa::kappa(x, method = "ipw", weight = "identity"),
             silent = TRUE)
  if (inherits(ac, "try-error") || inherits(ipw, "try-error")) {
    return(c(NA_real_, NA_real_))
  }
  c(as.numeric(coef(ac)[["Conger"]]),
    as.numeric(coef(ipw)[["Conger"]]))
}

run_uc4 <- function() {
  tic("UC4 start")
  R <- 6L
  p_truth <- c(0.05, 0.10, 0.20, 0.30, 0.35)
  guess_mat <- matrix(rep(p_truth, each = R), nrow = R, byrow = TRUE)

  # Cell 1 (size): exchangeable raters + uniform pi.
  # Cell 2 (power): non-exchangeable raters + varying pi (the experiment-03
  # cell that drives AC bias of ~0.085).
  cells <- list(
    null    = list(rho = rep(0.92, R),
                   pi  = rep(0.6, R)),
    nonexch = list(rho = c(0.98, 0.96, 0.92, 0.85, 0.72, 0.60),
                   pi  = 0.6 + c(0.35, 0.25, 0.15, -0.15, -0.25, -0.35))
  )

  z975 <- qnorm(0.975)
  per_rep <- list()
  for (cell_idx in seq_along(cells)) {
    cell_name <- names(cells)[cell_idx]
    cfg <- cells[[cell_idx]]
    for (b in seq_len(reps)) {
      set.seed(seed_base + 303000L * cell_idx + b)
      x_star <- simulate_truth_guess(n, R, cfg$rho, p_truth, guess_mat)
      x      <- apply_mcar_independent(x_star, cfg$pi)
      boot   <- bootstrap_fits(x, ac_ipw_kappas, boots, seed_offset = b)
      pt <- attr(boot, "point")
      ok <- apply(is.finite(boot), 2, all) & is.finite(pt)
      if (!all(ok)) next
      V <- cov(boot)
      delta <- pt[1] - pt[2]
      v <- V[1, 1] + V[2, 2] - 2 * V[1, 2]
      se <- sqrt(max(v, 0))
      z <- if (se > 0) delta / se else NA_real_
      per_rep[[length(per_rep) + 1L]] <- data.frame(
        cell = cell_name, b = b,
        kappa_AC = pt[1], kappa_IPW = pt[2],
        delta = delta, se = se, z = z,
        boot_var_AC  = V[1, 1], boot_var_IPW = V[2, 2], boot_cov = V[1, 2],
        reject_05 = if (is.na(z)) NA else abs(z) > z975
      )
    }
  }
  rep_df <- do.call(rbind, per_rep)
  write.csv(rep_df, file.path(results_dir, "uc4_replicates.csv"), row.names = FALSE)

  summ <- do.call(rbind, by(rep_df, rep_df$cell, function(d) {
    data.frame(
      cell        = unique(d$cell),
      reps        = nrow(d),
      mean_AC     = mean(d$kappa_AC),
      mean_IPW    = mean(d$kappa_IPW),
      mean_delta  = mean(d$delta),
      sd_delta_mc = sd(d$delta),
      mean_se_boot = mean(d$se, na.rm = TRUE),
      mean_boot_cor = mean(d$boot_cov / sqrt(pmax(d$boot_var_AC * d$boot_var_IPW, 0)),
                          na.rm = TRUE),
      reject_05   = mean(d$reject_05, na.rm = TRUE)
    )
  }))
  rownames(summ) <- NULL
  write.csv(summ, file.path(results_dir, "uc4_summary.csv"), row.names = FALSE)
  tic("UC4 done")
}

# ===================== UC2: weight-scheme sensitivity ====================
#
# Three weighted kappas on the same complete data. The three estimands
# differ in population, so there is no natural null. Report joint
# distribution and bootstrap correlations between schemes.

three_weight_kappas <- function(x) {
  k_id <- try(misskappa::kappa(x, method = "available", weight = "identity"),
              silent = TRUE)
  k_li <- try(misskappa::kappa(x, method = "available", weight = "linear"),
              silent = TRUE)
  k_qu <- try(misskappa::kappa(x, method = "available", weight = "quadratic"),
              silent = TRUE)
  if (inherits(k_id, "try-error") ||
      inherits(k_li, "try-error") ||
      inherits(k_qu, "try-error")) return(rep(NA_real_, 3L))
  c(as.numeric(coef(k_id)[["Conger"]]),
    as.numeric(coef(k_li)[["Conger"]]),
    as.numeric(coef(k_qu)[["Conger"]]))
}

run_uc2 <- function() {
  tic("UC2 start")
  R <- 6L
  # Ordinal-style skewed marginals so weighted schemes diverge.
  p_truth <- c(0.05, 0.10, 0.20, 0.30, 0.35)
  rho_vec <- rep(0.85, R)
  guess_mat <- matrix(rep(p_truth, each = R), nrow = R, byrow = TRUE)
  per_rep <- list()
  for (b in seq_len(reps)) {
    set.seed(seed_base + 404000L + b)
    x <- simulate_truth_guess(n, R, rho_vec, p_truth, guess_mat)
    boot <- bootstrap_fits(x, three_weight_kappas, boots, seed_offset = b)
    pt <- attr(boot, "point")
    ok <- apply(is.finite(boot), 2, all) & is.finite(pt)
    if (!all(ok)) next
    V <- cov(boot)
    cor_il <- V[1, 2] / sqrt(V[1, 1] * V[2, 2])
    cor_iq <- V[1, 3] / sqrt(V[1, 1] * V[3, 3])
    cor_lq <- V[2, 3] / sqrt(V[2, 2] * V[3, 3])
    per_rep[[length(per_rep) + 1L]] <- data.frame(
      b = b,
      kappa_id = pt[1], kappa_lin = pt[2], kappa_quad = pt[3],
      se_id = sqrt(V[1, 1]), se_lin = sqrt(V[2, 2]), se_quad = sqrt(V[3, 3]),
      cor_id_lin = cor_il, cor_id_quad = cor_iq, cor_lin_quad = cor_lq,
      delta_id_quad = pt[1] - pt[3],
      se_id_quad    = sqrt(max(V[1, 1] + V[3, 3] - 2 * V[1, 3], 0))
    )
  }
  rep_df <- do.call(rbind, per_rep)
  write.csv(rep_df, file.path(results_dir, "uc2_replicates.csv"), row.names = FALSE)
  summ <- data.frame(
    reps         = nrow(rep_df),
    mean_id      = mean(rep_df$kappa_id),
    mean_lin     = mean(rep_df$kappa_lin),
    mean_quad    = mean(rep_df$kappa_quad),
    mean_cor_id_lin = mean(rep_df$cor_id_lin),
    mean_cor_id_quad = mean(rep_df$cor_id_quad),
    mean_cor_lin_quad = mean(rep_df$cor_lin_quad),
    # Wald rate at 0.05 for identity vs quadratic (descriptive only - the
    # two have different population targets, so this is "are they distinct
    # enough that no single weight scheme is a good summary?").
    reject_id_eq_quad = mean(abs(rep_df$delta_id_quad / rep_df$se_id_quad) >
                             qnorm(0.975), na.rm = TRUE)
  )
  write.csv(summ, file.path(results_dir, "uc2_summary.csv"), row.names = FALSE)
  tic("UC2 done")
}

# ---- Dispatch -----------------------------------------------------------
if ("uc8" %in% ucs) run_uc8()
if ("uc1" %in% ucs) run_uc1()
if ("uc4" %in% ucs) run_uc4()
if ("uc2" %in% ucs) run_uc2()

# ---- Metadata -----------------------------------------------------------
meta <- data.frame(
  key = c("seed_base", "reps", "boots", "n", "C", "ucs",
          "R_version", "misskappa_version", "started_at", "elapsed_s"),
  value = c(
    as.character(seed_base),
    as.character(reps),
    as.character(boots),
    as.character(n),
    as.character(C),
    paste(ucs, collapse = ","),
    as.character(getRversion()),
    as.character(utils::packageVersion("misskappa")),
    format(t0, "%Y-%m-%dT%H:%M:%S"),
    sprintf("%.1f", as.numeric(difftime(Sys.time(), t0, units = "secs")))
  )
)
write.csv(meta, file.path(results_dir, "metadata.csv"), row.names = FALSE)

cat(sprintf("Wrote %s/{metadata,uc*_replicates,uc*_summary}.csv\n", results_dir))
