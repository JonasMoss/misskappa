#!/usr/bin/env Rscript
#
# 09-joint-vcov-pilot
#
# Joint inference across non-independent kappa estimates on the same
# data ("Case 3" of the equality-of-two-kappas question). All four
# sub-experiments use the analytical `misskappa::joint_vcov()` helper,
# which assembles the joint asymptotic covariance from the per-subject
# influence functions exposed by the categorical raw estimators
# (`available`, `ipw`, `gwet`). No bootstrap.
#
# Sub-experiments:
#   UC8 - Fleiss vs BP same-fit contrast on a single available-case fit.
#         Uses the existing 3x3 vcov; included for completeness.
#   UC1 - All choose(R, 2) pairwise Cohen kappas + their joint Wald
#         test of homogeneity. Size under exchangeable raters, power
#         under one miscalibrated rater.
#   UC4 - Hausman test of kappa_AC = kappa_IPW on the same incomplete
#         data. Size under exchangeable + MCAR, power under the
#         experiment-03 non-exchangeable cell.
#   UC2 - Joint inference across three weighted kappas (identity,
#         linear, quadratic) on the same complete data.
#
# Outputs (under results/):
#   uc8_replicates.csv   per-replicate Fleiss/BP estimates + contrast z
#   uc8_summary.csv      cell-level rejection rate at 0.05
#   uc1_replicates.csv   per-replicate pairwise kappas + Wald stat
#   uc1_summary.csv      cell-level size / power at 0.05
#   uc4_replicates.csv   per-replicate AC + IPW + Wald z for contrast
#   uc4_summary.csv      cell-level size / power at 0.05
#   uc2_replicates.csv   per-replicate (identity, linear, quadratic)
#                        kappas with joint vcov entries
#   uc2_summary.csv      cell-level: means, SDs, correlations across
#                        the three weight schemes
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
      " --smoke         Fast smoke run (reps=20).\n",
      " --reps N        MC replicates per cell (default 200).\n",
      " --n N           Subjects per replicate (default 500).\n",
      " --seed-base K   Seed base (default 1).\n",
      " --only UC       Run only one UC: uc8|uc1|uc4|uc2 (default all).\n",
      " --help, -h      This help.\n", sep = "")
  quit("no", status = 0)
}

smoke      <- has_flag("--smoke")
seed_base  <- get_val("--seed-base", 1L, as.integer)
reps_user  <- get_val("--reps", NA_integer_, as.integer)
n_user     <- get_val("--n", NA_integer_, as.integer)
only_user  <- get_val("--only", NA_character_, as.character)

reps  <- if (smoke) 20L else 200L
n     <- 500L
if (!is.na(reps_user))  reps  <- reps_user
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

# Wald quadratic form with SVD-based pseudoinverse so small singular
# values don't blow up the test statistic. Returns (stat, df).
wald_quadform <- function(theta, V, A = diag(length(theta)),
                          rcond_tol = 1e-8) {
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

tic <- function(tag) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), tag))

# ===================== UC8: Fleiss vs BP free contrast ===================
#
# Same fit, two coefficients with known covariance. delta = kappa_F - BP
# is zero in population iff the rater marginals are uniform (under
# identity loss); under skewed marginals the contrast is non-zero and
# the test has power.

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
        reject_05 = if (is.na(z)) NA else abs(z) > z975
      )
    }
  }
  rep_df <- do.call(rbind, per_rep)
  write.csv(rep_df, file.path(results_dir, "uc8_replicates.csv"),
            row.names = FALSE)

  summ <- do.call(rbind, by(rep_df, rep_df$marg, function(d) {
    data.frame(
      marg          = unique(d$marg),
      reps          = nrow(d),
      mean_kappa_F  = mean(d$kappa_F),
      mean_kappa_BP = mean(d$kappa_BP),
      mean_delta    = mean(d$delta),
      sd_delta      = sd(d$delta),
      mean_se       = mean(d$se),
      reject_05     = mean(d$reject_05, na.rm = TRUE)
    )
  }))
  rownames(summ) <- NULL
  write.csv(summ, file.path(results_dir, "uc8_summary.csv"),
            row.names = FALSE)
  tic("UC8 done")
}

# ===================== UC1: pairwise Cohen joint vcov ====================
#
# All choose(R, 2) pairwise Cohen kappas with analytical joint vcov
# from the per-subject Conger influence functions on each 2-rater
# submatrix. Wald test of pairwise homogeneity uses contrasts
# kappa_p - kappa_1 for p = 2..K (chi-sq df=14).

run_uc1 <- function() {
  tic("UC1 start")
  R <- 6L
  p_truth   <- c(0.05, 0.10, 0.20, 0.30, 0.35)
  guess_mat <- matrix(rep(p_truth, each = R), nrow = R, byrow = TRUE)

  cells <- list(
    null    = rep(0.92, R),
    one_bad = c(0.92, 0.92, 0.92, 0.92, 0.92, 0.60)
  )
  pairs <- combn(R, 2)
  n_pairs <- ncol(pairs)
  pair_names <- apply(pairs, 2,
                      function(rs) sprintf("k_%d_%d", rs[1], rs[2]))
  A <- cbind(rep(-1, n_pairs - 1L), diag(n_pairs - 1L))

  per_rep <- list()
  for (cell_idx in seq_along(cells)) {
    cell_name <- names(cells)[cell_idx]
    rho_vec <- cells[[cell_idx]]
    for (b in seq_len(reps)) {
      set.seed(seed_base + 202000L * cell_idx + b)
      x <- simulate_truth_guess(n, R, rho_vec, p_truth, guess_mat)

      # Fit each pair; collect Conger point estimate and IF column.
      point <- numeric(n_pairs)
      psi_cols <- vector("list", n_pairs)
      failed <- FALSE
      for (p_idx in seq_len(n_pairs)) {
        r <- pairs[1, p_idx]; s <- pairs[2, p_idx]
        fit <- try(misskappa::kappa(x[, c(r, s)], method = "available",
                                    weight = "identity"),
                   silent = TRUE)
        if (inherits(fit, "try-error")) { failed <- TRUE; break }
        point[p_idx] <- as.numeric(coef(fit)[["Conger"]])
        psi_cols[[p_idx]] <- fit$psi[, "Conger"]
      }
      if (failed || !all(is.finite(point))) next

      psi_stack <- do.call(cbind, psi_cols)
      V <- crossprod(psi_stack) / (as.numeric(n) * as.numeric(n))
      w <- wald_quadform(point, V, A)
      p_chi <- pchisq(w$stat, df = w$df, lower.tail = FALSE)

      row <- as.data.frame(t(setNames(point, pair_names)))
      row$cell      <- cell_name
      row$b         <- b
      row$wald_stat <- w$stat
      row$wald_df   <- w$df
      row$p_chi     <- p_chi
      row$reject_05 <- p_chi < 0.05
      per_rep[[length(per_rep) + 1L]] <- row
    }
  }
  rep_df <- do.call(rbind, per_rep)
  write.csv(rep_df, file.path(results_dir, "uc1_replicates.csv"),
            row.names = FALSE)

  summ <- do.call(rbind, by(rep_df, rep_df$cell, function(d) {
    pair_cols <- grep("^k_", colnames(d), value = TRUE)
    means <- sapply(pair_cols, function(cn) mean(d[[cn]], na.rm = TRUE))
    data.frame(
      cell           = unique(d$cell),
      reps           = nrow(d),
      mean_min_kappa = min(means),
      mean_max_kappa = max(means),
      mean_range     = max(means) - min(means),
      mean_wald      = mean(d$wald_stat, na.rm = TRUE),
      median_wald    = median(d$wald_stat, na.rm = TRUE),
      reject_05      = mean(d$reject_05, na.rm = TRUE)
    )
  }))
  rownames(summ) <- NULL
  write.csv(summ, file.path(results_dir, "uc1_summary.csv"),
            row.names = FALSE)
  tic("UC1 done")
}

# ===================== UC4: Hausman AC vs IPW ============================
#
# joint_vcov(ac, ipw) returns a 6x6 block matrix; the Conger contrast
# c = (1, 0, 0, -1, 0, 0) gives a chi-sq df=1 Wald test of
# kappa_AC = kappa_IPW.

run_uc4 <- function() {
  tic("UC4 start")
  R <- 6L
  p_truth   <- c(0.05, 0.10, 0.20, 0.30, 0.35)
  guess_mat <- matrix(rep(p_truth, each = R), nrow = R, byrow = TRUE)

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

      ac  <- try(misskappa::kappa(x, method = "available",
                                  weight = "identity"),
                 silent = TRUE)
      ipw <- try(misskappa::kappa(x, method = "ipw",
                                  weight = "identity"),
                 silent = TRUE)
      if (inherits(ac, "try-error") || inherits(ipw, "try-error")) next

      V <- misskappa::joint_vcov(ac = ac, ipw = ipw)
      delta <- as.numeric(coef(ac)[["Conger"]] - coef(ipw)[["Conger"]])
      v <- V["ac.Conger", "ac.Conger"] +
           V["ipw.Conger", "ipw.Conger"] -
           2 * V["ac.Conger", "ipw.Conger"]
      se <- sqrt(max(v, 0))
      z  <- if (se > 0) delta / se else NA_real_

      per_rep[[length(per_rep) + 1L]] <- data.frame(
        cell = cell_name, b = b,
        kappa_AC = as.numeric(coef(ac)[["Conger"]]),
        kappa_IPW = as.numeric(coef(ipw)[["Conger"]]),
        delta = delta, se = se, z = z,
        var_AC  = V["ac.Conger", "ac.Conger"],
        var_IPW = V["ipw.Conger", "ipw.Conger"],
        cov_AC_IPW = V["ac.Conger", "ipw.Conger"],
        reject_05 = if (is.na(z)) NA else abs(z) > z975
      )
    }
  }
  rep_df <- do.call(rbind, per_rep)
  write.csv(rep_df, file.path(results_dir, "uc4_replicates.csv"),
            row.names = FALSE)

  summ <- do.call(rbind, by(rep_df, rep_df$cell, function(d) {
    data.frame(
      cell        = unique(d$cell),
      reps        = nrow(d),
      mean_AC     = mean(d$kappa_AC),
      mean_IPW    = mean(d$kappa_IPW),
      mean_delta  = mean(d$delta),
      sd_delta_mc = sd(d$delta),
      mean_se     = mean(d$se, na.rm = TRUE),
      mean_cor_AC_IPW = mean(d$cov_AC_IPW /
                               sqrt(pmax(d$var_AC * d$var_IPW, 0)),
                             na.rm = TRUE),
      reject_05   = mean(d$reject_05, na.rm = TRUE)
    )
  }))
  rownames(summ) <- NULL
  write.csv(summ, file.path(results_dir, "uc4_summary.csv"),
            row.names = FALSE)
  tic("UC4 done")
}

# ===================== UC2: weight-scheme sensitivity ====================
#
# Three weighted kappas on the same complete data. joint_vcov() returns
# the 9x9 block matrix over (Conger, Fleiss, BP) x (identity, linear,
# quadratic); we report the 3x3 Conger sub-block plus a one-df Wald
# test of identity == quadratic.

run_uc2 <- function() {
  tic("UC2 start")
  R <- 6L
  p_truth   <- c(0.05, 0.10, 0.20, 0.30, 0.35)
  rho_vec   <- rep(0.85, R)
  guess_mat <- matrix(rep(p_truth, each = R), nrow = R, byrow = TRUE)

  z975 <- qnorm(0.975)
  per_rep <- list()
  for (b in seq_len(reps)) {
    set.seed(seed_base + 404000L + b)
    x <- simulate_truth_guess(n, R, rho_vec, p_truth, guess_mat)

    k_id <- try(misskappa::kappa(x, method = "available",
                                 weight = "identity"), silent = TRUE)
    k_li <- try(misskappa::kappa(x, method = "available",
                                 weight = "linear"), silent = TRUE)
    k_qu <- try(misskappa::kappa(x, method = "available",
                                 weight = "quadratic"), silent = TRUE)
    if (inherits(k_id, "try-error") ||
        inherits(k_li, "try-error") ||
        inherits(k_qu, "try-error")) next

    V <- misskappa::joint_vcov(id = k_id, lin = k_li, quad = k_qu)
    keys <- c("id.Conger", "lin.Conger", "quad.Conger")
    V_C <- V[keys, keys]
    pt <- c(coef(k_id)[["Conger"]], coef(k_li)[["Conger"]],
            coef(k_qu)[["Conger"]])
    cor_il <- V_C[1, 2] / sqrt(V_C[1, 1] * V_C[2, 2])
    cor_iq <- V_C[1, 3] / sqrt(V_C[1, 1] * V_C[3, 3])
    cor_lq <- V_C[2, 3] / sqrt(V_C[2, 2] * V_C[3, 3])
    delta_iq <- pt[1] - pt[3]
    var_iq <- V_C[1, 1] + V_C[3, 3] - 2 * V_C[1, 3]
    se_iq <- sqrt(max(var_iq, 0))
    z_iq <- if (se_iq > 0) delta_iq / se_iq else NA_real_

    per_rep[[length(per_rep) + 1L]] <- data.frame(
      b = b,
      kappa_id = pt[1], kappa_lin = pt[2], kappa_quad = pt[3],
      se_id = sqrt(V_C[1, 1]), se_lin = sqrt(V_C[2, 2]),
      se_quad = sqrt(V_C[3, 3]),
      cor_id_lin   = cor_il,
      cor_id_quad  = cor_iq,
      cor_lin_quad = cor_lq,
      delta_id_quad = delta_iq,
      se_id_quad    = se_iq,
      reject_id_eq_quad = if (is.na(z_iq)) NA else abs(z_iq) > z975
    )
  }
  rep_df <- do.call(rbind, per_rep)
  write.csv(rep_df, file.path(results_dir, "uc2_replicates.csv"),
            row.names = FALSE)
  summ <- data.frame(
    reps              = nrow(rep_df),
    mean_id           = mean(rep_df$kappa_id),
    mean_lin          = mean(rep_df$kappa_lin),
    mean_quad         = mean(rep_df$kappa_quad),
    mean_cor_id_lin   = mean(rep_df$cor_id_lin),
    mean_cor_id_quad  = mean(rep_df$cor_id_quad),
    mean_cor_lin_quad = mean(rep_df$cor_lin_quad),
    reject_id_eq_quad = mean(rep_df$reject_id_eq_quad, na.rm = TRUE)
  )
  write.csv(summ, file.path(results_dir, "uc2_summary.csv"),
            row.names = FALSE)
  tic("UC2 done")
}

# ---- Dispatch -----------------------------------------------------------
if ("uc8" %in% ucs) run_uc8()
if ("uc1" %in% ucs) run_uc1()
if ("uc4" %in% ucs) run_uc4()
if ("uc2" %in% ucs) run_uc2()

# ---- Metadata -----------------------------------------------------------
meta <- data.frame(
  key = c("seed_base", "reps", "n", "C", "ucs",
          "R_version", "misskappa_version", "started_at", "elapsed_s"),
  value = c(
    as.character(seed_base),
    as.character(reps),
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

cat(sprintf("Wrote %s/{metadata,uc*_replicates,uc*_summary}.csv\n",
            results_dir))
