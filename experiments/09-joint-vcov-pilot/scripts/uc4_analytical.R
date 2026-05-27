#!/usr/bin/env Rscript
#
# UC4 analytical pass: rerun the Hausman test using the new IF-based
# joint_vcov() instead of the bootstrap from run_experiment.R. Same
# DGPs and seeds; produces results/uc4_analytical_summary.csv and
# results/uc4_analytical_replicates.csv.
#
# Companion to run_experiment.R; designed to be runnable after the new
# C++ field has landed (misskappa::joint_vcov and stats::influence on
# misskappa_estimate must work).

suppressPackageStartupMessages({
  library(misskappa)
})

args <- commandArgs(trailingOnly = TRUE)
has_flag <- function(name) name %in% args
get_val  <- function(name, default, parser = as.character) {
  i <- which(args == name)
  if (length(i) == 0L) return(default)
  parser(args[i + 1L])
}

smoke     <- has_flag("--smoke")
seed_base <- get_val("--seed-base", 1L, as.integer)
reps_user <- get_val("--reps", NA_integer_, as.integer)
n_user    <- get_val("--n", NA_integer_, as.integer)

reps <- if (smoke) 20L else 200L
n    <- 500L
if (!is.na(reps_user)) reps <- reps_user
if (!is.na(n_user))    n    <- n_user

C <- 5L
R <- 6L
results_dir <- "results"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# ---- DGP (matches run_experiment.R) -------------------------------------
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
t0 <- Sys.time()

for (cell_idx in seq_along(cells)) {
  cell_name <- names(cells)[cell_idx]
  cfg <- cells[[cell_idx]]
  for (b in seq_len(reps)) {
    # Same seed schedule as run_experiment.R::run_uc4 so the rows match
    # MC replicate by MC replicate.
    set.seed(seed_base + 303000L * cell_idx + b)
    x_star <- simulate_truth_guess(n, R, cfg$rho, p_truth, guess_mat)
    x      <- apply_mcar_independent(x_star, cfg$pi)

    ac  <- try(misskappa::kappa(x, method = "available", weight = "identity"),
               silent = TRUE)
    ipw <- try(misskappa::kappa(x, method = "ipw", weight = "identity"),
               silent = TRUE)
    if (inherits(ac, "try-error") || inherits(ipw, "try-error")) next

    V <- joint_vcov(ac = ac, ipw = ipw)
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
  cat(sprintf("[%s] cell=%s done (%.1fs)\n",
              format(Sys.time(), "%H:%M:%S"), cell_name,
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

rep_df <- do.call(rbind, per_rep)
write.csv(rep_df,
          file.path(results_dir, "uc4_analytical_replicates.csv"),
          row.names = FALSE)

summ <- do.call(rbind, by(rep_df, rep_df$cell, function(d) {
  data.frame(
    cell         = unique(d$cell),
    reps         = nrow(d),
    mean_AC      = mean(d$kappa_AC),
    mean_IPW     = mean(d$kappa_IPW),
    mean_delta   = mean(d$delta),
    sd_delta_mc  = stats::sd(d$delta),
    mean_se_analytical = mean(d$se, na.rm = TRUE),
    mean_cor_AC_IPW = mean(d$cov_AC_IPW /
                             sqrt(pmax(d$var_AC * d$var_IPW, 0)),
                           na.rm = TRUE),
    reject_05    = mean(d$reject_05, na.rm = TRUE)
  )
}))
rownames(summ) <- NULL
write.csv(summ,
          file.path(results_dir, "uc4_analytical_summary.csv"),
          row.names = FALSE)
cat(sprintf("Wrote %s/uc4_analytical_*.csv\n", results_dir))
