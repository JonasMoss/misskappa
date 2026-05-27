#!/usr/bin/env Rscript
#
# 03-ac-vs-ipw-efficiency
#
# When does IPW's MSE beat AC's for Conger's kappa? Factorial sweep over
#   (i)   rater exchangeability        (2 levels: exch, nonexch)
#   (ii)  variance of pi_j over raters (3 levels: zero, mid, high)
#   (iii) within-pair correlation of M (3 levels: 0, 0.4, 0.8)
# at one moderate n. Reports bias, MC SD, MSE for AC and IPW.
#
# Rater model: latent truth + guess (matches paper/scripts). Missingness:
# Gaussian-copula latent with equicorrelated within-subject structure. The
# mechanism is independent of the latent truth, so IPW is consistent in
# every cell; AC is consistent only under exchangeable raters.
#
# Outputs (under results/):
#   mse.csv         cell-level: bias, sd, mse for AC and IPW
#   replicates.csv  per-replicate point estimates (long form)
#   metadata.csv    run metadata

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
      " --smoke         Fast smoke run (reps=10, n_truth=20000).\n",
      " --reps N        Override replicate count (default 500).\n",
      " --n N           Sample size per replicate (default 1000).\n",
      " --seed-base K   Seed base (default 1).\n",
      " --help, -h      This help.\n", sep = "")
  quit("no", status = 0)
}

smoke     <- has_flag("--smoke")
seed_base <- get_val("--seed-base", 1L, as.integer)
reps_user <- get_val("--reps", NA_integer_, as.integer)
n_user    <- get_val("--n", NA_integer_, as.integer)

reps <- if (smoke) 10L else 500L
n    <- 1000L
if (!is.na(reps_user)) reps <- reps_user
if (!is.na(n_user))    n    <- n_user

# ---- Rating model -------------------------------------------------------
C <- 5L
R <- 6L
p <- colSums(as.matrix(dat.fleiss1971)) / sum(as.matrix(dat.fleiss1971))
guess_mat <- matrix(rep(p, each = R), nrow = R, byrow = TRUE)

simulate_ratings <- function(n, rho_base, guess) {
  truth <- sample.int(C, n, replace = TRUE, prob = p)
  x_star <- matrix(NA_integer_, nrow = n, ncol = R)
  for (j in seq_len(R)) {
    correct <- stats::rbinom(n, 1L, prob = rho_base[j]) == 1L
    guessed <- sample.int(C, n, replace = TRUE, prob = guess[j, ])
    x_star[, j] <- ifelse(correct, truth, guessed)
  }
  x_star
}

# Equicorrelated Gaussian copula. Latent L_{i,j} = sqrt(rho_M) Z_i +
# sqrt(1-rho_M) V_{i,j}, with Z_i, V_{i,j} ~ N(0,1) iid. Set
# M_{i,j} = 1 iff L_{i,j} < Phi^{-1}(pi_j). Marginal P(M_{i,j} = 1) = pi_j;
# tetrachoric correlation across raters within subject = rho_M.
apply_missing_corr <- function(x_star, pi_rater, rho_M = 0) {
  n <- nrow(x_star); R <- ncol(x_star)
  if (rho_M < 0 || rho_M >= 1) stop("rho_M must be in [0, 1).")
  Z <- if (rho_M > 0) stats::rnorm(n) else rep(0, n)
  V <- matrix(stats::rnorm(n * R), nrow = n, ncol = R)
  L <- sqrt(rho_M) * Z + sqrt(1 - rho_M) * V
  cutoff <- matrix(stats::qnorm(pi_rater), nrow = n, ncol = R, byrow = TRUE)
  M <- L < cutoff
  x <- x_star
  x[!M] <- NA_integer_
  x
}

# ---- Levels -------------------------------------------------------------
rho_levels <- list(
  exch    = rep(0.92, R),
  nonexch = c(0.98, 0.96, 0.92, 0.85, 0.72, 0.60)
)

# All pi schedules share mean 0.6, so the expected observed cell count is
# constant across pi_var levels.
pi_levels <- list(
  zero = rep(0.6, R),
  mid  = 0.6 + c(0.20, 0.15, 0.10, -0.10, -0.15, -0.20),
  high = 0.6 + c(0.35, 0.25, 0.15, -0.15, -0.25, -0.35)
)

corr_M_levels <- c(none = 0.0, mod = 0.4, high = 0.8)

# ---- Truth: huge complete-data sample, no missingness ------------------
n_truth <- if (smoke) 20000L else 200000L
estimate_truth <- function(rho_base) {
  set.seed(seed_base + 999983L)
  x_star <- simulate_ratings(n_truth, rho_base, guess_mat)
  k <- misskappa::kappa(x_star, method = "available", weight = "identity")
  as.numeric(k$estimates[["Conger"]])
}
truths <- vapply(rho_levels, estimate_truth, numeric(1))

# ---- Fitter -------------------------------------------------------------
fit_one <- function(x, method) {
  res <- try(misskappa::kappa(x, method = method, weight = "identity"),
             silent = TRUE)
  if (inherits(res, "try-error")) return(NA_real_)
  as.numeric(res$estimates[["Conger"]])
}

# ---- Seeds --------------------------------------------------------------
make_seed <- function(seed_base, ex_idx, pi_idx, rho_idx, b) {
  as.integer(seed_base) +
    as.integer(b) +
    1000L  * (as.integer(rho_idx) - 1L) +
    10000L * (as.integer(pi_idx)  - 1L) +
    100000L * (as.integer(ex_idx) - 1L)
}

# ---- Sweep --------------------------------------------------------------
methods     <- c("available", "ipw")
ex_names    <- names(rho_levels)
pi_names    <- names(pi_levels)
rho_names_M <- names(corr_M_levels)

per_rep_rows <- list()
summary_rows <- list()

t0 <- Sys.time()
for (ex_idx in seq_along(rho_levels)) {
  ex_name  <- ex_names[ex_idx]
  rho_base <- rho_levels[[ex_idx]]
  truth    <- truths[[ex_name]]
  for (pi_idx in seq_along(pi_levels)) {
    pi_name  <- pi_names[pi_idx]
    pi_rater <- pi_levels[[pi_idx]]
    pi_var_pop <- mean((pi_rater - mean(pi_rater))^2)
    for (rho_idx in seq_along(corr_M_levels)) {
      rho_name <- rho_names_M[rho_idx]
      rho_M    <- corr_M_levels[rho_idx]
      cells <- matrix(NA_real_, nrow = reps, ncol = length(methods),
                      dimnames = list(NULL, methods))
      for (b in seq_len(reps)) {
        set.seed(make_seed(seed_base, ex_idx, pi_idx, rho_idx, b))
        x_star <- simulate_ratings(n, rho_base, guess_mat)
        x      <- apply_missing_corr(x_star, pi_rater, rho_M = rho_M)
        for (m in methods) {
          cells[b, m] <- fit_one(x, m)
        }
      }
      for (m in methods) {
        est  <- cells[, m]
        ok   <- is.finite(est)
        n_ok <- sum(ok)
        bias <- mean(est[ok]) - truth
        sd_e <- stats::sd(est[ok])
        mse  <- bias^2 + sd_e^2
        summary_rows[[length(summary_rows) + 1L]] <- data.frame(
          exch_label         = ex_name,
          pi_var_label       = pi_name,
          corr_M_label       = rho_name,
          method             = m,
          n                  = n,
          B                  = reps,
          truth              = truth,
          pi_var_pop         = pi_var_pop,
          tetrachoric_corr_M = rho_M,
          n_ok               = n_ok,
          mean_est           = mean(est[ok]),
          bias               = bias,
          sd_est             = sd_e,
          mse                = mse,
          mc_se_bias         = sd_e / sqrt(n_ok),
          mc_se_sd           = sd_e / sqrt(2 * n_ok)
        )
        for (b in seq_len(reps)) {
          per_rep_rows[[length(per_rep_rows) + 1L]] <- data.frame(
            exch_label   = ex_name,
            pi_var_label = pi_name,
            corr_M_label = rho_name,
            method       = m,
            n            = n,
            b            = b,
            truth        = truth,
            est          = est[b]
          )
        }
      }
      cat(sprintf("[%s] exch=%s pi_var=%s corr_M=%s -- done (%.1fs)\n",
                  format(Sys.time(), "%H:%M:%S"),
                  ex_name, pi_name, rho_name,
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    }
  }
}

mse_table  <- do.call(rbind, summary_rows)
replicates <- do.call(rbind, per_rep_rows)

# ---- Write --------------------------------------------------------------
results_dir <- "results"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(mse_table,  file.path(results_dir, "mse.csv"),        row.names = FALSE)
write.csv(replicates, file.path(results_dir, "replicates.csv"), row.names = FALSE)

meta <- data.frame(
  key = c("seed_base", "reps", "n", "methods", "smoke",
          "C", "R", "n_truth", "weight", "kappa",
          "exch_levels", "pi_levels", "corr_M_levels",
          "R_version", "misskappa_version", "started_at", "elapsed_s"),
  value = c(
    as.character(seed_base),
    as.character(reps),
    as.character(n),
    paste(methods, collapse = ","),
    as.character(smoke),
    as.character(C),
    as.character(R),
    as.character(n_truth),
    "identity",
    "Conger",
    paste(ex_names, collapse = ","),
    paste(pi_names, collapse = ","),
    paste(rho_names_M, collapse = ","),
    as.character(getRversion()),
    as.character(utils::packageVersion("misskappa")),
    format(t0, "%Y-%m-%dT%H:%M:%S"),
    sprintf("%.1f", as.numeric(difftime(Sys.time(), t0, units = "secs")))
  )
)
write.csv(meta, file.path(results_dir, "metadata.csv"), row.names = FALSE)

cat(sprintf("Wrote %s/mse.csv (%d rows)\n",        results_dir, nrow(mse_table)))
cat(sprintf("Wrote %s/replicates.csv (%d rows)\n", results_dir, nrow(replicates)))
cat(sprintf("Wrote %s/metadata.csv\n",             results_dir))
