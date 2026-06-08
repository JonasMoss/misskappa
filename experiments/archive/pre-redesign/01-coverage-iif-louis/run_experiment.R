#!/usr/bin/env Rscript
#
# 01-coverage-iif-louis
#
# Empirical coverage of Wald confidence intervals for Conger's kappa
# under the IPW influence-function variance and the FIML Louis observed-
# information variance. Three DGPs (A / B / C) carried over from
# papers/combined/scripts/simulations_raw_three_estimators.R, three sample sizes,
# Conger x identity loss.
#
# Outputs (under results/):
#   coverage.csv     summary by (DGP, method, n): coverage rates,
#                    mean width, mean / SD of standardised residual
#   replicates.csv   per-replicate estimates and SEs (for plots)
#   metadata.csv     runtime arguments, package versions, design

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
      " --smoke         Fast smoke run (reps=5, single n=200). Default off.\n",
      " --reps N        Override replicate count.\n",
      " --info-rcond X  Override FIML Louis relative eigenvalue cutoff.\n",
      " --seed-base K   Seed base (default 1).\n",
      " --help, -h      This help.\n", sep = "")
  quit("no", status = 0)
}

smoke     <- has_flag("--smoke")
seed_base <- get_val("--seed-base", 1L, as.integer)
reps_user <- get_val("--reps", NA_integer_, as.integer)
info_rcond <- get_val("--info-rcond", NA_real_, as.numeric)

if (smoke) {
  ns   <- c(200L)
  reps <- 5L
} else {
  ns   <- c(500L, 2000L, 8000L)
  reps <- 200L
}
if (!is.na(reps_user)) reps <- reps_user

# ---- DGPs ---------------------------------------------------------------
C <- 5L
R <- 6L
p <- colSums(as.matrix(dat.fleiss1971)) / sum(as.matrix(dat.fleiss1971))
guess_mat <- matrix(rep(p, each = R), nrow = R, byrow = TRUE)

simulate_ratings <- function(n, rho_base, rho_truth_mult, guess) {
  truth <- sample.int(C, n, replace = TRUE, prob = p)
  x_star <- matrix(NA_integer_, nrow = n, ncol = R)
  for (j in seq_len(R)) {
    rho_i <- pmin(pmax(rho_base[j] * rho_truth_mult[truth], 0), 1)
    correct <- stats::rbinom(n, 1, prob = rho_i) == 1
    guessed <- sample.int(C, n, replace = TRUE, prob = guess[j, ])
    x_star[, j] <- ifelse(correct, truth, guessed)
  }
  list(x_star = x_star, truth = truth)
}

apply_missing_mcar <- function(x_star, pi_rater) {
  n <- nrow(x_star)
  M <- matrix(FALSE, nrow = n, ncol = R)
  for (j in seq_len(R)) M[, j] <- stats::runif(n) < pi_rater[j]
  x <- x_star
  x[!M] <- NA_integer_
  x
}

apply_missing_mar_truth <- function(x_star, truth, pi_truth) {
  n <- nrow(x_star)
  p_i <- pi_truth[truth]
  M <- matrix(stats::runif(n * R) < rep(p_i, times = R), nrow = n, ncol = R)
  x <- x_star
  x[!M] <- NA_integer_
  x
}

simulate_one <- function(dgp, n) {
  d <- simulate_ratings(n, dgp$rho_base, dgp$rho_truth_mult, dgp$guess)
  switch(dgp$missing,
         "mcar"      = apply_missing_mcar(d$x_star, dgp$pi_rater),
         "mar_truth" = apply_missing_mar_truth(d$x_star, d$truth, dgp$pi_truth),
         stop("unknown missing mechanism: ", dgp$missing))
}

dgpA <- list(label = "A",
             rho_base = rep(0.92, R), rho_truth_mult = rep(1, C),
             guess = guess_mat,
             missing = "mcar",
             pi_rater = c(0.95, 0.90, 0.80, 0.65, 0.50, 0.35))
dgpB <- list(label = "B",
             rho_base = c(0.98, 0.96, 0.92, 0.85, 0.72, 0.60),
             rho_truth_mult = rep(1, C),
             guess = guess_mat,
             missing = "mcar",
             pi_rater = c(0.95, 0.85, 0.70, 0.45, 0.25, 0.15))
dgpC <- list(label = "C",
             rho_base = rep(0.92, R),
             rho_truth_mult = c(1.05, 1.00, 0.95, 0.85, 0.70),
             guess = guess_mat,
             missing = "mar_truth",
             pi_truth = c(0.95, 0.90, 0.80, 0.55, 0.25))
dgps <- list(dgpA, dgpB, dgpC)

# ---- Truth: large complete-data sample, no missingness ------------------
n_truth <- if (smoke) 20000L else 200000L
estimate_truth <- function(dgp) {
  set.seed(seed_base + 999983L)  # truth seed independent of replicate seeds
  d <- simulate_ratings(n_truth, dgp$rho_base, dgp$rho_truth_mult, dgp$guess)
  k <- misskappa::kappa(d$x_star, method = "available", weight = "identity")
  as.numeric(k$estimates[["Conger"]])
}
truths <- vapply(dgps, estimate_truth, numeric(1))
names(truths) <- vapply(dgps, `[[`, character(1), "label")

# ---- Fitters ------------------------------------------------------------
fit_one <- function(x, method) {
  call_args <- list(x = x, method = method, weight = "identity")
  if (method == "fiml") {
    call_args$em_options <- list(max_iter = 50000L, tol = 1e-7)
    if (is.finite(info_rcond)) call_args$em_options$info_rcond <- info_rcond
  }
  res <- try(do.call(misskappa::kappa, call_args), silent = TRUE)
  if (inherits(res, "try-error")) return(c(est = NA_real_, se = NA_real_))
  est <- as.numeric(res$estimates[["Conger"]])
  v   <- res$vcov["Conger", "Conger"]
  se  <- if (is.finite(v) && v >= 0) sqrt(v) else NA_real_
  c(est = est, se = se)
}

# ---- Seeds (deterministic, non-colliding) -------------------------------
make_seed <- function(seed_base, dgp_idx, n_idx, b) {
  as.integer(seed_base) +
    as.integer(b) +
    1000L * (as.integer(n_idx) - 1L) +
    1000000L * (as.integer(dgp_idx) - 1L)
}

# ---- Run grid -----------------------------------------------------------
methods <- c("ipw", "fiml")
z975 <- stats::qnorm(0.975)
z95  <- stats::qnorm(0.95)
z995 <- stats::qnorm(0.995)

per_rep_rows <- list()
summary_rows <- list()

t0 <- Sys.time()
for (dgp_idx in seq_along(dgps)) {
  dgp <- dgps[[dgp_idx]]
  truth <- truths[[dgp$label]]
  for (n_idx in seq_along(ns)) {
    n <- ns[n_idx]
    cells <- vector("list", length(methods))
    names(cells) <- methods
    for (m in methods) {
      cells[[m]] <- matrix(NA_real_, nrow = reps, ncol = 2,
                           dimnames = list(NULL, c("est", "se")))
    }
    for (b in seq_len(reps)) {
      set.seed(make_seed(seed_base, dgp_idx, n_idx, b))
      x <- simulate_one(dgp, n)
      for (m in methods) {
        fr <- fit_one(x, m)
        cells[[m]][b, ] <- fr
        per_rep_rows[[length(per_rep_rows) + 1L]] <- data.frame(
          dgp = dgp$label, method = m, n = n, b = b, truth = truth,
          est = fr[["est"]], se = fr[["se"]]
        )
      }
    }
    for (m in methods) {
      est <- cells[[m]][, "est"]
      se  <- cells[[m]][, "se"]
      ok  <- is.finite(est) & is.finite(se) & se > 0
      z   <- (est[ok] - truth) / se[ok]
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(
        dgp           = dgp$label,
        method        = m,
        n             = n,
        B             = reps,
        truth         = truth,
        n_ok          = sum(ok),
        cov_90        = mean(abs(z) <= z95),
        cov_95        = mean(abs(z) <= z975),
        cov_99        = mean(abs(z) <= z995),
        mean_width_95 = 2 * z975 * mean(se[ok]),
        mean_z        = mean(z),
        sd_z          = stats::sd(z),
        mean_se       = mean(se[ok]),
        mc_sd_est     = stats::sd(est[ok])
      )
    }
    cat(sprintf("[%s] DGP %s, n=%d, reps=%d -- done (%.1fs)\n",
                format(Sys.time(), "%H:%M:%S"),
                dgp$label, n, reps,
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
}

coverage <- do.call(rbind, summary_rows)
replicates <- do.call(rbind, per_rep_rows)

# ---- Write --------------------------------------------------------------
results_dir <- "results"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(coverage,   file.path(results_dir, "coverage.csv"),   row.names = FALSE)
write.csv(replicates, file.path(results_dir, "replicates.csv"), row.names = FALSE)

meta <- data.frame(
  key = c("seed_base", "reps", "ns", "methods", "smoke",
          "C", "R", "n_truth", "weight", "kappa", "info_rcond",
          "R_version", "misskappa_version", "started_at", "elapsed_s"),
  value = c(
    as.character(seed_base),
    as.character(reps),
    paste(ns, collapse = ","),
    paste(methods, collapse = ","),
    as.character(smoke),
    as.character(C),
    as.character(R),
    as.character(n_truth),
    "identity",
    "Conger",
    if (is.finite(info_rcond)) as.character(info_rcond) else "default",
    as.character(getRversion()),
    as.character(utils::packageVersion("misskappa")),
    format(t0, "%Y-%m-%dT%H:%M:%S"),
    sprintf("%.1f", as.numeric(difftime(Sys.time(), t0, units = "secs")))
  )
)
write.csv(meta, file.path(results_dir, "metadata.csv"), row.names = FALSE)

cat(sprintf("Wrote %s/coverage.csv (%d rows)\n",   results_dir, nrow(coverage)))
cat(sprintf("Wrote %s/replicates.csv (%d rows)\n", results_dir, nrow(replicates)))
cat(sprintf("Wrote %s/metadata.csv\n",             results_dir))
