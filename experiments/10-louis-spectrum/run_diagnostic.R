#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(misskappa)
})

args <- commandArgs(trailingOnly = TRUE)
has_flag <- function(name) name %in% args
get_val <- function(name, default, parser = as.character) {
  i <- which(args == name)
  if (length(i) == 0L) return(default)
  parser(args[i + 1L])
}

if (has_flag("--help") || has_flag("-h")) {
  cat("Usage: Rscript run_diagnostic.R [options]\n",
      " --n N            Sample size (default 8000).\n",
      " --seed K         RNG seed (default 20260527).\n",
      " --info-rcond X   Rank cutoff to evaluate (default 5e-5).\n",
      " --help, -h       This help.\n", sep = "")
  quit("no", status = 0)
}

n <- get_val("--n", 8000L, as.integer)
seed <- get_val("--seed", 20260527L, as.integer)
info_rcond <- get_val("--info-rcond", 5e-5, as.numeric)

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

dgpA <- list(label = "A",
             rho_base = rep(0.92, R),
             rho_truth_mult = rep(1, C),
             guess = guess_mat,
             pi_rater = c(0.95, 0.90, 0.80, 0.65, 0.50, 0.35))

set.seed(seed)
d <- simulate_ratings(n, dgpA$rho_base, dgpA$rho_truth_mult, dgpA$guess)
x <- apply_missing_mcar(d$x_star, dgpA$pi_rater)

base_opts <- list(max_iter = 50000L, tol = 1e-7)
diag_all <- misskappa:::rcpp_fiml_louis_spectrum(
  x, weight_type = "identity", values = NULL,
  em_options = c(base_opts, list(info_rcond = 0))
)
diag_cut <- misskappa:::rcpp_fiml_louis_spectrum(
  x, weight_type = "identity", values = NULL,
  em_options = c(base_opts, list(info_rcond = info_rcond))
)
fit_cut <- misskappa::kappa(
  x, method = "fiml", weight = "identity",
  em_options = c(base_opts, list(info_rcond = info_rcond))
)

eigenvalues <- diag_all$eigenvalues
relative <- eigenvalues / diag_all$lambda_max
retained <- eigenvalues > diag_cut$threshold
contrib_all <- diag_all$variance_contribution
contrib_cut <- ifelse(retained, contrib_all, 0)
total_all <- sum(contrib_all)
total_cut <- sum(contrib_cut)

out_dir <- "results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

singular_values <- data.frame(
  index = seq_along(eigenvalues),
  eigenvalue = eigenvalues,
  relative = relative,
  retained_at_info_rcond = retained
)
write.csv(singular_values, file.path(out_dir, "singular_values.csv"), row.names = FALSE)

gradient_decomposition <- data.frame(
  index = seq_along(eigenvalues),
  eigenvalue = eigenvalues,
  relative = relative,
  gradient_projection = diag_all$gradient_projection
)
write.csv(gradient_decomposition,
          file.path(out_dir, "gradient_decomposition.csv"),
          row.names = FALSE)

variance_contributions <- data.frame(
  index = seq_along(eigenvalues),
  eigenvalue = eigenvalues,
  relative = relative,
  contribution_untruncated = contrib_all,
  contribution_truncated = contrib_cut,
  share_untruncated = contrib_all / total_all,
  cumulative_share_untruncated = cumsum(contrib_all) / total_all,
  retained_at_info_rcond = retained
)
write.csv(variance_contributions,
          file.path(out_dir, "variance_contributions.csv"),
          row.names = FALSE)

share_below_1e3 <- sum(contrib_all[relative <= 1e-3]) / total_all
share_dropped <- sum(contrib_all[!retained]) / total_all

summary_lines <- c(
  "# Louis Spectrum Diagnostic",
  "",
  sprintf("- DGP: A, homogeneous raters, MCAR"),
  sprintf("- n: %d", n),
  sprintf("- seed: %d", seed),
  sprintf("- kappa Conger: %.10f", diag_cut$kappa_conger),
  sprintf("- package vcov Conger: %.10g", fit_cut$vcov["Conger", "Conger"]),
  sprintf("- diagnostic truncated variance: %.10g", total_cut),
  sprintf("- untruncated variance: %.10g", total_all),
  sprintf("- lambda_max: %.10g", diag_all$lambda_max),
  sprintf("- evaluated info_rcond: %.3g", info_rcond),
  sprintf("- cutoff threshold: %.10g", diag_cut$threshold),
  sprintf("- retained rank: %d / %d", diag_cut$retained_rank, length(eigenvalues)),
  sprintf("- pruned theta support size: %d", diag_cut$n_patterns),
  sprintf("- share of untruncated variance from relative eigenvalues <= 1e-3: %.4f",
          share_below_1e3),
  sprintf("- share of untruncated variance dropped by info_rcond: %.4f",
          share_dropped)
)
writeLines(summary_lines, file.path(out_dir, "summary.md"))

cat(paste(summary_lines, collapse = "\n"), "\n")
