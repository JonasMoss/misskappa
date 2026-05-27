#!/usr/bin/env Rscript
#
# UC1 analytical pass: rerun the pairwise Cohen homogeneity test using
# the new IF-based joint_vcov() instead of bootstrap. Same DGPs and
# seeds as run_uc1() in run_experiment.R. Tests whether the bootstrap
# size inflation in the original 15-dim bootstrap collapses to nominal
# 0.05 once the analytical joint vcov is used.
#
# Outputs results/uc1_analytical_{summary,replicates}.csv.

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

p_truth   <- c(0.05, 0.10, 0.20, 0.30, 0.35)
guess_mat <- matrix(rep(p_truth, each = R), nrow = R, byrow = TRUE)

cells <- list(
  null    = rep(0.92, R),
  one_bad = c(0.92, 0.92, 0.92, 0.92, 0.92, 0.60)
)

pairs <- combn(R, 2)
n_pairs <- ncol(pairs)
pair_names <- apply(pairs, 2, function(rs) sprintf("p%d_%d", rs[1], rs[2]))

# Contrast matrix: zero out kappa_1; rows give kappa_p - kappa_1.
A <- cbind(rep(-1, n_pairs - 1L), diag(n_pairs - 1L))

per_rep <- list()
t0 <- Sys.time()

for (cell_idx in seq_along(cells)) {
  cell_name <- names(cells)[cell_idx]
  rho_vec <- cells[[cell_idx]]
  for (b in seq_len(reps)) {
    # Same seed schedule as run_experiment.R::run_uc1.
    set.seed(seed_base + 202000L * cell_idx + b)
    x <- simulate_truth_guess(n, R, rho_vec, p_truth, guess_mat)

    # Per-pair fit; the per-subject IF on Conger for R=2 IS the pairwise
    # Cohen IF (Conger collapses to Cohen for R=2).
    fits <- vector("list", n_pairs)
    names(fits) <- pair_names
    point <- numeric(n_pairs)
    failed <- FALSE
    for (p_idx in seq_len(n_pairs)) {
      r <- pairs[1, p_idx]; s <- pairs[2, p_idx]
      fit <- try(misskappa::kappa(x[, c(r, s)], method = "available",
                                  weight = "identity"),
                 silent = TRUE)
      if (inherits(fit, "try-error")) { failed <- TRUE; break }
      fits[[p_idx]] <- fit
      point[p_idx] <- as.numeric(coef(fit)[["Conger"]])
    }
    if (failed || !all(is.finite(point))) next

    # Stack per-subject Conger IFs into an n x n_pairs matrix; joint vcov
    # is (1 / n^2) * crossprod. We can use joint_vcov() but it returns a
    # 3*n_pairs x 3*n_pairs block over (Conger, Fleiss, BP); easier to
    # build the n x n_pairs Conger-only stack directly here.
    psi_stack <- do.call(cbind, lapply(fits, function(f) {
      stats::influence(f)[, "Conger", drop = FALSE]
    }))
    V <- crossprod(psi_stack) / (as.numeric(n) * as.numeric(n))

    delta <- A %*% point
    M <- A %*% V %*% t(A)
    # Wald with pseudoinverse for safety (M is 14 x 14, generally full rank
    # but small singular values still possible).
    sv <- svd(M)
    keep <- sv$d > sv$d[1] * 1e-8
    M_pinv <- sv$v[, keep, drop = FALSE] %*%
              diag(1 / sv$d[keep], nrow = sum(keep)) %*%
              t(sv$u[, keep, drop = FALSE])
    stat <- as.numeric(t(delta) %*% M_pinv %*% delta)
    df <- sum(keep)
    p_chi <- pchisq(stat, df = df, lower.tail = FALSE)

    row <- as.data.frame(t(setNames(point, pair_names)))
    row$cell <- cell_name
    row$b <- b
    row$wald_stat <- stat
    row$wald_df <- df
    row$p_chi <- p_chi
    row$reject_05 <- p_chi < 0.05
    per_rep[[length(per_rep) + 1L]] <- row
  }
  cat(sprintf("[%s] cell=%s done (%.1fs)\n",
              format(Sys.time(), "%H:%M:%S"), cell_name,
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

rep_df <- do.call(rbind, per_rep)
write.csv(rep_df,
          file.path(results_dir, "uc1_analytical_replicates.csv"),
          row.names = FALSE)

summ <- do.call(rbind, by(rep_df, rep_df$cell, function(d) {
  pair_cols <- grep("^p[0-9]", colnames(d), value = TRUE)
  means <- sapply(pair_cols, function(cn) mean(d[[cn]], na.rm = TRUE))
  data.frame(
    cell           = unique(d$cell),
    reps           = nrow(d),
    mean_min_kappa = min(means),
    mean_max_kappa = max(means),
    mean_range     = max(means) - min(means),
    mean_wald      = mean(d$wald_stat, na.rm = TRUE),
    median_wald    = stats::median(d$wald_stat, na.rm = TRUE),
    reject_05      = mean(d$reject_05, na.rm = TRUE)
  )
}))
rownames(summ) <- NULL
write.csv(summ,
          file.path(results_dir, "uc1_analytical_summary.csv"),
          row.names = FALSE)
cat(sprintf("Wrote %s/uc1_analytical_*.csv\n", results_dir))
