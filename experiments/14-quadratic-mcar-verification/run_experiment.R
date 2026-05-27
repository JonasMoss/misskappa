#!/usr/bin/env Rscript
#
# Experiment 14: MCAR verification sims for the closed-form quadratic
# kappa estimator. Perreault-Leigh judge-skill DGP, MCAR missingness,
# empirical and normal Gamma. Reports bias, SE accuracy, Wald coverage,
# and RMSE for Fleiss and Conger. Includes a paired comparison between
# the pairwise-available estimator and listwise deletion.

suppressPackageStartupMessages({
  library(misskappa)
})

usage <- function(status = 0L) {
  cat(paste0(
    "Usage: Rscript run_experiment.R [options]\n\n",
    "Options:\n",
    "  --help              Show this help and exit.\n",
    "  --smoke             Cheap mechanical check (small grid, few reps).\n",
    "  --n-grid CSV        Sample sizes. Default: 10,50,200,500.\n",
    "  --r-grid CSV        Rater counts. Default: 2,4.\n",
    "  --miss-grid CSV     Pairwise observation rates in (0,1]. Default: 0.6.\n",
    "  --reps N            Monte Carlo replicates per cell. Default: 1000.\n",
    "  --seed-base N       Base seed. Default: 141414.\n",
    "  --out-dir PATH      Output directory. Default: script-local results/.\n",
    "  --progress          One line per design cell.\n\n",
    "DGP: Perreault-Leigh judge-skill model with knowledge probability\n",
    "p_k = 0.8 and integer category scores {-2,-1,0,1,2}. Implies\n",
    "kappa_F = kappa_C = 0.8 for any R.\n"
  ))
  quit(save = "no", status = status)
}

parse_int_csv <- function(x, min_value, arg_name) {
  out <- as.integer(strsplit(x, ",", fixed = TRUE)[[1L]])
  if (any(is.na(out)) || any(out < min_value)) {
    stop(arg_name, " must contain integers >= ", min_value, ".", call. = FALSE)
  }
  unique(out)
}

parse_num_csv <- function(x, arg_name) {
  out <- as.numeric(strsplit(x, ",", fixed = TRUE)[[1L]])
  if (any(is.na(out)) || any(out <= 0) || any(out > 1)) {
    stop(arg_name, " must contain numbers in (0, 1].", call. = FALSE)
  }
  unique(out)
}

parse_args <- function(argv) {
  opts <- list(
    smoke = FALSE,
    n_grid = c(10L, 50L, 200L, 500L),
    r_grid = c(2L, 4L),
    miss_grid = 0.6,
    reps = 1000L,
    seed_base = 141414L,
    out_dir = NULL,
    progress = FALSE
  )
  i <- 1L
  while (i <= length(argv)) {
    arg <- argv[[i]]
    if (arg == "--help" || arg == "-h") usage(0L)
    if (arg == "--smoke") { opts$smoke <- TRUE; i <- i + 1L; next }
    if (arg == "--progress") { opts$progress <- TRUE; i <- i + 1L; next }
    needs_value <- c("--n-grid", "--r-grid", "--miss-grid", "--reps",
                     "--seed-base", "--out-dir")
    if (arg %in% needs_value) {
      if (i == length(argv)) stop(arg, " needs a value.", call. = FALSE)
      val <- argv[[i + 1L]]
      if (arg == "--n-grid") opts$n_grid <- parse_int_csv(val, 5L, "--n-grid")
      if (arg == "--r-grid") opts$r_grid <- parse_int_csv(val, 2L, "--r-grid")
      if (arg == "--miss-grid") opts$miss_grid <- parse_num_csv(val, "--miss-grid")
      if (arg == "--reps") opts$reps <- as.integer(val)
      if (arg == "--seed-base") opts$seed_base <- as.integer(val)
      if (arg == "--out-dir") opts$out_dir <- val
      i <- i + 2L; next
    }
    stop("Unknown argument: ", arg, call. = FALSE)
  }
  if (opts$smoke) {
    opts$n_grid <- c(10L, 50L)
    opts$r_grid <- 2L
    opts$miss_grid <- 0.6
    opts$reps <- 50L
  }
  if (is.na(opts$reps) || opts$reps < 1L) stop("--reps must be >= 1.", call. = FALSE)
  if (is.na(opts$seed_base)) stop("--seed-base must be an integer.", call. = FALSE)
  opts
}

# Perreault-Leigh DGP --------------------------------------------------------
#
# Categories C = {-2, -1, 0, 1, 2}, knowledge probability p_k = 0.9 (matches
# Moss 2024 Example 3). For each subject draw a truth uniformly from C; each
# rater independently reports the truth with probability p_k or a uniform
# guess otherwise. The marginal is uniform on C, so the (mu, Sigma)-based
# Conger and Fleiss kappa both equal p_k^2 under exchangeable raters:
# Cov(X_j, X_k) = p_k^2 Var(T) and Var(X_j) = Var(T). At p_k = 0.9 the truth
# is therefore 0.81 (Moss 2024 Example 3 reports 0.816 following Perreault &
# Leigh 1989; the small gap is a normalization detail of their derivation).

VALUES <- c(-2, -1, 0, 1, 2)
KNOWLEDGE_PROB <- 0.9
TRUTH <- KNOWLEDGE_PROB^2  # exchangeable PL: kappa_F = kappa_C = p_k^2

simulate_perreault_leigh <- function(n, R) {
  C <- length(VALUES)
  truth_idx <- sample.int(C, n, replace = TRUE)
  x <- matrix(NA_integer_, nrow = n, ncol = R)
  for (j in seq_len(R)) {
    knows <- stats::runif(n) < KNOWLEDGE_PROB
    guesses <- sample.int(C, n, replace = TRUE)
    x[, j] <- ifelse(knows, truth_idx, guesses)
  }
  matrix(VALUES[x], nrow = n, ncol = R)
}

apply_mcar <- function(X, p_pair) {
  if (p_pair >= 1.0) return(X)
  n <- nrow(X); R <- ncol(X)
  p_cell <- sqrt(p_pair)
  M <- matrix(stats::rbinom(n * R, 1L, p_cell), n, R)
  X[M == 0L] <- NA_real_
  X
}

# Fit pair (raw matrix) and listwise (complete.cases subset) under both Gamma
# variants. Returns a data.frame with one row per (estimator, vcov_kind,
# coefficient) successful fit, NULL if all fits fail.

fit_paired <- function(X, vcov_kinds = c("empirical", "normal")) {
  rows <- list(); pos <- 1L
  # Pair: keep every row with at least one observed entry. Singleton-observed
  # rows still contribute to rater means under the pairwise-available
  # estimator. Listwise: drop any row with at least one NA.
  X_pair <- X[rowSums(!is.na(X)) >= 1L, , drop = FALSE]
  X_lw <- X[stats::complete.cases(X), , drop = FALSE]
  fit_at <- function(Xfit, estimator) {
    if (nrow(Xfit) < 3L) return(NULL)
    out <- list()
    for (vk in vcov_kinds) {
      fit <- tryCatch(
        misskappa::kappa_quadratic(Xfit, values = VALUES, vcov = vk),
        error = function(e) NULL
      )
      if (is.null(fit)) next
      coefs <- fit$estimates[c("Fleiss", "Conger")]
      ses <- sqrt(diag(fit$vcov)[c("Fleiss", "Conger")])
      for (k in c("Fleiss", "Conger")) {
        if (!is.finite(coefs[k]) || !is.finite(ses[k])) next
        out[[length(out) + 1L]] <- data.frame(
          estimator = estimator, vcov_kind = vk, coefficient = k,
          estimate = unname(coefs[k]), se = unname(ses[k]),
          n_used = nrow(Xfit),
          stringsAsFactors = FALSE
        )
      }
    }
    if (length(out) == 0L) NULL else do.call(rbind, out)
  }
  pair_rows <- fit_at(X_pair, "pair")
  lw_rows <- fit_at(X_lw, "listwise")
  if (is.null(pair_rows) && is.null(lw_rows)) return(NULL)
  do.call(rbind, list(pair_rows, lw_rows))
}

# Driver ---------------------------------------------------------------------

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE))
} else {
  getwd()
}
opts <- parse_args(commandArgs(trailingOnly = TRUE))
if (is.null(opts$out_dir)) opts$out_dir <- file.path(script_dir, "results")

z975 <- stats::qnorm(0.975)
rows <- vector("list", 0L); pos <- 1L

for (R in opts$r_grid) {
  for (n in opts$n_grid) {
    for (p_pair in opts$miss_grid) {
      if (opts$progress) {
        message(sprintf("R=%d n=%d p_pair=%.2f reps=%d", R, n, p_pair, opts$reps))
      }
      for (rep in seq_len(opts$reps)) {
        seed <- opts$seed_base + 1e7 * R + 1e4 * n + 1e2 * round(100 * p_pair) + rep
        set.seed(seed)
        X <- simulate_perreault_leigh(n, R)
        X <- apply_mcar(X, p_pair)
        fits <- fit_paired(X)
        if (is.null(fits)) next
        for (i in seq_len(nrow(fits))) {
          est <- fits$estimate[[i]]
          se <- fits$se[[i]]
          rows[[pos]] <- data.frame(
            R = R, n = n, p_pair = p_pair, rep = rep,
            estimator = fits$estimator[[i]],
            vcov_kind = fits$vcov_kind[[i]],
            coefficient = fits$coefficient[[i]],
            n_used = fits$n_used[[i]],
            truth = TRUTH,
            estimate = est,
            se = se,
            lower = est - z975 * se,
            upper = est + z975 * se,
            stringsAsFactors = FALSE
          )
          pos <- pos + 1L
        }
      }
    }
  }
}

estimates <- do.call(rbind, rows)
estimates$covered <- estimates$lower <= estimates$truth &
  estimates$truth <= estimates$upper

key <- with(estimates, interaction(R, n, p_pair, estimator, vcov_kind,
                                   coefficient, drop = TRUE, sep = "|"))
pieces <- split(estimates, key)
summary <- do.call(rbind, lapply(pieces, function(d) {
  data.frame(
    R = d$R[1L], n = d$n[1L], p_pair = d$p_pair[1L],
    estimator = d$estimator[1L], vcov_kind = d$vcov_kind[1L],
    coefficient = d$coefficient[1L],
    truth = d$truth[1L],
    reps = nrow(d),
    mean_n_used = mean(d$n_used),
    bias = mean(d$estimate) - d$truth[1L],
    mc_sd = stats::sd(d$estimate),
    mean_se = mean(d$se),
    se_over_sd = mean(d$se) / stats::sd(d$estimate),
    rmse = sqrt(mean((d$estimate - d$truth[1L])^2)),
    cov95 = mean(d$covered),
    stringsAsFactors = FALSE
  )
}))
summary <- summary[order(summary$coefficient, summary$estimator,
                         summary$vcov_kind, summary$R, summary$n,
                         summary$p_pair), ]

dir.create(opts$out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(estimates, file.path(opts$out_dir, "estimates.csv"), row.names = FALSE)
write.csv(summary, file.path(opts$out_dir, "summary.csv"), row.names = FALSE)

metadata <- data.frame(
  key = c("experiment", "dgp", "knowledge_prob", "values",
          "n_grid", "r_grid", "miss_grid", "reps", "seed_base",
          "estimators", "vcov_kinds",
          "R_version", "misskappa_version"),
  value = c(
    "14-quadratic-mcar-verification",
    "perreault-leigh (Moss 2024 Sec 5.1.1)",
    as.character(KNOWLEDGE_PROB),
    paste(VALUES, collapse = ","),
    paste(opts$n_grid, collapse = ","),
    paste(opts$r_grid, collapse = ","),
    paste(opts$miss_grid, collapse = ","),
    as.character(opts$reps),
    as.character(opts$seed_base),
    "pair,listwise",
    "empirical,normal",
    R.version.string,
    as.character(utils::packageVersion("misskappa"))
  ),
  stringsAsFactors = FALSE
)
write.csv(metadata, file.path(opts$out_dir, "metadata.csv"), row.names = FALSE)

cat("Wrote:\n")
cat(" ", file.path(opts$out_dir, "estimates.csv"), "\n")
cat(" ", file.path(opts$out_dir, "summary.csv"), "\n")
cat(" ", file.path(opts$out_dir, "metadata.csv"), "\n")
