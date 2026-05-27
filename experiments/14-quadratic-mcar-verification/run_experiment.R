#!/usr/bin/env Rscript
#
# Experiment 14: MCAR verification sims for the closed-form quadratic
# kappa estimator. Heterogeneous judge-skill DGP with per-rater
# knowledge probabilities and per-rater guessing distributions, MCAR
# missingness, empirical and normal Gamma. Reports bias, SE accuracy,
# Wald coverage, and RMSE for Fleiss and Conger. Includes a paired
# comparison between the pairwise-available estimator and listwise
# deletion. The per-rater guessing distributions are intentionally
# asymmetric across raters so that the rater means differ; this
# separates kappa_F from kappa_C, which would otherwise collapse under
# the exchangeable Perreault-Leigh special case.

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
    "DGP: heterogeneous judge-skill model on category scores {-2,-1,0,1,2}.\n",
    "Per-rater knowledge probabilities and per-rater guessing distributions\n",
    "produce different rater means, separating kappa_F from kappa_C.\n"
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

# Heterogeneous judge-skill DGP ----------------------------------------------
#
# Categories C = {-2, -1, 0, 1, 2}, uniform truth distribution. Each rater
# j has its own knowledge probability s_j and guessing distribution G[j, ].
# With prob s_j rater j reports the truth; otherwise rater j draws from
# G[j, ]. The bias rows below are intentionally asymmetric across raters
# so that rater means differ -- this gives kappa_F < kappa_C.
#
# For R = 4 we use all four bias rows (raters span low-bias to high-bias).
# For R = 2 we use rows 2 and 3 -- the two mildly biased raters -- so that
# the kappa values stay in a comparable range to the R = 4 cells. This
# yields the population truths kappa_F = 0.503, kappa_C = 0.525 at R = 2
# and kappa_F = 0.428, kappa_C = 0.445 at R = 4 (computed in closed form
# from the population (mu, Sigma) below).

VALUES <- c(-2, -1, 0, 1, 2)
TRUE_DIST <- rep(1 / length(VALUES), length(VALUES))

GUESS_FULL <- rbind(
  c(0.70, 0.30, 0.00, 0.00, 0.00),
  c(0.30, 0.50, 0.20, 0.00, 0.00),
  c(0.00, 0.00, 0.20, 0.50, 0.30),
  c(0.00, 0.00, 0.00, 0.30, 0.70)
)
SKILL_FULL <- c(0.65, 0.70, 0.75, 0.70)
SELECT_FOR_R <- list(`2` = c(2L, 3L), `4` = 1:4)

dgp_params <- function(R) {
  key <- as.character(R)
  if (is.null(SELECT_FOR_R[[key]])) {
    stop("Unsupported R: ", R, ". Add an entry to SELECT_FOR_R.", call. = FALSE)
  }
  idx <- SELECT_FOR_R[[key]]
  list(s = SKILL_FULL[idx], G = GUESS_FULL[idx, , drop = FALSE])
}

simulate_jsm <- function(n, R) {
  pars <- dgp_params(R)
  x_idx <- misskappa::sim$jsm(n = n, s = pars$s, model = "general",
                              true_dist = TRUE_DIST, guessing_dist = pars$G)
  matrix(VALUES[x_idx], nrow = n, ncol = R)
}

true_kappas <- function(R) {
  pars <- dgp_params(R); s <- pars$s; G <- pars$G
  ET <- sum(VALUES * TRUE_DIST)
  ET2 <- sum(VALUES^2 * TRUE_DIST)
  EG <- as.numeric(G %*% VALUES)
  EG2 <- as.numeric(G %*% (VALUES^2))
  mu <- s * ET + (1 - s) * EG
  EV <- s * (ET2 - ET^2) + (1 - s) * (EG2 - EG^2)
  VE <- s * (ET - mu)^2 + (1 - s) * (EG - mu)^2
  sigma_jj <- EV + VE
  Sigma <- matrix(0, R, R); diag(Sigma) <- sigma_jj
  for (j in 1:R) for (k in 1:R) if (j != k) {
    EXX <- s[j] * s[k] * ET2 +
           s[j] * (1 - s[k]) * ET * EG[k] +
           (1 - s[j]) * s[k] * EG[j] * ET +
           (1 - s[j]) * (1 - s[k]) * EG[j] * EG[k]
    Sigma[j, k] <- EXX - mu[j] * mu[k]
  }
  t1 <- sum(Sigma); t2 <- sum(diag(Sigma)); t3 <- sum((mu - mean(mu))^2)
  c(Fleiss = (t1 - t2 - t3) / ((R - 1) * (t2 + t3)),
    Conger = (t1 - t2) / ((R - 1) * t2 + R * t3))
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
truth_by_R <- lapply(opts$r_grid, function(R) true_kappas(R))
names(truth_by_R) <- as.character(opts$r_grid)

for (R in opts$r_grid) {
  truth <- truth_by_R[[as.character(R)]]
  for (n in opts$n_grid) {
    for (p_pair in opts$miss_grid) {
      if (opts$progress) {
        message(sprintf("R=%d n=%d p_pair=%.2f reps=%d", R, n, p_pair, opts$reps))
      }
      for (rep in seq_len(opts$reps)) {
        seed <- opts$seed_base + 1e7 * R + 1e4 * n + 1e2 * round(100 * p_pair) + rep
        set.seed(seed)
        X <- simulate_jsm(n, R)
        X <- apply_mcar(X, p_pair)
        fits <- fit_paired(X)
        if (is.null(fits)) next
        for (i in seq_len(nrow(fits))) {
          coef <- fits$coefficient[[i]]
          est <- fits$estimate[[i]]
          se <- fits$se[[i]]
          rows[[pos]] <- data.frame(
            R = R, n = n, p_pair = p_pair, rep = rep,
            estimator = fits$estimator[[i]],
            vcov_kind = fits$vcov_kind[[i]],
            coefficient = coef,
            n_used = fits$n_used[[i]],
            truth = unname(truth[coef]),
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
  key = c("experiment", "dgp", "skills_full", "values",
          "n_grid", "r_grid", "miss_grid", "reps", "seed_base",
          "estimators", "vcov_kinds",
          "R_version", "misskappa_version"),
  value = c(
    "14-quadratic-mcar-verification",
    "heterogeneous judge-skill (sim$jsm 'general'); per-rater bias",
    paste(SKILL_FULL, collapse = ","),
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
