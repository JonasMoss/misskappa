#!/usr/bin/env Rscript
#
# Experiment 14: MCAR verification sims for the closed-form quadratic
# kappa estimator. Continuous (multivariate normal) ratings, MCAR
# missingness, empirical and normal-Gamma Wald inference. Reports bias,
# SE accuracy, and Wald coverage for Fleiss and Conger on a small grid.

suppressPackageStartupMessages({
  library(misskappa)
  library(MASS)
})

usage <- function(status = 0L) {
  cat(paste0(
    "Usage: Rscript run_experiment.R [options]\n\n",
    "Options:\n",
    "  --help              Show this help and exit.\n",
    "  --smoke             Cheap mechanical check (small grid, few reps).\n",
    "  --n-grid CSV        Sample sizes. Default: 50,200,500.\n",
    "  --r-grid CSV        Rater counts. Default: 2,4.\n",
    "  --miss-grid CSV     Pairwise observation rates in (0,1]. Default: 0.6,1.0.\n",
    "  --reps N            Monte Carlo replicates per cell. Default: 1000.\n",
    "  --seed-base N       Base seed. Default: 141414.\n",
    "  --out-dir PATH      Output directory. Default: script-local results/.\n",
    "  --progress          One line per design cell.\n\n",
    "Examples:\n",
    "  Rscript run_experiment.R --smoke\n",
    "  Rscript run_experiment.R --reps 1000 --progress\n"
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
    n_grid = c(50L, 200L, 500L),
    r_grid = c(2L, 4L),
    miss_grid = c(0.6, 1.0),
    reps = 1000L,
    seed_base = 141414L,
    out_dir = NULL,
    progress = FALSE
  )
  i <- 1L
  while (i <= length(argv)) {
    arg <- argv[[i]]
    if (arg == "--help" || arg == "-h") usage(0L)
    if (arg == "--smoke") {
      opts$smoke <- TRUE
      i <- i + 1L; next
    }
    if (arg == "--progress") {
      opts$progress <- TRUE
      i <- i + 1L; next
    }
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
    opts$n_grid <- c(50L, 200L)
    opts$r_grid <- 2L
    opts$miss_grid <- c(0.6, 1.0)
    opts$reps <- 50L
  }
  if (is.na(opts$reps) || opts$reps < 1L) stop("--reps must be >= 1.", call. = FALSE)
  if (is.na(opts$seed_base)) stop("--seed-base must be an integer.", call. = FALSE)
  opts
}

# DGP definitions ------------------------------------------------------------
#
# DGP-A (exchangeable): all rater means 5, unit variance, common correlation
# rho = 0.7. Implies kappa_F = kappa_C = 0.7 for any R.
#
# DGP-B (heterogeneous means): rater means linearly spread by +/- 0.5 about 5,
# unit variance, common correlation 0.85. Yields kappa_F < kappa_C.

dgp_A <- function(R) {
  Sigma <- matrix(0.7, R, R)
  diag(Sigma) <- 1.0
  list(mu = rep(5.0, R), Sigma = Sigma)
}

dgp_B <- function(R) {
  Sigma <- matrix(0.85, R, R)
  diag(Sigma) <- 1.0
  list(mu = 5 + 0.5 * (seq_len(R) - (R + 1) / 2), Sigma = Sigma)
}

true_kappas <- function(mu, Sigma) {
  R <- length(mu)
  t1 <- sum(Sigma)
  t2 <- sum(diag(Sigma))
  t3 <- sum((mu - mean(mu))^2)
  c(
    Fleiss = (t1 - t2 - t3) / ((R - 1) * (t2 + t3)),
    Conger = (t1 - t2) / ((R - 1) * t2 + R * t3)
  )
}

# Per-replicate fit ----------------------------------------------------------

VALUES_PLACEHOLDER <- c(0, 1)  # unused by Fleiss / Conger; only c1 (BP) needs it

draw <- function(n, mu, Sigma, p_pair) {
  R <- length(mu)
  X <- MASS::mvrnorm(n, mu, Sigma)
  if (p_pair < 1.0) {
    p_cell <- sqrt(p_pair)
    M <- matrix(stats::rbinom(n * R, 1L, p_cell), n, R)
    X[M == 0L] <- NA_real_
    # Drop rows with < 2 observed entries (pairwise estimator drops them anyway,
    # but kappa_quadratic() does not currently handle all-NA rows gracefully).
    keep <- rowSums(!is.na(X)) >= 2L
    X <- X[keep, , drop = FALSE]
  }
  X
}

fit_one <- function(X, vcov_kind) {
  est <- tryCatch(
    misskappa::kappa_quadratic(X, values = VALUES_PLACEHOLDER, vcov = vcov_kind),
    error = function(e) NULL
  )
  if (is.null(est)) return(NULL)
  coefs <- est$estimates[c("Fleiss", "Conger")]
  ses <- sqrt(diag(est$vcov)[c("Fleiss", "Conger")])
  data.frame(
    coefficient = c("Fleiss", "Conger"),
    estimate = as.numeric(coefs),
    se = as.numeric(ses),
    stringsAsFactors = FALSE
  )
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
dgp_specs <- list(A_exch = dgp_A, B_het = dgp_B)
vcov_kinds <- c("empirical", "normal")

rows <- vector("list", 0L)
pos <- 1L

for (dgp_name in names(dgp_specs)) {
  dgp_fn <- dgp_specs[[dgp_name]]
  for (R in opts$r_grid) {
    spec <- dgp_fn(R)
    truth <- true_kappas(spec$mu, spec$Sigma)
    for (n in opts$n_grid) {
      for (p_pair in opts$miss_grid) {
        if (opts$progress) {
          message(sprintf(
            "dgp=%s R=%d n=%d p_pair=%.2f reps=%d",
            dgp_name, R, n, p_pair, opts$reps
          ))
        }
        for (rep in seq_len(opts$reps)) {
          seed <- opts$seed_base +
            1e7 * match(dgp_name, names(dgp_specs)) +
            1e5 * R +
            1e2 * round(100 * p_pair) +
            rep
          set.seed(seed)
          X <- draw(n, spec$mu, spec$Sigma, p_pair)
          if (nrow(X) < 5L) next  # degenerate; skip
          for (vk in vcov_kinds) {
            fit <- fit_one(X, vk)
            if (is.null(fit)) next
            for (i in seq_len(nrow(fit))) {
              co <- fit$coefficient[[i]]
              est <- fit$estimate[[i]]
              se <- fit$se[[i]]
              if (!is.finite(est) || !is.finite(se)) next
              rows[[pos]] <- data.frame(
                dgp = dgp_name,
                R = R,
                n = n,
                p_pair = p_pair,
                rep = rep,
                vcov_kind = vk,
                coefficient = co,
                truth = unname(truth[co]),
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
  }
}

estimates <- do.call(rbind, rows)
estimates$covered <- estimates$lower <= estimates$truth &
  estimates$truth <= estimates$upper

key <- with(estimates, interaction(dgp, R, n, p_pair, vcov_kind, coefficient,
                                   drop = TRUE, sep = "|"))
pieces <- split(estimates, key)
summary <- do.call(rbind, lapply(pieces, function(d) {
  data.frame(
    dgp = d$dgp[1L],
    R = d$R[1L],
    n = d$n[1L],
    p_pair = d$p_pair[1L],
    vcov_kind = d$vcov_kind[1L],
    coefficient = d$coefficient[1L],
    truth = d$truth[1L],
    reps = nrow(d),
    bias = mean(d$estimate) - d$truth[1L],
    mc_sd = stats::sd(d$estimate),
    mean_se = mean(d$se),
    se_over_sd = mean(d$se) / stats::sd(d$estimate),
    cov95 = mean(d$covered),
    stringsAsFactors = FALSE
  )
}))
summary <- summary[order(summary$dgp, summary$coefficient, summary$vcov_kind,
                         summary$R, summary$n, summary$p_pair), ]

dir.create(opts$out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(estimates, file.path(opts$out_dir, "estimates.csv"), row.names = FALSE)
write.csv(summary, file.path(opts$out_dir, "summary.csv"), row.names = FALSE)

metadata <- data.frame(
  key = c("experiment", "n_grid", "r_grid", "miss_grid", "reps",
          "seed_base", "dgps", "vcov_kinds", "values_placeholder",
          "R_version", "misskappa_version"),
  value = c(
    "14-quadratic-mcar-verification",
    paste(opts$n_grid, collapse = ","),
    paste(opts$r_grid, collapse = ","),
    paste(opts$miss_grid, collapse = ","),
    as.character(opts$reps),
    as.character(opts$seed_base),
    paste(names(dgp_specs), collapse = ","),
    paste(vcov_kinds, collapse = ","),
    paste(VALUES_PLACEHOLDER, collapse = ","),
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
