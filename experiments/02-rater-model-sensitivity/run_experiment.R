#!/usr/bin/env Rscript
#
# Experiment 02: rater model sensitivity.
#
# Compares the current latent-truth-plus-guess simulation model used by the
# paper with Dawid-Skene per-rater confusion-matrix variants of the same A/B/C
# design. Writes rectangular CSVs under results/.

usage <- function(status = 0L) {
  cat(paste0(
    "Usage: Rscript run_experiment.R [options]\n\n",
    "Options:\n",
    "  --help              Show this help and exit.\n",
    "  --smoke             Cheap mechanical check: small n, n_truth, reps.\n",
    "  --reps N            Monte Carlo replicates per cell.\n",
    "  --n N               Sample size per replicate.\n",
    "  --n-truth N         Complete-data sample size for truth estimates.\n",
    "  --seed-base N       Base seed for deterministic runs.\n",
    "  --out-dir PATH      Output directory. Default: script-local results/.\n\n",
    "Examples:\n",
    "  Rscript run_experiment.R --smoke\n",
    "  Rscript run_experiment.R --reps 200 --n 4000 --n-truth 200000\n"
  ))
  quit(save = "no", status = status)
}

parse_args <- function(argv) {
  opts <- list(
    smoke = FALSE,
    reps = 50L,
    n = 1000L,
    n_truth = 50000L,
    seed_base = 200200L,
    out_dir = NULL
  )
  i <- 1L
  while (i <= length(argv)) {
    arg <- argv[[i]]
    if (arg == "--help" || arg == "-h") usage(0L)
    if (arg == "--smoke") {
      opts$smoke <- TRUE
      i <- i + 1L
      next
    }
    needs_value <- c("--reps", "--n", "--n-truth", "--seed-base", "--out-dir")
    if (arg %in% needs_value) {
      if (i == length(argv)) stop(arg, " needs a value.", call. = FALSE)
      val <- argv[[i + 1L]]
      if (arg == "--reps") opts$reps <- as.integer(val)
      if (arg == "--n") opts$n <- as.integer(val)
      if (arg == "--n-truth") opts$n_truth <- as.integer(val)
      if (arg == "--seed-base") opts$seed_base <- as.integer(val)
      if (arg == "--out-dir") opts$out_dir <- val
      i <- i + 2L
      next
    }
    stop("Unknown argument: ", arg, call. = FALSE)
  }
  if (opts$smoke) {
    opts$reps <- min(opts$reps, 2L)
    opts$n <- min(opts$n, 120L)
    opts$n_truth <- min(opts$n_truth, 1200L)
  }
  if (is.na(opts$reps) || opts$reps < 1L) stop("--reps must be >= 1.", call. = FALSE)
  if (is.na(opts$n) || opts$n < 2L) stop("--n must be >= 2.", call. = FALSE)
  if (is.na(opts$n_truth) || opts$n_truth < 2L) stop("--n-truth must be >= 2.", call. = FALSE)
  if (is.na(opts$seed_base)) stop("--seed-base must be an integer.", call. = FALSE)
  opts
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE))
} else {
  getwd()
}
if (is.null(opts$out_dir)) opts$out_dir <- file.path(script_dir, "results")

suppressPackageStartupMessages({
  library(misskappa)
})

C <- 5L
R <- 6L
methods <- c("available", "ipw", "fiml", "gwet")
coefficients <- c("Conger", "Fleiss", "Brennan-Prediger")
weight_main <- "identity"

p <- colSums(as.matrix(dat.fleiss1971)) / sum(as.matrix(dat.fleiss1971))

clamp <- function(x, lo = 0, hi = 1) pmin(pmax(x, lo), hi)

paper_guess_mat <- matrix(rep(p, each = R), nrow = R, byrow = TRUE)

make_confusion <- function(diagonal, base = p, direction = 0, temperature = 1.1) {
  if (length(diagonal) == 1L) diagonal <- rep(diagonal, C)
  out <- matrix(0, nrow = C, ncol = C)
  for (truth in seq_len(C)) {
    cats <- seq_len(C)
    off <- base * exp(-abs(cats - truth) / temperature + direction * (cats - truth))
    off[truth] <- 0
    if (sum(off) <= 0) off[-truth] <- 1
    off <- off / sum(off)
    out[truth, ] <- (1 - diagonal[truth]) * off
    out[truth, truth] <- diagonal[truth]
  }
  out
}

shared_confusion <- make_confusion(0.92, direction = 0.0, temperature = 0.95)
difficult_confusion <- make_confusion(
  c(0.96, 0.93, 0.88, 0.78, 0.64),
  direction = 0.0,
  temperature = 0.95
)

ds_confusions_b <- lapply(seq_len(R), function(j) {
  make_confusion(
    diagonal = c(0.98, 0.96, 0.92, 0.86, 0.74, 0.62)[j],
    direction = c(-0.45, -0.25, -0.10, 0.10, 0.25, 0.45)[j],
    temperature = c(0.85, 0.90, 0.95, 1.05, 1.15, 1.25)[j]
  )
})

dgps <- list(
  list(
    label = "A",
    name = "Exchangeable + MCAR",
    guess = list(
      rho_base = rep(0.92, R),
      rho_truth_mult = rep(1, C),
      guess = paper_guess_mat,
      missing = "mcar",
      pi_rater = c(0.95, 0.90, 0.80, 0.65, 0.50, 0.35)
    ),
    dawid_skene = list(
      confusion = rep(list(shared_confusion), R),
      missing = "mcar",
      pi_rater = c(0.95, 0.90, 0.80, 0.65, 0.50, 0.35)
    )
  ),
  list(
    label = "B",
    name = "Non-exchangeable + MCAR",
    guess = list(
      rho_base = c(0.98, 0.96, 0.92, 0.85, 0.72, 0.60),
      rho_truth_mult = rep(1, C),
      guess = paper_guess_mat,
      missing = "mcar",
      pi_rater = c(0.95, 0.85, 0.70, 0.45, 0.25, 0.15)
    ),
    dawid_skene = list(
      confusion = ds_confusions_b,
      missing = "mcar",
      pi_rater = c(0.95, 0.85, 0.70, 0.45, 0.25, 0.15)
    )
  ),
  list(
    label = "C",
    name = "Difficulty + MAR",
    guess = list(
      rho_base = rep(0.92, R),
      rho_truth_mult = c(1.05, 1.00, 0.95, 0.85, 0.70),
      guess = paper_guess_mat,
      missing = "mar_truth",
      pi_truth = c(0.95, 0.90, 0.80, 0.55, 0.25)
    ),
    dawid_skene = list(
      confusion = rep(list(difficult_confusion), R),
      missing = "mar_previous_observed",
      pi_prev_rating = c(0.95, 0.90, 0.80, 0.55, 0.30),
      pi_first = 0.95,
      pi_fallback = 0.65
    )
  )
)

simulate_guess_complete <- function(n, spec) {
  truth <- sample.int(C, n, replace = TRUE, prob = p)
  x_star <- matrix(NA_integer_, nrow = n, ncol = R)
  for (j in seq_len(R)) {
    rho_i <- clamp(spec$rho_base[j] * spec$rho_truth_mult[truth])
    correct <- stats::rbinom(n, 1, prob = rho_i) == 1L
    guessed <- sample.int(C, n, replace = TRUE, prob = spec$guess[j, ])
    x_star[, j] <- ifelse(correct, truth, guessed)
  }
  list(x_star = x_star, truth = truth)
}

sample_from_confusion <- function(truth, confusion) {
  x <- integer(length(truth))
  for (cat in seq_len(C)) {
    idx <- which(truth == cat)
    if (length(idx) > 0L) {
      x[idx] <- sample.int(C, length(idx), replace = TRUE, prob = confusion[cat, ])
    }
  }
  x
}

simulate_ds_complete <- function(n, spec) {
  truth <- sample.int(C, n, replace = TRUE, prob = p)
  x_star <- matrix(NA_integer_, nrow = n, ncol = R)
  for (j in seq_len(R)) {
    x_star[, j] <- sample_from_confusion(truth, spec$confusion[[j]])
  }
  list(x_star = x_star, truth = truth)
}

apply_missing_mcar <- function(x_star, pi_rater) {
  x <- x_star
  for (j in seq_len(R)) {
    observed <- stats::runif(nrow(x_star)) < pi_rater[j]
    x[!observed, j] <- NA_integer_
  }
  x
}

apply_missing_mar_truth <- function(x_star, truth, pi_truth) {
  x <- x_star
  p_i <- pi_truth[truth]
  observed <- matrix(stats::runif(nrow(x_star) * R) < rep(p_i, times = R),
                     nrow = nrow(x_star), ncol = R)
  x[!observed] <- NA_integer_
  x
}

apply_missing_mar_previous_observed <- function(x_star, pi_prev_rating,
                                                pi_first, pi_fallback) {
  n <- nrow(x_star)
  x <- x_star
  observed_first <- stats::runif(n) < pi_first
  x[!observed_first, 1L] <- NA_integer_
  if (R == 1L) return(x)

  for (j in 2:R) {
    prev <- x[, j - 1L]
    p_obs <- ifelse(is.na(prev), pi_fallback, pi_prev_rating[prev])
    observed <- stats::runif(n) < p_obs
    x[!observed, j] <- NA_integer_
  }
  x
}

simulate_complete <- function(n, family, spec) {
  if (family == "guess") return(simulate_guess_complete(n, spec))
  if (family == "dawid_skene") return(simulate_ds_complete(n, spec))
  stop("Unknown model family: ", family, call. = FALSE)
}

apply_missing <- function(dat, spec) {
  switch(
    spec$missing,
    mcar = apply_missing_mcar(dat$x_star, spec$pi_rater),
    mar_truth = apply_missing_mar_truth(dat$x_star, dat$truth, spec$pi_truth),
    mar_previous_observed = apply_missing_mar_previous_observed(
      dat$x_star, spec$pi_prev_rating, spec$pi_first, spec$pi_fallback
    ),
    stop("Unknown missingness mechanism: ", spec$missing, call. = FALSE)
  )
}

kappa_estimates <- function(x, method) {
  args <- list(x = x, method = method, weight = weight_main)
  if (method == "fiml") {
    args$em_options <- list(max_iter = 50000L, tol = 1e-7)
  }
  coef(do.call(misskappa::kappa, args))
}

truth_estimate <- function(n_truth, family, spec, seed) {
  set.seed(seed)
  dat <- simulate_complete(n_truth, family, spec)
  kappa_estimates(dat$x_star, method = "available")
}

fit_one <- function(x, method) {
  tryCatch(
    list(estimates = kappa_estimates(x, method), error = NA_character_),
    error = function(e) {
      est <- setNames(rep(NA_real_, length(coefficients)), coefficients)
      list(estimates = est, error = conditionMessage(e))
    }
  )
}

replicate_cell <- function(n, family, dgp, spec, truth, rep_id, seed) {
  set.seed(seed)
  dat <- simulate_complete(n, family, spec)
  x <- apply_missing(dat, spec)
  rows <- vector("list", length(methods) * length(coefficients))
  pos <- 1L
  for (method in methods) {
    fit <- fit_one(x, method)
    for (coefficient in coefficients) {
      rows[[pos]] <- data.frame(
        model_family = family,
        dgp = dgp$label,
        dgp_name = dgp$name,
        n = n,
        rep = rep_id,
        seed = seed,
        method = method,
        coefficient = coefficient,
        estimate = unname(fit$estimates[[coefficient]]),
        truth = unname(truth[[coefficient]]),
        error = ifelse(is.na(fit$error), "", fit$error),
        stringsAsFactors = FALSE
      )
      pos <- pos + 1L
    }
  }
  do.call(rbind, rows)
}

split_keys <- function(data, keys) {
  interaction(data[, keys], drop = TRUE, lex.order = TRUE)
}

summarize_results <- function(per_rep) {
  keys <- c("model_family", "dgp", "dgp_name", "method", "coefficient")
  groups <- split(per_rep, split_keys(per_rep, keys))
  out <- lapply(groups, function(g) {
    valid <- is.finite(g$estimate)
    truth <- unique(g$truth)
    if (length(truth) != 1L) stop("Non-unique truth in summary group.", call. = FALSE)
    est <- g$estimate[valid]
    bias <- if (length(est) > 0L) mean(est) - truth else NA_real_
    data.frame(
      model_family = g$model_family[[1L]],
      dgp = g$dgp[[1L]],
      dgp_name = g$dgp_name[[1L]],
      method = g$method[[1L]],
      coefficient = g$coefficient[[1L]],
      n = unique(g$n)[[1L]],
      reps = length(unique(g$rep)),
      n_valid = sum(valid),
      n_error = sum(!valid),
      truth = truth,
      mean = if (length(est) > 0L) mean(est) else NA_real_,
      sd = if (length(est) > 1L) stats::sd(est) else NA_real_,
      bias = bias,
      abs_bias = abs(bias),
      mse = if (length(est) > 0L) mean((est - truth)^2) else NA_real_,
      error_rate = mean(!valid),
      stringsAsFactors = FALSE
    )
  })
  summary <- do.call(rbind, out)
  summary <- summary[order(summary$model_family, summary$dgp,
                           summary$coefficient, summary$method), ]
  rownames(summary) <- NULL

  summary$rank_abs_bias <- NA_integer_
  summary$rank_sd <- NA_integer_
  rank_keys <- c("model_family", "dgp", "coefficient")
  for (key in unique(split_keys(summary, rank_keys))) {
    idx <- which(split_keys(summary, rank_keys) == key)
    summary$rank_abs_bias[idx] <- rank(summary$abs_bias[idx], ties.method = "min", na.last = "keep")
    summary$rank_sd[idx] <- rank(summary$sd[idx], ties.method = "min", na.last = "keep")
  }
  summary
}

metadata <- data.frame(
  key = c(
    "timestamp", "smoke", "reps", "n", "n_truth", "seed_base",
    "C", "R", "methods", "weight", "misskappa_version",
    "r_version", "model_families", "dgp_labels"
  ),
  value = c(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
    as.character(opts$smoke),
    as.character(opts$reps),
    as.character(opts$n),
    as.character(opts$n_truth),
    as.character(opts$seed_base),
    as.character(C),
    as.character(R),
    paste(methods, collapse = ";"),
    weight_main,
    as.character(utils::packageVersion("misskappa")),
    paste(R.version$major, R.version$minor, sep = "."),
    "guess;dawid_skene",
    paste(vapply(dgps, `[[`, character(1), "label"), collapse = ";")
  ),
  stringsAsFactors = FALSE
)

families <- c("guess", "dawid_skene")
truth_lookup <- list()
per_rep_rows <- list()
cell_id <- 0L

for (family in families) {
  truth_lookup[[family]] <- list()
  for (dgp in dgps) {
    spec <- dgp[[family]]
    cell_id <- cell_id + 1L
    truth_seed <- opts$seed_base + 1000000L + cell_id
    truth <- truth_estimate(opts$n_truth, family, spec, truth_seed)
    truth_lookup[[family]][[dgp$label]] <- truth
    for (rep_id in seq_len(opts$reps)) {
      seed <- opts$seed_base + cell_id * 10000L + rep_id
      per_rep_rows[[length(per_rep_rows) + 1L]] <- replicate_cell(
        opts$n, family, dgp, spec, truth, rep_id, seed
      )
    }
  }
}

per_rep <- do.call(rbind, per_rep_rows)
summary <- summarize_results(per_rep)

dir.create(opts$out_dir, recursive = TRUE, showWarnings = FALSE)
metadata_path <- file.path(opts$out_dir, "metadata.csv")
per_rep_path <- file.path(opts$out_dir, "per_rep.csv")
summary_path <- file.path(opts$out_dir, "summary.csv")

write.csv(metadata, metadata_path, row.names = FALSE)
write.csv(per_rep, per_rep_path, row.names = FALSE)
write.csv(summary, summary_path, row.names = FALSE)

cat("Wrote ", metadata_path, "\n", sep = "")
cat("Wrote ", per_rep_path, "\n", sep = "")
cat("Wrote ", summary_path, "\n", sep = "")
