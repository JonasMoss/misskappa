#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(misskappa)
})

usage <- function(status = 0L) {
  cat(paste0(
    "Usage: Rscript run_experiment.R [options]\n\n",
    "Options:\n",
    "  --help, -h             Show this help and exit.\n",
    "  --smoke                Cheap check: P2 only, n=120,400, reps=2.\n",
    "  --reps N               Replicates per scaling cell. Default: 8.\n",
    "  --n-grid LIST          Comma-separated n grid. Default: 500,2000,8000,32000.\n",
    "  --profiles LIST        Comma-separated profiles: p1,p2. Default: p1,p2.\n",
    "  --prune-reps N         Replicates per prune cell. Default: same as --reps.\n",
    "  --prune-n N            Sample size for prune sweep. Default: 2000.\n",
    "  --prune-tols LIST      Comma-separated prune_tols. Default: 1e-12,1e-9,1e-6.\n",
    "  --seed-base N          Base seed. Default: 505000.\n",
    "  --out-dir PATH         Output directory. Default: script-local results/.\n",
    "  --progress-every N     Print progress every N replicates. Default: 4.\n\n",
    "Examples:\n",
    "  Rscript run_experiment.R --smoke\n",
    "  Rscript run_experiment.R --reps 100 --prune-reps 40\n"
  ))
  quit(save = "no", status = status)
}

script_dir <- local({
  cmd <- commandArgs(FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) == 0L) getwd()
  else dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE))
})

parse_csv_int <- function(x, arg) {
  out <- as.integer(strsplit(x, ",", fixed = TRUE)[[1L]])
  if (any(is.na(out))) stop(arg, " must be a comma-separated integer list.", call. = FALSE)
  out
}

parse_csv_num <- function(x, arg) {
  out <- as.numeric(strsplit(x, ",", fixed = TRUE)[[1L]])
  if (any(is.na(out))) stop(arg, " must be a comma-separated numeric list.", call. = FALSE)
  out
}

parse_csv_chr <- function(x) {
  trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
}

parse_args <- function(argv) {
  opts <- list(
    smoke = FALSE,
    reps = 8L,
    n_grid = c(500L, 2000L, 8000L, 32000L),
    profiles = c("p1", "p2"),
    prune_reps = NA_integer_,
    prune_n = 2000L,
    prune_tols = c(1e-12, 1e-9, 1e-6),
    seed_base = 505000L,
    out_dir = file.path(script_dir, "results"),
    progress_every = 4L
  )
  explicit <- list(reps = FALSE, n_grid = FALSE, profiles = FALSE,
                   prune_reps = FALSE, prune_n = FALSE, prune_tols = FALSE)

  i <- 1L
  while (i <= length(argv)) {
    arg <- argv[[i]]
    if (arg %in% c("--help", "-h")) usage(0L)
    if (arg == "--smoke") {
      opts$smoke <- TRUE
      i <- i + 1L
      next
    }
    needs_value <- c(
      "--reps", "--n-grid", "--profiles", "--prune-reps", "--prune-n",
      "--prune-tols", "--seed-base", "--out-dir", "--progress-every"
    )
    if (arg %in% needs_value) {
      if (i == length(argv)) stop(arg, " needs a value.", call. = FALSE)
      val <- argv[[i + 1L]]
      if (arg == "--reps") {
        opts$reps <- as.integer(val)
        explicit$reps <- TRUE
      }
      if (arg == "--n-grid") {
        opts$n_grid <- parse_csv_int(val, arg)
        explicit$n_grid <- TRUE
      }
      if (arg == "--profiles") {
        opts$profiles <- parse_csv_chr(val)
        explicit$profiles <- TRUE
      }
      if (arg == "--prune-reps") {
        opts$prune_reps <- as.integer(val)
        explicit$prune_reps <- TRUE
      }
      if (arg == "--prune-n") {
        opts$prune_n <- as.integer(val)
        explicit$prune_n <- TRUE
      }
      if (arg == "--prune-tols") {
        opts$prune_tols <- parse_csv_num(val, arg)
        explicit$prune_tols <- TRUE
      }
      if (arg == "--seed-base") opts$seed_base <- as.integer(val)
      if (arg == "--out-dir") opts$out_dir <- val
      if (arg == "--progress-every") opts$progress_every <- as.integer(val)
      i <- i + 2L
      next
    }
    stop("Unknown argument: ", arg, call. = FALSE)
  }

  if (opts$smoke) {
    if (!explicit$reps) opts$reps <- 2L
    if (!explicit$n_grid) opts$n_grid <- c(120L, 400L)
    if (!explicit$profiles) opts$profiles <- "p2"
    if (!explicit$prune_reps) opts$prune_reps <- opts$reps
    if (!explicit$prune_n) opts$prune_n <- max(opts$n_grid)
    if (!explicit$prune_tols) opts$prune_tols <- c(1e-9, 1e-6)
  }
  if (is.na(opts$prune_reps)) opts$prune_reps <- opts$reps

  if (any(is.na(opts$n_grid)) || any(opts$n_grid < 2L)) {
    stop("--n-grid must contain integers >= 2.", call. = FALSE)
  }
  if (is.na(opts$reps) || opts$reps < 1L) stop("--reps must be >= 1.", call. = FALSE)
  if (is.na(opts$prune_reps) || opts$prune_reps < 1L) {
    stop("--prune-reps must be >= 1.", call. = FALSE)
  }
  if (is.na(opts$prune_n) || opts$prune_n < 2L) stop("--prune-n must be >= 2.", call. = FALSE)
  if (any(is.na(opts$prune_tols)) || any(opts$prune_tols <= 0)) {
    stop("--prune-tols must contain positive values.", call. = FALSE)
  }
  if (!all(opts$profiles %in% c("p1", "p2"))) {
    stop("--profiles must contain only p1 and/or p2.", call. = FALSE)
  }
  if (is.na(opts$seed_base)) stop("--seed-base must be an integer.", call. = FALSE)
  if (is.na(opts$progress_every) || opts$progress_every < 1L) {
    stop("--progress-every must be >= 1.", call. = FALSE)
  }
  opts
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))

profiles <- list(
  p1 = list(
    profile = "p1",
    label = "P1: C = 5, R = 6",
    C = 5L,
    R = 6L,
    p = colSums(as.matrix(dat.fleiss1971)) / sum(as.matrix(dat.fleiss1971)),
    rho = 0.92,
    pi_rater = c(0.95, 0.90, 0.80, 0.65, 0.50, 0.35)
  ),
  p2 = list(
    profile = "p2",
    label = "P2: C = 3, R = 4",
    C = 3L,
    R = 4L,
    p = c(0.25, 0.50, 0.25),
    rho = 0.92,
    pi_rater = c(0.95, 0.80, 0.60, 0.40)
  )
)
profiles <- profiles[opts$profiles]

coefficients <- c("Conger", "Fleiss", "Brennan-Prediger")
default_prune_tol <- 1e-9
em_base <- list(max_iter = 50000L, tol = 1e-7, start_alpha = 0.1)

loss_identity <- function(C) {
  L <- matrix(1, nrow = C, ncol = C)
  diag(L) <- 0
  L
}

population_truth <- function(profile) {
  C <- profile$C
  R <- profile$R
  p <- profile$p
  rho <- profile$rho
  L <- loss_identity(C)

  cond <- lapply(seq_len(R), function(j) {
    out <- matrix(0, nrow = C, ncol = C)
    for (truth in seq_len(C)) {
      out[truth, ] <- (1 - rho) * p
      out[truth, truth] <- out[truth, truth] + rho
    }
    out
  })
  marg <- lapply(cond, function(mat) as.numeric(colSums(mat * p)))

  loss_between <- function(a, b) as.numeric(t(a) %*% L %*% b)

  d_obs <- 0
  d_conger <- 0
  n_pairs <- R * (R - 1) / 2
  for (j in seq_len(R - 1L)) {
    for (k in (j + 1L):R) {
      for (truth in seq_len(C)) {
        d_obs <- d_obs + p[[truth]] * loss_between(cond[[j]][truth, ], cond[[k]][truth, ])
      }
      d_conger <- d_conger + loss_between(marg[[j]], marg[[k]])
    }
  }
  d_obs <- d_obs / n_pairs
  d_conger <- d_conger / n_pairs

  pooled <- Reduce(`+`, marg) / R
  d_fleiss <- loss_between(pooled, pooled)
  d_bp <- sum(L) / (C * C)

  setNames(
    c(
      1 - d_obs / d_conger,
      1 - d_obs / d_fleiss,
      1 - d_obs / d_bp
    ),
    coefficients
  )
}

simulate_complete <- function(n, profile) {
  C <- profile$C
  R <- profile$R
  truth <- sample.int(C, n, replace = TRUE, prob = profile$p)
  x_star <- matrix(NA_integer_, nrow = n, ncol = R)
  for (j in seq_len(R)) {
    correct <- stats::runif(n) < profile$rho
    guessed <- sample.int(C, n, replace = TRUE, prob = profile$p)
    x_star[, j] <- ifelse(correct, truth, guessed)
  }
  x_star
}

apply_missing_mcar <- function(x_star, pi_rater) {
  x <- x_star
  for (j in seq_len(ncol(x_star))) {
    observed <- stats::runif(nrow(x_star)) < pi_rater[[j]]
    x[!observed, j] <- NA_integer_
  }
  x
}

fit_fiml <- function(x, prune_tol) {
  tryCatch({
    fit <- misskappa::kappa(
      x,
      method = "fiml",
      weight = "identity",
      em_options = c(em_base, list(prune_tol = prune_tol))
    )
    se <- sqrt(diag(vcov(fit)))
    list(estimates = coef(fit), se = se, error = "")
  }, error = function(e) {
    empty <- setNames(rep(NA_real_, length(coefficients)), coefficients)
    list(estimates = empty, se = empty, error = conditionMessage(e))
  })
}

run_one <- function(phase, profile, n, rep_id, prune_tol, truth, seed) {
  set.seed(seed)
  x_star <- simulate_complete(n, profile)
  x <- apply_missing_mcar(x_star, profile$pi_rater)
  fit <- fit_fiml(x, prune_tol)
  rows <- vector("list", length(coefficients))
  for (i in seq_along(coefficients)) {
    coef_name <- coefficients[[i]]
    rows[[i]] <- data.frame(
      phase = phase,
      profile = profile$profile,
      profile_label = profile$label,
      C = profile$C,
      R = profile$R,
      n = n,
      rep = rep_id,
      seed = seed,
      prune_tol = prune_tol,
      coefficient = coef_name,
      estimate = unname(fit$estimates[[coef_name]]),
      se = unname(fit$se[[coef_name]]),
      truth = unname(truth[[coef_name]]),
      error = fit$error,
      observed_fraction = mean(!is.na(x)),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

split_keys <- function(data, keys) interaction(data[, keys], drop = TRUE, lex.order = TRUE)

summarize_replicates <- function(df) {
  keys <- c("phase", "profile", "profile_label", "C", "R", "n", "prune_tol", "coefficient")
  groups <- split(df, split_keys(df, keys))
  out <- lapply(groups, function(g) {
    ok <- is.finite(g$estimate)
    est <- g$estimate[ok]
    bias <- est - g$truth[[1L]]
    data.frame(
      phase = g$phase[[1L]],
      profile = g$profile[[1L]],
      profile_label = g$profile_label[[1L]],
      C = g$C[[1L]],
      R = g$R[[1L]],
      n = g$n[[1L]],
      prune_tol = g$prune_tol[[1L]],
      coefficient = g$coefficient[[1L]],
      reps = length(unique(g$rep)),
      n_valid = sum(ok),
      failures = sum(!ok),
      truth = g$truth[[1L]],
      mean_estimate = if (length(est)) mean(est) else NA_real_,
      bias = if (length(bias)) mean(bias) else NA_real_,
      sd = if (length(est) > 1L) stats::sd(est) else NA_real_,
      mc_se_bias = if (length(est) > 1L) stats::sd(est) / sqrt(length(est)) else NA_real_,
      rmse = if (length(bias)) sqrt(mean(bias^2)) else NA_real_,
      mean_se = if (any(is.finite(g$se))) mean(g$se, na.rm = TRUE) else NA_real_,
      mean_observed_fraction = mean(g$observed_fraction),
      bias_times_n = if (length(bias)) mean(bias) * g$n[[1L]] else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  summary <- do.call(rbind, out)
  summary <- summary[order(
    summary$phase, summary$profile, summary$coefficient,
    summary$prune_tol, summary$n
  ), ]
  rownames(summary) <- NULL
  summary
}

log_progress <- function(...) {
  message(format(Sys.time(), "%H:%M:%S"), " ", sprintf(...))
}

truth_rows <- do.call(rbind, lapply(profiles, function(profile) {
  truth <- population_truth(profile)
  data.frame(
    profile = profile$profile,
    profile_label = profile$label,
    C = profile$C,
    R = profile$R,
    coefficient = names(truth),
    truth = as.numeric(truth),
    stringsAsFactors = FALSE
  )
}))

replicate_rows <- list()
cell <- 0L
total_cells <- length(profiles) * (length(opts$n_grid) + length(opts$prune_tols))
for (profile in profiles) {
  truth <- setNames(
    truth_rows$truth[truth_rows$profile == profile$profile],
    truth_rows$coefficient[truth_rows$profile == profile$profile]
  )
  for (n in opts$n_grid) {
    cell <- cell + 1L
    log_progress("scaling cell %d/%d: %s n=%d", cell, total_cells, profile$profile, n)
    for (rep_id in seq_len(opts$reps)) {
      if (rep_id == 1L || rep_id == opts$reps || rep_id %% opts$progress_every == 0L) {
        log_progress("  replicate %d/%d", rep_id, opts$reps)
      }
      seed <- opts$seed_base + 100000L * cell + rep_id
      replicate_rows[[length(replicate_rows) + 1L]] <- run_one(
        "scaling", profile, n, rep_id, default_prune_tol, truth, seed
      )
    }
  }
  for (prune_tol in opts$prune_tols) {
    cell <- cell + 1L
    log_progress(
      "prune cell %d/%d: %s n=%d prune_tol=%g",
      cell, total_cells, profile$profile, opts$prune_n, prune_tol
    )
    for (rep_id in seq_len(opts$prune_reps)) {
      if (rep_id == 1L || rep_id == opts$prune_reps || rep_id %% opts$progress_every == 0L) {
        log_progress("  replicate %d/%d", rep_id, opts$prune_reps)
      }
      seed <- opts$seed_base + 100000L * cell + rep_id
      replicate_rows[[length(replicate_rows) + 1L]] <- run_one(
        "prune", profile, opts$prune_n, rep_id, prune_tol, truth, seed
      )
    }
  }
}

replicates <- do.call(rbind, replicate_rows)
summary <- summarize_replicates(replicates)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) == 0L) {
  file.path(script_dir, "run_experiment.R")
} else {
  sub("^--file=", "", script_arg[[1L]])
}

metadata <- data.frame(
  key = c(
    "generated_at", "script", "smoke", "reps", "n_grid", "profiles",
    "prune_reps", "prune_n", "prune_tols", "seed_base", "default_prune_tol",
    "em_max_iter", "em_tol", "weight", "misskappa_version", "r_version"
  ),
  value = c(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
    normalizePath(script_path, mustWork = FALSE),
    as.character(opts$smoke),
    as.character(opts$reps),
    paste(opts$n_grid, collapse = ","),
    paste(opts$profiles, collapse = ","),
    as.character(opts$prune_reps),
    as.character(opts$prune_n),
    paste(opts$prune_tols, collapse = ","),
    as.character(opts$seed_base),
    as.character(default_prune_tol),
    as.character(em_base$max_iter),
    as.character(em_base$tol),
    "identity",
    as.character(utils::packageVersion("misskappa")),
    R.version.string
  ),
  stringsAsFactors = FALSE
)

dir.create(opts$out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(metadata, file.path(opts$out_dir, "metadata.csv"), row.names = FALSE)
write.csv(truth_rows, file.path(opts$out_dir, "truth.csv"), row.names = FALSE)
write.csv(replicates, file.path(opts$out_dir, "replicates.csv"), row.names = FALSE)
write.csv(summary, file.path(opts$out_dir, "summary.csv"), row.names = FALSE)

cat("Wrote:\n")
cat(sprintf("  %s\n", file.path(opts$out_dir, "metadata.csv")))
cat(sprintf("  %s\n", file.path(opts$out_dir, "truth.csv")))
cat(sprintf("  %s\n", file.path(opts$out_dir, "replicates.csv")))
cat(sprintf("  %s\n", file.path(opts$out_dir, "summary.csv")))
