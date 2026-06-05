#!/usr/bin/env Rscript
#
# Experiment 15: categorical coefficient-alpha smoke checks.
# The goal is mechanical validation, not paper-grade coverage. The runner
# compares the new alpha-available and alpha-FIML paths against each
# replicate's complete-data alpha, then records small C^R timing cells for the
# saturated categorical EM.

suppressPackageStartupMessages({
  library(misskappa)
})

usage <- function(status = 0L) {
  cat(paste0(
    "Usage: Rscript run_experiment.R [options]\n\n",
    "Options:\n",
    "  --help, -h             Show this help and exit.\n",
    "  --smoke                Cheap check: tiny profile, n=160, reps=5.\n",
    "  --reps N               Replicates per mechanism/profile cell. Default: 40.\n",
    "  --n N                  Sample size for estimate checks. Default: 600.\n",
    "  --profiles LIST        Comma-separated profiles. Default: tiny3x4,ordinal4x5.\n",
    "  --mechanisms LIST      Mechanisms: complete,mcar,mar. Default: all three.\n",
    "  --timing-profiles LIST C^R timing profiles. Default: tiny3x4,ordinal4x5,paper5x6.\n",
    "  --timing-n N           Sample size for timing cells. Default: 300.\n",
    "  --timing-reps N        Replicates per timing profile. Default: 2.\n",
    "  --seed-base N          Base seed. Default: 151500.\n",
    "  --out-dir PATH         Output directory. Default: script-local results/.\n",
    "  --progress             Print one line per design cell.\n\n",
    "Examples:\n",
    "  Rscript run_experiment.R --smoke\n",
    "  Rscript run_experiment.R --reps 100 --n 1000 --timing-reps 4\n"
  ))
  quit(save = "no", status = status)
}

script_dir <- local({
  cmd <- commandArgs(FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) == 0L) getwd()
  else dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE))
})

parse_csv_chr <- function(x) trimws(strsplit(x, ",", fixed = TRUE)[[1L]])

parse_args <- function(argv) {
  opts <- list(
    smoke = FALSE,
    reps = 40L,
    n = 600L,
    profiles = c("tiny3x4", "ordinal4x5"),
    mechanisms = c("complete", "mcar", "mar"),
    timing_profiles = c("tiny3x4", "ordinal4x5", "paper5x6"),
    timing_n = 300L,
    timing_reps = 2L,
    seed_base = 151500L,
    out_dir = file.path(script_dir, "results"),
    progress = FALSE
  )

  explicit <- list(
    reps = FALSE, n = FALSE, profiles = FALSE, mechanisms = FALSE,
    timing_profiles = FALSE, timing_n = FALSE, timing_reps = FALSE
  )

  i <- 1L
  while (i <= length(argv)) {
    arg <- argv[[i]]
    if (arg %in% c("--help", "-h")) usage(0L)
    if (arg == "--smoke") {
      opts$smoke <- TRUE
      i <- i + 1L
      next
    }
    if (arg == "--progress") {
      opts$progress <- TRUE
      i <- i + 1L
      next
    }

    needs_value <- c(
      "--reps", "--n", "--profiles", "--mechanisms", "--timing-profiles",
      "--timing-n", "--timing-reps", "--seed-base", "--out-dir"
    )
    if (arg %in% needs_value) {
      if (i == length(argv)) stop(arg, " needs a value.", call. = FALSE)
      val <- argv[[i + 1L]]
      if (arg == "--reps") {
        opts$reps <- as.integer(val)
        explicit$reps <- TRUE
      }
      if (arg == "--n") {
        opts$n <- as.integer(val)
        explicit$n <- TRUE
      }
      if (arg == "--profiles") {
        opts$profiles <- parse_csv_chr(val)
        explicit$profiles <- TRUE
      }
      if (arg == "--mechanisms") {
        opts$mechanisms <- parse_csv_chr(val)
        explicit$mechanisms <- TRUE
      }
      if (arg == "--timing-profiles") {
        opts$timing_profiles <- parse_csv_chr(val)
        explicit$timing_profiles <- TRUE
      }
      if (arg == "--timing-n") {
        opts$timing_n <- as.integer(val)
        explicit$timing_n <- TRUE
      }
      if (arg == "--timing-reps") {
        opts$timing_reps <- as.integer(val)
        explicit$timing_reps <- TRUE
      }
      if (arg == "--seed-base") opts$seed_base <- as.integer(val)
      if (arg == "--out-dir") opts$out_dir <- val
      i <- i + 2L
      next
    }
    stop("Unknown argument: ", arg, call. = FALSE)
  }

  if (opts$smoke) {
    if (!explicit$reps) opts$reps <- 5L
    if (!explicit$n) opts$n <- 160L
    if (!explicit$profiles) opts$profiles <- "tiny3x4"
    if (!explicit$mechanisms) opts$mechanisms <- c("complete", "mcar", "mar")
    if (!explicit$timing_profiles) opts$timing_profiles <- c("tiny3x4", "ordinal4x5")
    if (!explicit$timing_n) opts$timing_n <- 120L
    if (!explicit$timing_reps) opts$timing_reps <- 1L
  }

  profile_names <- c("tiny3x4", "ordinal4x5", "paper5x6")
  if (!all(opts$profiles %in% profile_names)) {
    stop("--profiles must contain only: ", paste(profile_names, collapse = ","), call. = FALSE)
  }
  if (!all(opts$timing_profiles %in% profile_names)) {
    stop("--timing-profiles must contain only: ", paste(profile_names, collapse = ","), call. = FALSE)
  }
  if (!all(opts$mechanisms %in% c("complete", "mcar", "mar"))) {
    stop("--mechanisms must contain only complete, mcar, mar.", call. = FALSE)
  }
  if (is.na(opts$reps) || opts$reps < 1L) stop("--reps must be >= 1.", call. = FALSE)
  if (is.na(opts$n) || opts$n < 20L) stop("--n must be >= 20.", call. = FALSE)
  if (is.na(opts$timing_n) || opts$timing_n < 20L) {
    stop("--timing-n must be >= 20.", call. = FALSE)
  }
  if (is.na(opts$timing_reps) || opts$timing_reps < 1L) {
    stop("--timing-reps must be >= 1.", call. = FALSE)
  }
  if (is.na(opts$seed_base)) stop("--seed-base must be an integer.", call. = FALSE)
  opts
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))

profiles_all <- list(
  tiny3x4 = list(
    profile = "tiny3x4",
    label = "C = 3, R = 4",
    C = 3L,
    R = 4L,
    loadings = c(0.85, 1.00, 0.90, 1.10),
    shifts = c(-0.20, 0.00, 0.15, -0.05),
    thresholds = c(-0.55, 0.55),
    error_sd = 0.70,
    mcar_obs = 0.78,
    mar_intercept = 0.35,
    mar_slope = -1.80
  ),
  ordinal4x5 = list(
    profile = "ordinal4x5",
    label = "C = 4, R = 5",
    C = 4L,
    R = 5L,
    loadings = c(0.75, 0.90, 1.05, 0.85, 1.10),
    shifts = c(-0.20, 0.00, 0.10, -0.10, 0.15),
    thresholds = c(-0.85, -0.05, 0.80),
    error_sd = 0.75,
    mcar_obs = 0.78,
    mar_intercept = 0.35,
    mar_slope = -1.80
  ),
  paper5x6 = list(
    profile = "paper5x6",
    label = "C = 5, R = 6",
    C = 5L,
    R = 6L,
    loadings = c(0.70, 0.85, 1.00, 0.90, 1.05, 0.80),
    shifts = c(-0.20, -0.05, 0.10, 0.00, 0.15, -0.10),
    thresholds = c(-1.05, -0.40, 0.35, 1.00),
    error_sd = 0.80,
    mcar_obs = 0.78,
    mar_intercept = 0.35,
    mar_slope = -1.80
  )
)

profiles <- profiles_all[opts$profiles]
timing_profiles <- profiles_all[opts$timing_profiles]

em_options <- list(
  max_iter = 20000L,
  tol = 1e-7,
  prune_tol = 1e-9,
  start_alpha = 0.1,
  info_rcond = 5e-5
)

score_values <- function(profile) seq.int(0, profile$C - 1L)

pattern_count <- function(profile) profile$C ^ profile$R

log_progress <- function(...) {
  if (opts$progress) message(format(Sys.time(), "%H:%M:%S"), " ", sprintf(...))
}

simulate_complete <- function(n, profile) {
  theta <- stats::rnorm(n)
  x <- matrix(0L, nrow = n, ncol = profile$R)
  for (j in seq_len(profile$R)) {
    y <- profile$loadings[[j]] * theta + profile$shifts[[j]] +
      stats::rnorm(n, sd = profile$error_sd)
    x[, j] <- as.integer(findInterval(y, profile$thresholds))
  }

  # Stabilise category support for small smoke runs. The adjustment is tiny
  # relative to n and avoids comparing estimators on silently different C.
  for (k in seq_len(profile$C)) x[k, 1L] <- k - 1L
  x
}

apply_missing <- function(x_star, mechanism, profile) {
  if (mechanism == "complete") return(x_star)

  x <- x_star
  n <- nrow(x)
  R <- ncol(x)
  C <- profile$C
  support_rows <- seq_len(C)

  if (mechanism == "mcar") {
    observed <- matrix(stats::runif(n * R) < profile$mcar_obs, nrow = n, ncol = R)
    x[!observed] <- NA_integer_
  } else if (mechanism == "mar") {
    anchor <- x_star[, 1L] / max(1L, C - 1L)
    x[, 1L] <- x_star[, 1L]
    for (j in 2:R) {
      item_shift <- (j - (R + 1) / 2) * 0.08
      p_obs <- stats::plogis(
        profile$mar_intercept + item_shift +
          profile$mar_slope * (anchor - 0.5)
      )
      observed <- stats::runif(n) < p_obs
      x[!observed, j] <- NA_integer_
    }
  } else {
    stop("Unknown mechanism: ", mechanism, call. = FALSE)
  }

  x[support_rows, 1L] <- seq.int(0L, C - 1L)
  x
}

min_pair_count <- function(x) {
  R <- ncol(x)
  if (R < 2L) return(NA_integer_)
  counts <- integer(R * (R - 1L) / 2L)
  pos <- 1L
  for (j in seq_len(R - 1L)) {
    for (k in (j + 1L):R) {
      counts[pos] <- sum(!is.na(x[, j]) & !is.na(x[, k]))
      pos <- pos + 1L
    }
  }
  min(counts)
}

fit_alpha <- function(x, method, values) {
  start <- proc.time()[["elapsed"]]
  out <- tryCatch({
    fit <- misskappa::alpha(
      x,
      method = method,
      values = values,
      em_options = if (method == "fiml") em_options else list()
    )
    elapsed_ms <- 1000 * (proc.time()[["elapsed"]] - start)
    vc <- stats::vcov(fit)
    psi <- stats::influence(fit)
    psi_vcov_max_abs_error <- NA_real_
    if (length(psi) > 0L && nrow(psi) > 0L) {
      psi_vc <- crossprod(psi) / (nrow(psi)^2)
      psi_vcov_max_abs_error <- max(abs(psi_vc - vc))
    }
    v <- as.numeric(vc[1L, 1L])
    list(
      estimate = as.numeric(stats::coef(fit)[["alpha"]]),
      vcov = v,
      se = if (is.finite(v) && v >= 0) sqrt(v) else NA_real_,
      elapsed_ms = elapsed_ms,
      psi_vcov_max_abs_error = psi_vcov_max_abs_error,
      error = ""
    )
  }, error = function(e) {
    elapsed_ms <- 1000 * (proc.time()[["elapsed"]] - start)
    list(
      estimate = NA_real_,
      vcov = NA_real_,
      se = NA_real_,
      elapsed_ms = elapsed_ms,
      psi_vcov_max_abs_error = NA_real_,
      error = conditionMessage(e)
    )
  })
  out
}

run_estimate_cell <- function(profile, mechanism, rep_id, seed) {
  set.seed(seed)
  x_star <- simulate_complete(opts$n, profile)
  values <- score_values(profile)
  ref_fit <- fit_alpha(x_star, "available", values)
  x_obs <- apply_missing(x_star, mechanism, profile)
  observed_fraction <- mean(!is.na(x_obs))
  min_pair <- min_pair_count(x_obs)

  rows <- list()
  for (method in c("available", "fiml")) {
    fit <- fit_alpha(x_obs, method, values)
    rows[[method]] <- data.frame(
      profile = profile$profile,
      profile_label = profile$label,
      C = profile$C,
      R = profile$R,
      patterns = pattern_count(profile),
      mechanism = mechanism,
      n = opts$n,
      rep = rep_id,
      seed = seed,
      method = method,
      estimate = fit$estimate,
      se = fit$se,
      vcov = fit$vcov,
      complete_reference = ref_fit$estimate,
      delta_complete = fit$estimate - ref_fit$estimate,
      elapsed_ms = fit$elapsed_ms,
      observed_fraction = observed_fraction,
      min_pair_count = min_pair,
      psi_vcov_max_abs_error = fit$psi_vcov_max_abs_error,
      error = fit$error,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

run_timing_cell <- function(profile, rep_id, seed) {
  set.seed(seed)
  x_star <- simulate_complete(opts$timing_n, profile)
  x_obs <- apply_missing(x_star, "mar", profile)
  fit <- fit_alpha(x_obs, "fiml", score_values(profile))
  data.frame(
    profile = profile$profile,
    profile_label = profile$label,
    C = profile$C,
    R = profile$R,
    patterns = pattern_count(profile),
    n = opts$timing_n,
    rep = rep_id,
    seed = seed,
    method = "fiml",
    mechanism = "mar",
    estimate = fit$estimate,
    se = fit$se,
    elapsed_ms = fit$elapsed_ms,
    observed_fraction = mean(!is.na(x_obs)),
    min_pair_count = min_pair_count(x_obs),
    error = fit$error,
    stringsAsFactors = FALSE
  )
}

split_keys <- function(data, keys) interaction(data[, keys], drop = TRUE, lex.order = TRUE)

summarize_estimates <- function(df) {
  keys <- c("profile", "profile_label", "C", "R", "patterns", "mechanism", "n", "method")
  groups <- split(df, split_keys(df, keys))
  out <- lapply(groups, function(g) {
    ok <- is.finite(g$estimate)
    finite_se <- ok & is.finite(g$se) & g$se > 0
    delta <- g$delta_complete[ok]
    cover <- finite_se & abs(g$estimate - g$complete_reference) <= 1.96 * g$se
    delta_sd_for_se <- if (sum(finite_se) > 1L) stats::sd(g$delta_complete[finite_se]) else NA_real_
    psi_errors <- g$psi_vcov_max_abs_error[is.finite(g$psi_vcov_max_abs_error)]
    data.frame(
      profile = g$profile[[1L]],
      profile_label = g$profile_label[[1L]],
      C = g$C[[1L]],
      R = g$R[[1L]],
      patterns = g$patterns[[1L]],
      mechanism = g$mechanism[[1L]],
      n = g$n[[1L]],
      method = g$method[[1L]],
      reps = length(unique(g$rep)),
      n_valid = sum(ok),
      failures = sum(!ok),
      complete_reference_mean = mean(g$complete_reference, na.rm = TRUE),
      mean_estimate = if (any(ok)) mean(g$estimate[ok]) else NA_real_,
      mean_delta_complete = if (length(delta)) mean(delta) else NA_real_,
      sd_delta_complete = if (length(delta) > 1L) stats::sd(delta) else NA_real_,
      mc_se_delta = if (length(delta) > 1L) stats::sd(delta) / sqrt(length(delta)) else NA_real_,
      rmse_delta_complete = if (length(delta)) sqrt(mean(delta^2)) else NA_real_,
      mean_se = if (any(finite_se)) mean(g$se[finite_se]) else NA_real_,
      se_over_sd_delta = if (is.finite(delta_sd_for_se) && delta_sd_for_se > 1e-10) {
        mean(g$se[finite_se]) / delta_sd_for_se
      } else {
        NA_real_
      },
      cover_complete_ref95 = if (any(finite_se)) mean(cover[finite_se]) else NA_real_,
      mean_elapsed_ms = mean(g$elapsed_ms, na.rm = TRUE),
      median_elapsed_ms = stats::median(g$elapsed_ms, na.rm = TRUE),
      mean_observed_fraction = mean(g$observed_fraction, na.rm = TRUE),
      min_pair_count_min = min(g$min_pair_count, na.rm = TRUE),
      max_psi_vcov_abs_error = if (length(psi_errors)) max(psi_errors) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  summary <- do.call(rbind, out)
  summary <- summary[order(summary$profile, summary$mechanism, summary$method), ]
  rownames(summary) <- NULL
  summary
}

summarize_timing <- function(df) {
  keys <- c("profile", "profile_label", "C", "R", "patterns", "n", "method", "mechanism")
  groups <- split(df, split_keys(df, keys))
  out <- lapply(groups, function(g) {
    ok <- is.finite(g$estimate)
    data.frame(
      profile = g$profile[[1L]],
      profile_label = g$profile_label[[1L]],
      C = g$C[[1L]],
      R = g$R[[1L]],
      patterns = g$patterns[[1L]],
      n = g$n[[1L]],
      method = g$method[[1L]],
      mechanism = g$mechanism[[1L]],
      reps = length(unique(g$rep)),
      successes = sum(ok),
      failures = sum(!ok),
      mean_elapsed_ms = mean(g$elapsed_ms, na.rm = TRUE),
      median_elapsed_ms = stats::median(g$elapsed_ms, na.rm = TRUE),
      max_elapsed_ms = max(g$elapsed_ms, na.rm = TRUE),
      mean_observed_fraction = mean(g$observed_fraction, na.rm = TRUE),
      min_pair_count_min = min(g$min_pair_count, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  summary <- do.call(rbind, out)
  summary <- summary[order(summary$patterns, summary$n), ]
  rownames(summary) <- NULL
  summary
}

estimate_rows <- list()
cell <- 0L
total_cells <- length(profiles) * length(opts$mechanisms)
for (profile in profiles) {
  for (mechanism in opts$mechanisms) {
    cell <- cell + 1L
    log_progress(
      "estimate cell %d/%d: %s %s n=%d reps=%d",
      cell, total_cells, profile$profile, mechanism, opts$n, opts$reps
    )
    for (rep_id in seq_len(opts$reps)) {
      seed <- opts$seed_base + 100000L * cell + rep_id
      estimate_rows[[length(estimate_rows) + 1L]] <-
        run_estimate_cell(profile, mechanism, rep_id, seed)
    }
  }
}

timing_rows <- list()
timing_cell <- 0L
for (profile in timing_profiles) {
  timing_cell <- timing_cell + 1L
  log_progress(
    "timing cell %d/%d: %s n=%d reps=%d",
    timing_cell, length(timing_profiles), profile$profile, opts$timing_n, opts$timing_reps
  )
  for (rep_id in seq_len(opts$timing_reps)) {
    seed <- opts$seed_base + 90000000L + 100000L * timing_cell + rep_id
    timing_rows[[length(timing_rows) + 1L]] <- run_timing_cell(profile, rep_id, seed)
  }
}

estimates <- do.call(rbind, estimate_rows)
summary <- summarize_estimates(estimates)
timing <- do.call(rbind, timing_rows)
timing_summary <- summarize_timing(timing)

dir.create(opts$out_dir, recursive = TRUE, showWarnings = FALSE)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) == 0L) {
  file.path(script_dir, "run_experiment.R")
} else {
  sub("^--file=", "", script_arg[[1L]])
}

metadata <- data.frame(
  key = c(
    "generated_at", "script", "smoke", "reps", "n", "profiles",
    "mechanisms", "timing_profiles", "timing_n", "timing_reps",
    "seed_base", "em_max_iter", "em_tol", "em_prune_tol",
    "em_start_alpha", "em_info_rcond", "misskappa_version", "r_version"
  ),
  value = c(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
    normalizePath(script_path, mustWork = FALSE),
    as.character(opts$smoke),
    as.character(opts$reps),
    as.character(opts$n),
    paste(opts$profiles, collapse = ","),
    paste(opts$mechanisms, collapse = ","),
    paste(opts$timing_profiles, collapse = ","),
    as.character(opts$timing_n),
    as.character(opts$timing_reps),
    as.character(opts$seed_base),
    as.character(em_options$max_iter),
    as.character(em_options$tol),
    as.character(em_options$prune_tol),
    as.character(em_options$start_alpha),
    as.character(em_options$info_rcond),
    as.character(utils::packageVersion("misskappa")),
    R.version.string
  ),
  stringsAsFactors = FALSE
)

write.csv(estimates, file.path(opts$out_dir, "estimates.csv"), row.names = FALSE)
write.csv(summary, file.path(opts$out_dir, "summary.csv"), row.names = FALSE)
write.csv(timing, file.path(opts$out_dir, "timing.csv"), row.names = FALSE)
write.csv(timing_summary, file.path(opts$out_dir, "timing_summary.csv"), row.names = FALSE)
write.csv(metadata, file.path(opts$out_dir, "metadata.csv"), row.names = FALSE)

cat("Wrote:\n")
cat("  ", normalizePath(file.path(opts$out_dir, "estimates.csv"), mustWork = FALSE), "\n", sep = "")
cat("  ", normalizePath(file.path(opts$out_dir, "summary.csv"), mustWork = FALSE), "\n", sep = "")
cat("  ", normalizePath(file.path(opts$out_dir, "timing.csv"), mustWork = FALSE), "\n", sep = "")
cat("  ", normalizePath(file.path(opts$out_dir, "timing_summary.csv"), mustWork = FALSE), "\n", sep = "")
cat("  ", normalizePath(file.path(opts$out_dir, "metadata.csv"), mustWork = FALSE), "\n", sep = "")
