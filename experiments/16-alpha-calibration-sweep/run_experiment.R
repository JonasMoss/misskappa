#!/usr/bin/env Rscript
#
# Experiment 16: categorical coefficient-alpha calibration sweep.
# Fixed population alpha is computed by one-dimensional quadrature under the
# latent ordinal item DGP. Simulations compare alpha-available and saturated
# categorical FIML under complete data, MCAR, and an anchor-item MAR mechanism.

suppressPackageStartupMessages({
  library(misskappa)
})

usage <- function(status = 0L) {
  cat(paste0(
    "Usage: Rscript run_experiment.R [options]\n\n",
    "Options:\n",
    "  --help, -h             Show this help and exit.\n",
    "  --smoke                Cheap check: tiny3x4 only, n=120, reps=4.\n",
    "  --reps N               Replicates per cell. Default: 40.\n",
    "  --n-grid LIST          Comma-separated n grid. Default: 250,1000.\n",
    "  --paper-n-grid LIST    Optional n grid for paper5x6 only. Default: --n-grid.\n",
    "  --profiles LIST        Profiles: tiny3x4,ordinal4x5,paper5x6. Default: all.\n",
    "  --mechanisms LIST      Mechanisms: complete,mcar,mar. Default: all three.\n",
    "  --seed-base N          Base seed. Default: 161600.\n",
    "  --out-dir PATH         Output directory. Default: script-local results/.\n",
    "  --progress             Print one line per design cell.\n\n",
    "Examples:\n",
    "  Rscript run_experiment.R --smoke\n",
    "  Rscript run_experiment.R --reps 200 --n-grid 250,1000,4000 --paper-n-grid 250,1000\n"
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

parse_csv_int <- function(x, arg) {
  out <- as.integer(strsplit(x, ",", fixed = TRUE)[[1L]])
  if (any(is.na(out))) stop(arg, " must be a comma-separated integer list.", call. = FALSE)
  unique(out)
}

parse_args <- function(argv) {
  opts <- list(
    smoke = FALSE,
    reps = 40L,
    n_grid = c(250L, 1000L),
    paper_n_grid = NULL,
    profiles = c("tiny3x4", "ordinal4x5", "paper5x6"),
    mechanisms = c("complete", "mcar", "mar"),
    seed_base = 161600L,
    out_dir = file.path(script_dir, "results"),
    progress = FALSE
  )
  explicit <- list(
    reps = FALSE, n_grid = FALSE, paper_n_grid = FALSE,
    profiles = FALSE, mechanisms = FALSE
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

    needs_value <- c("--reps", "--n-grid", "--paper-n-grid", "--profiles", "--mechanisms",
                     "--seed-base", "--out-dir")
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
      if (arg == "--paper-n-grid") {
        opts$paper_n_grid <- parse_csv_int(val, arg)
        explicit$paper_n_grid <- TRUE
      }
      if (arg == "--profiles") {
        opts$profiles <- parse_csv_chr(val)
        explicit$profiles <- TRUE
      }
      if (arg == "--mechanisms") {
        opts$mechanisms <- parse_csv_chr(val)
        explicit$mechanisms <- TRUE
      }
      if (arg == "--seed-base") opts$seed_base <- as.integer(val)
      if (arg == "--out-dir") opts$out_dir <- val
      i <- i + 2L
      next
    }
    stop("Unknown argument: ", arg, call. = FALSE)
  }

  if (opts$smoke) {
    if (!explicit$reps) opts$reps <- 4L
    if (!explicit$n_grid) opts$n_grid <- 120L
    if (!explicit$paper_n_grid) opts$paper_n_grid <- NULL
    if (!explicit$profiles) opts$profiles <- "tiny3x4"
    if (!explicit$mechanisms) opts$mechanisms <- c("complete", "mcar", "mar")
  }

  profile_names <- c("tiny3x4", "ordinal4x5", "paper5x6")
  if (!all(opts$profiles %in% profile_names)) {
    stop("--profiles must contain only: ", paste(profile_names, collapse = ","), call. = FALSE)
  }
  if (!all(opts$mechanisms %in% c("complete", "mcar", "mar"))) {
    stop("--mechanisms must contain only complete, mcar, mar.", call. = FALSE)
  }
  if (is.na(opts$reps) || opts$reps < 1L) stop("--reps must be >= 1.", call. = FALSE)
  if (any(is.na(opts$n_grid)) || any(opts$n_grid < 20L)) {
    stop("--n-grid must contain integers >= 20.", call. = FALSE)
  }
  if (!is.null(opts$paper_n_grid) &&
      (any(is.na(opts$paper_n_grid)) || any(opts$paper_n_grid < 20L))) {
    stop("--paper-n-grid must contain integers >= 20.", call. = FALSE)
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

em_options <- list(
  max_iter = 50000L,
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

category_probs <- function(theta, profile, item) {
  eta <- profile$shifts[[item]] + profile$loadings[[item]] * theta
  cuts <- c(-Inf, profile$thresholds, Inf)
  upper <- (cuts[-1L] - eta) / profile$error_sd
  lower <- (cuts[-length(cuts)] - eta) / profile$error_sd
  stats::pnorm(upper) - stats::pnorm(lower)
}

item_mean_given_theta <- function(theta, profile, item) {
  sum(score_values(profile) * category_probs(theta, profile, item))
}

item_second_given_theta <- function(theta, profile, item) {
  scores <- score_values(profile)
  sum(scores * scores * category_probs(theta, profile, item))
}

normal_expectation <- function(f) {
  integrand <- function(z) {
    vapply(z, function(zz) f(zz), numeric(1)) * stats::dnorm(z)
  }
  stats::integrate(
    integrand,
    lower = -Inf,
    upper = Inf,
    rel.tol = 1e-10,
    subdivisions = 200L
  )$value
}

population_truth <- function(profile) {
  R <- profile$R
  mu <- numeric(R)
  second <- numeric(R)
  exy <- matrix(0, nrow = R, ncol = R)

  for (j in seq_len(R)) {
    mu[j] <- normal_expectation(function(theta) item_mean_given_theta(theta, profile, j))
    second[j] <- normal_expectation(function(theta) item_second_given_theta(theta, profile, j))
    exy[j, j] <- second[j]
  }
  for (j in seq_len(R - 1L)) {
    for (k in (j + 1L):R) {
      exy[j, k] <- normal_expectation(function(theta) {
        item_mean_given_theta(theta, profile, j) *
          item_mean_given_theta(theta, profile, k)
      })
      exy[k, j] <- exy[j, k]
    }
  }

  sigma <- exy - outer(mu, mu)
  t1 <- sum(sigma)
  t2 <- sum(diag(sigma))
  alpha <- (R / (R - 1.0)) * (1.0 - t2 / t1)
  list(alpha = alpha, total_variance = t1, item_variance_sum = t2)
}

simulate_complete <- function(n, profile) {
  theta <- stats::rnorm(n)
  x <- matrix(0L, nrow = n, ncol = profile$R)
  for (j in seq_len(profile$R)) {
    y <- profile$loadings[[j]] * theta + profile$shifts[[j]] +
      stats::rnorm(n, sd = profile$error_sd)
    x[, j] <- as.integer(findInterval(y, profile$thresholds))
  }
  x
}

apply_missing <- function(x_star, mechanism, profile) {
  if (mechanism == "complete") return(x_star)

  x <- x_star
  n <- nrow(x)
  R <- ncol(x)
  C <- profile$C

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
  x
}

min_pair_count <- function(x) {
  R <- ncol(x)
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

observed_category_count <- function(x) length(unique(as.vector(x[!is.na(x)])))

fit_alpha <- function(x, method) {
  start <- proc.time()[["elapsed"]]
  tryCatch({
    fit <- misskappa::alpha(
      x,
      method = method,
      values = NULL,
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
      se = if (is.finite(v) && v >= 0.0) sqrt(v) else NA_real_,
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
}

truth_rows <- do.call(rbind, lapply(profiles, function(profile) {
  truth <- population_truth(profile)
  data.frame(
    profile = profile$profile,
    profile_label = profile$label,
    C = profile$C,
    R = profile$R,
    patterns = pattern_count(profile),
    truth_alpha = truth$alpha,
    total_variance = truth$total_variance,
    item_variance_sum = truth$item_variance_sum,
    stringsAsFactors = FALSE
  )
}))

truth_for_profile <- function(profile_name) {
  truth_rows$truth_alpha[match(profile_name, truth_rows$profile)]
}

run_one <- function(profile, mechanism, n, rep_id, seed) {
  set.seed(seed)
  x_star <- simulate_complete(n, profile)
  x <- apply_missing(x_star, mechanism, profile)
  truth <- truth_for_profile(profile$profile)

  rows <- list()
  for (method in c("available", "fiml")) {
    fit <- fit_alpha(x, method)
    rows[[method]] <- data.frame(
      profile = profile$profile,
      profile_label = profile$label,
      C = profile$C,
      R = profile$R,
      patterns = pattern_count(profile),
      mechanism = mechanism,
      n = n,
      rep = rep_id,
      seed = seed,
      method = method,
      truth = truth,
      estimate = fit$estimate,
      bias = fit$estimate - truth,
      se = fit$se,
      vcov = fit$vcov,
      elapsed_ms = fit$elapsed_ms,
      observed_fraction = mean(!is.na(x)),
      min_pair_count = min_pair_count(x),
      observed_categories = observed_category_count(x),
      psi_vcov_max_abs_error = fit$psi_vcov_max_abs_error,
      error = fit$error,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

split_keys <- function(data, keys) interaction(data[, keys], drop = TRUE, lex.order = TRUE)

summarize_replicates <- function(df) {
  keys <- c("profile", "profile_label", "C", "R", "patterns", "mechanism", "n", "method")
  groups <- split(df, split_keys(df, keys))
  out <- lapply(groups, function(g) {
    ok <- is.finite(g$estimate)
    finite_se <- ok & is.finite(g$se) & g$se > 0
    z <- (g$estimate[finite_se] - g$truth[finite_se]) / g$se[finite_se]
    coverage <- abs(g$estimate[finite_se] - g$truth[finite_se]) <= 1.96 * g$se[finite_se]
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
      n_finite_se = sum(finite_se),
      failures = sum(!ok),
      truth = g$truth[[1L]],
      mean_estimate = if (any(ok)) mean(g$estimate[ok]) else NA_real_,
      bias = if (any(ok)) mean(g$bias[ok]) else NA_real_,
      sd_estimate = if (sum(ok) > 1L) stats::sd(g$estimate[ok]) else NA_real_,
      mc_se_bias = if (sum(ok) > 1L) stats::sd(g$estimate[ok]) / sqrt(sum(ok)) else NA_real_,
      rmse = if (any(ok)) sqrt(mean(g$bias[ok]^2)) else NA_real_,
      mean_se = if (any(finite_se)) mean(g$se[finite_se]) else NA_real_,
      se_over_sd = if (sum(finite_se) > 1L && stats::sd(g$estimate[finite_se]) > 1e-10) {
        mean(g$se[finite_se]) / stats::sd(g$estimate[finite_se])
      } else {
        NA_real_
      },
      coverage95 = if (any(finite_se)) mean(coverage) else NA_real_,
      mean_z = if (length(z)) mean(z) else NA_real_,
      sd_z = if (length(z) > 1L) stats::sd(z) else NA_real_,
      mean_elapsed_ms = mean(g$elapsed_ms, na.rm = TRUE),
      median_elapsed_ms = stats::median(g$elapsed_ms, na.rm = TRUE),
      mean_observed_fraction = mean(g$observed_fraction, na.rm = TRUE),
      min_pair_count_min = min(g$min_pair_count, na.rm = TRUE),
      min_observed_categories = min(g$observed_categories, na.rm = TRUE),
      max_psi_vcov_abs_error = if (length(psi_errors)) max(psi_errors) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  summary <- do.call(rbind, out)
  summary <- summary[order(summary$profile, summary$mechanism, summary$n, summary$method), ]
  rownames(summary) <- NULL
  summary
}

replicate_rows <- list()
cell <- 0L
n_grid_for_profile <- function(profile) {
  if (identical(profile$profile, "paper5x6") && !is.null(opts$paper_n_grid)) {
    opts$paper_n_grid
  } else {
    opts$n_grid
  }
}

total_cells <- sum(vapply(profiles, function(profile) {
  length(opts$mechanisms) * length(n_grid_for_profile(profile))
}, integer(1)))
for (profile in profiles) {
  for (mechanism in opts$mechanisms) {
    for (n in n_grid_for_profile(profile)) {
      cell <- cell + 1L
      log_progress(
        "cell %d/%d: %s %s n=%d reps=%d",
        cell, total_cells, profile$profile, mechanism, n, opts$reps
      )
      for (rep_id in seq_len(opts$reps)) {
        seed <- opts$seed_base + 1000000L * cell + rep_id
        replicate_rows[[length(replicate_rows) + 1L]] <- run_one(
          profile, mechanism, n, rep_id, seed
        )
      }
    }
  }
}

replicates <- do.call(rbind, replicate_rows)
summary <- summarize_replicates(replicates)

dir.create(opts$out_dir, recursive = TRUE, showWarnings = FALSE)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) == 0L) {
  file.path(script_dir, "run_experiment.R")
} else {
  sub("^--file=", "", script_arg[[1L]])
}

metadata <- data.frame(
  key = c(
    "generated_at", "script", "smoke", "reps", "n_grid", "paper_n_grid",
    "profiles", "mechanisms", "seed_base", "em_max_iter", "em_tol",
    "em_prune_tol", "em_start_alpha", "em_info_rcond", "truth_method",
    "misskappa_version", "r_version"
  ),
  value = c(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
    normalizePath(script_path, mustWork = FALSE),
    as.character(opts$smoke),
    as.character(opts$reps),
    paste(opts$n_grid, collapse = ","),
    if (is.null(opts$paper_n_grid)) "" else paste(opts$paper_n_grid, collapse = ","),
    paste(opts$profiles, collapse = ","),
    paste(opts$mechanisms, collapse = ","),
    as.character(opts$seed_base),
    as.character(em_options$max_iter),
    as.character(em_options$tol),
    as.character(em_options$prune_tol),
    as.character(em_options$start_alpha),
    as.character(em_options$info_rcond),
    "one-dimensional quadrature over latent normal factor",
    as.character(utils::packageVersion("misskappa")),
    R.version.string
  ),
  stringsAsFactors = FALSE
)

write.csv(truth_rows, file.path(opts$out_dir, "truth.csv"), row.names = FALSE)
write.csv(replicates, file.path(opts$out_dir, "replicates.csv"), row.names = FALSE)
write.csv(summary, file.path(opts$out_dir, "summary.csv"), row.names = FALSE)
write.csv(metadata, file.path(opts$out_dir, "metadata.csv"), row.names = FALSE)

cat("Wrote:\n")
cat("  ", normalizePath(file.path(opts$out_dir, "truth.csv"), mustWork = FALSE), "\n", sep = "")
cat("  ", normalizePath(file.path(opts$out_dir, "replicates.csv"), mustWork = FALSE), "\n", sep = "")
cat("  ", normalizePath(file.path(opts$out_dir, "summary.csv"), mustWork = FALSE), "\n", sep = "")
cat("  ", normalizePath(file.path(opts$out_dir, "metadata.csv"), mustWork = FALSE), "\n", sep = "")
