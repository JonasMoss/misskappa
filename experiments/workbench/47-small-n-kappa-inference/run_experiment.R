#!/usr/bin/env Rscript

# Experiment 47: small-n interval methods for complete-data quadratic Fleiss
# kappa in one deliberately difficult high-agreement regime.

suppressPackageStartupMessages({
  library(misskappa)
})

usage <- function(status = 0L) {
  cat(paste0(
    "Usage: Rscript run_experiment.R [options]\n\n",
    "Options:\n",
    "  --smoke            Mechanical check: n=20, reps=4.\n",
    "  --reps N           Replications per n. Default: 2000.\n",
    "  --n-grid CSV       Sample sizes. Default: 20,30.\n",
    "  --raters N         Number of raters. Default: 5.\n",
    "  --categories N     Ordered categories. Default: 5.\n",
    "  --target K         True quadratic Fleiss kappa. Default: 0.90.\n",
    "  --seed-base N      Deterministic seed base. Default: 474700.\n",
    "  --cores N          Fork workers on non-Windows systems. Default: 1.\n",
    "  --out-dir PATH     Output directory. Default: script-local results/.\n",
    "  --progress         Print one line per sample-size cell.\n",
    "  --help, -h         Show this help and exit.\n\n",
    "Examples:\n",
    "  Rscript run_experiment.R --smoke --progress\n",
    "  Rscript run_experiment.R --reps 2000 --cores 4 --progress\n"
  ))
  quit(save = "no", status = status)
}

script_dir <- local({
  cmd <- commandArgs(FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) == 0L) getwd()
  else dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE))
})
source(file.path(script_dir, "R", "interval_methods.R"), local = TRUE)

parse_int_csv <- function(x, flag) {
  out <- as.integer(strsplit(x, ",", fixed = TRUE)[[1L]])
  if (any(is.na(out))) stop(flag, " must be a comma-separated integer list.",
                            call. = FALSE)
  unique(out)
}

parse_args <- function(argv) {
  opts <- list(smoke = FALSE, reps = 2000L, n_grid = c(20L, 30L),
               raters = 5L, categories = 5L, target = 0.90,
               seed_base = 474700L, cores = 1L,
               out_dir = file.path(script_dir, "results"), progress = FALSE)
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
    needs_value <- c("--reps", "--n-grid", "--raters", "--categories",
                     "--target", "--seed-base", "--cores", "--out-dir")
    if (arg %in% needs_value) {
      if (i == length(argv)) stop(arg, " needs a value.", call. = FALSE)
      value <- argv[[i + 1L]]
      if (arg == "--reps") opts$reps <- as.integer(value)
      else if (arg == "--n-grid") opts$n_grid <- parse_int_csv(value, arg)
      else if (arg == "--raters") opts$raters <- as.integer(value)
      else if (arg == "--categories") opts$categories <- as.integer(value)
      else if (arg == "--target") opts$target <- as.numeric(value)
      else if (arg == "--seed-base") opts$seed_base <- as.integer(value)
      else if (arg == "--cores") opts$cores <- as.integer(value)
      else if (arg == "--out-dir") opts$out_dir <- value
      i <- i + 2L
      next
    }
    stop("Unknown argument: ", arg, call. = FALSE)
  }
  if (opts$smoke) {
    opts$reps <- 4L
    opts$n_grid <- 20L
    opts$cores <- 1L
  }
  if (is.na(opts$reps) || opts$reps < 1L) stop("--reps must be >= 1.", call. = FALSE)
  if (any(is.na(opts$n_grid)) || any(opts$n_grid < 8L)) {
    stop("--n-grid entries must be >= 8.", call. = FALSE)
  }
  if (is.na(opts$raters) || opts$raters < 2L) stop("--raters must be >= 2.", call. = FALSE)
  if (is.na(opts$categories) || opts$categories < 2L) {
    stop("--categories must be >= 2.", call. = FALSE)
  }
  if (is.na(opts$target) || opts$target <= 0 || opts$target >= 1) {
    stop("--target must be in (0, 1).", call. = FALSE)
  }
  if (is.na(opts$seed_base)) stop("--seed-base must be an integer.", call. = FALSE)
  if (is.na(opts$cores) || opts$cores < 1L) stop("--cores must be >= 1.", call. = FALSE)
  opts
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))
values <- seq_len(opts$categories)
knowledge_probability <- sqrt(opts$target)
method_labels <- c(
  wald_z = "Wald (normal)",
  wald_fisher = "Wald (Fisher transform)",
  profile_lr = "Normal profile LR",
  profile_lr_robust = "Robust-scaled profile LR",
  score_normal = "Normal score",
  score_robust = "Robust score",
  score_robust_t = "Robust score, t(n-3)",
  moment_profile = "Robust moment profile",
  moment_profile_f = "Robust moment profile, F(1,n-3)"
)

simulate_perreault_leigh <- function(n, seed) {
  set.seed(seed)
  truth <- sample.int(opts$categories, n, replace = TRUE)
  x <- matrix(NA_integer_, nrow = n, ncol = opts$raters)
  for (j in seq_len(opts$raters)) {
    knows <- stats::runif(n) < knowledge_probability
    guesses <- sample.int(opts$categories, n, replace = TRUE)
    x[, j] <- ifelse(knows, truth, guesses)
  }
  x
}

one_replicate <- function(n, rep) {
  seed <- opts$seed_base + 100000L * n + rep
  started <- proc.time()[["elapsed"]]
  tryCatch({
    x <- simulate_perreault_leigh(n, seed)
    counts <- ratings_to_counts(x, categories = values)
    fit <- kappa_counts(counts, estimator = "fleiss_cuzick",
                        weight = "quadratic", values = values)
    estimate <- unname(stats::coef(fit)[["Fleiss"]])
    se <- sqrt(unname(stats::vcov(fit)["Fleiss", "Fleiss"]))
    focused <- .focused_intervals(x, estimate, se)
    intervals <- focused$intervals
    data.frame(
      n = n, rep = rep, seed = seed, method = rownames(intervals),
      estimate = estimate, se = se,
      lower = intervals[, "lower"], upper = intervals[, "upper"],
      truth = opts$target,
      package_delta_estimate = focused$package_delta_estimate,
      package_delta_se = focused$package_delta_se,
      cs_center_scale = focused$cs_center_scale,
      gmm_rank = focused$gmm_rank,
      gmm_min_eigenvalue = focused$gmm_min_eigenvalue,
      gmm_max_eigenvalue = focused$gmm_max_eigenvalue,
      elapsed_seconds = proc.time()[["elapsed"]] - started,
      status = "ok", stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      n = n, rep = rep, seed = seed, method = names(method_labels),
      estimate = NA_real_, se = NA_real_, lower = NA_real_, upper = NA_real_,
      truth = opts$target, package_delta_estimate = NA_real_,
      package_delta_se = NA_real_, cs_center_scale = NA_real_,
      gmm_rank = NA_integer_, gmm_min_eigenvalue = NA_real_,
      gmm_max_eigenvalue = NA_real_,
      elapsed_seconds = proc.time()[["elapsed"]] - started,
      status = conditionMessage(e), stringsAsFactors = FALSE
    )
  })
}

run_cell <- function(n) {
  if (opts$progress) message(sprintf("n=%d, reps=%d", n, opts$reps))
  reps <- seq_len(opts$reps)
  if (opts$cores > 1L && .Platform$OS.type != "windows") {
    pieces <- parallel::mclapply(reps, function(rep) one_replicate(n, rep),
                                 mc.cores = opts$cores, mc.preschedule = TRUE)
  } else {
    pieces <- lapply(reps, function(rep) one_replicate(n, rep))
  }
  do.call(rbind, pieces)
}

started_all <- proc.time()[["elapsed"]]
replicates <- do.call(rbind, lapply(opts$n_grid, run_cell))
elapsed_all <- proc.time()[["elapsed"]] - started_all
replicates$method_label <- unname(method_labels[replicates$method])
replicates$valid <- replicates$status == "ok" &
  is.finite(replicates$lower) & is.finite(replicates$upper)
replicates$covered <- replicates$valid &
  replicates$lower <= replicates$truth & replicates$truth <= replicates$upper
replicates$miss_below <- replicates$valid & replicates$upper < replicates$truth
replicates$miss_above <- replicates$valid & replicates$lower > replicates$truth
replicates$length <- replicates$upper - replicates$lower

summary_pieces <- split(replicates, interaction(replicates$n, replicates$method,
                                                 drop = TRUE))
summary <- do.call(rbind, lapply(summary_pieces, function(d) {
  valid <- d$valid
  estimates <- d$estimate[d$status == "ok"]
  ses <- d$se[d$status == "ok"]
  data.frame(
    n = d$n[[1L]], method = d$method[[1L]],
    method_label = d$method_label[[1L]], reps = nrow(d),
    valid_reps = sum(valid), failures = sum(!valid),
    coverage = mean(d$covered),
    coverage_mcse = sqrt(mean(d$covered) * (1 - mean(d$covered)) / nrow(d)),
    coverage_valid = if (any(valid)) mean(d$covered[valid]) else NA_real_,
    miss_below = if (any(valid)) mean(d$miss_below[valid]) else NA_real_,
    miss_above = if (any(valid)) mean(d$miss_above[valid]) else NA_real_,
    mean_length = if (any(valid)) mean(d$length[valid]) else NA_real_,
    median_length = if (any(valid)) stats::median(d$length[valid]) else NA_real_,
    bias = mean(estimates, na.rm = TRUE) - opts$target,
    empirical_sd = stats::sd(estimates, na.rm = TRUE),
    mean_se = mean(ses, na.rm = TRUE),
    se_ratio = mean(ses, na.rm = TRUE) / stats::sd(estimates, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
summary <- summary[order(summary$n, match(summary$method, names(method_labels))), ]
rownames(summary) <- NULL

diagnostics <- do.call(rbind, lapply(split(replicates, replicates$n), function(d) {
  one <- d[d$method == "wald_z" & d$status == "ok", ]
  data.frame(
    n = d$n[[1L]],
    max_abs_package_delta_estimate = max(abs(one$package_delta_estimate), na.rm = TRUE),
    max_abs_package_delta_se = max(abs(one$package_delta_se), na.rm = TRUE),
    min_gmm_rank = min(one$gmm_rank, na.rm = TRUE),
    min_gmm_eigenvalue = min(one$gmm_min_eigenvalue, na.rm = TRUE),
    median_cs_scale = stats::median(one$cs_center_scale, na.rm = TRUE),
    median_seconds_per_replicate = stats::median(one$elapsed_seconds, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))

dir.create(opts$out_dir, recursive = TRUE, showWarnings = FALSE)
replicate_path <- file.path(opts$out_dir, "replicates.csv")
summary_path <- file.path(opts$out_dir, "summary.csv")
diagnostic_path <- file.path(opts$out_dir, "diagnostics.csv")
metadata_path <- file.path(opts$out_dir, "metadata.csv")
write.csv(replicates, replicate_path, row.names = FALSE)
write.csv(summary, summary_path, row.names = FALSE)
write.csv(diagnostics, diagnostic_path, row.names = FALSE)

metadata <- data.frame(
  key = c("experiment", "mode", "reps", "n_grid", "raters", "categories",
          "target", "knowledge_probability", "seed_base", "cores",
          "methods", "elapsed_seconds", "R_version", "misskappa_version",
          "misskappa_library", "generated_at"),
  value = c("47-small-n-kappa-inference", if (opts$smoke) "smoke" else "focused",
            opts$reps, paste(opts$n_grid, collapse = ","), opts$raters,
            opts$categories, opts$target, knowledge_probability, opts$seed_base,
            opts$cores, paste(names(method_labels), collapse = ","), elapsed_all,
            R.version.string, as.character(utils::packageVersion("misskappa")),
            find.package("misskappa"), format(Sys.time(), tz = "UTC", usetz = TRUE)),
  stringsAsFactors = FALSE
)
write.csv(metadata, metadata_path, row.names = FALSE)

cat("Wrote:\n")
for (path in c(replicate_path, summary_path, diagnostic_path, metadata_path)) {
  cat("  ", path, "\n", sep = "")
}
