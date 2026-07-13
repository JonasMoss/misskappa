#!/usr/bin/env Rscript

# Experiment 48: delete-one jackknife studentization and an explicit
# exact-agreement boundary fallback for small-n quadratic Fleiss kappa.

suppressPackageStartupMessages({
  library(misskappa)
})

usage <- function(status = 0L) {
  cat(paste0(
    "Usage: Rscript run_experiment.R [options]\n\n",
    "Options:\n",
    "  --smoke                 Mechanical check: n=20, reps=4, bootstrap B=29.\n",
    "  --reps N                Replications per n. Default: 1000.\n",
    "  --bootstrap-reps N      Subject bootstrap draws per replication. Default: 399.\n",
    "  --n-grid CSV            Sample sizes. Default: 20,30.\n",
    "  --raters N              Number of raters. Default: 5.\n",
    "  --categories N          Ordered categories. Default: 5.\n",
    "  --target K              True quadratic Fleiss kappa. Default: 0.90.\n",
    "  --seed-base N           Deterministic seed base. Default: 474700.\n",
    "  --cores N               Fork workers on non-Windows systems. Default: 1.\n",
    "  --out-dir PATH          Output directory. Default: script-local results/.\n",
    "  --progress              Print one line per sample-size cell.\n",
    "  --help, -h              Show this help and exit.\n\n",
    "Examples:\n",
    "  Rscript run_experiment.R --smoke --progress\n",
    "  Rscript run_experiment.R --reps 1000 --bootstrap-reps 399 --progress\n"
  ))
  quit(save = "no", status = status)
}

script_dir <- local({
  cmd <- commandArgs(FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) == 0L) getwd()
  else dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]),
                             mustWork = FALSE))
})
source(file.path(script_dir, "R", "interval_methods.R"), local = TRUE)

parse_int_csv <- function(x, flag) {
  out <- as.integer(strsplit(x, ",", fixed = TRUE)[[1L]])
  if (any(is.na(out))) {
    stop(flag, " must be a comma-separated integer list.", call. = FALSE)
  }
  unique(out)
}

parse_args <- function(argv) {
  opts <- list(
    smoke = FALSE, reps = 1000L, bootstrap_reps = 399L,
    n_grid = c(20L, 30L), raters = 5L, categories = 5L, target = 0.90,
    seed_base = 474700L, cores = 1L,
    out_dir = file.path(script_dir, "results"), progress = FALSE
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
      "--reps", "--bootstrap-reps", "--n-grid", "--raters",
      "--categories", "--target", "--seed-base", "--cores", "--out-dir"
    )
    if (arg %in% needs_value) {
      if (i == length(argv)) stop(arg, " needs a value.", call. = FALSE)
      value <- argv[[i + 1L]]
      if (arg == "--reps") opts$reps <- as.integer(value)
      else if (arg == "--bootstrap-reps") opts$bootstrap_reps <- as.integer(value)
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
    opts$bootstrap_reps <- 29L
    opts$n_grid <- 20L
    opts$cores <- 1L
  }
  if (is.na(opts$reps) || opts$reps < 1L) stop("--reps must be >= 1.", call. = FALSE)
  if (is.na(opts$bootstrap_reps) || opts$bootstrap_reps < 20L) {
    stop("--bootstrap-reps must be >= 20.", call. = FALSE)
  }
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
min_valid_bootstrap_fraction <- 0.70
method_labels <- c(
  wald_if_z = "IF Wald (normal)",
  wald_if_fisher = "IF Wald (Fisher)",
  jack_t = "Jackknife SE, t",
  jack_fisher_delta_t = "Jackknife SE, delta Fisher-t",
  jack_fisher_direct_t = "Direct Fisher jackknife-t",
  jack_fisher_direct_boundary = "Direct Fisher jackknife-t + boundary",
  bca = "BCa bootstrap",
  bca_boundary = "BCa + boundary",
  bootstrap_t_jack = "Bootstrap-t, jackknife SE",
  bootstrap_t_jack_boundary = "Bootstrap-t, jackknife SE + boundary"
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
  bootstrap_seed <- seed + 480000000L
  started <- proc.time()[["elapsed"]]
  tryCatch({
    x <- simulate_perreault_leigh(n, seed)
    counts <- ratings_to_counts(x, categories = values)
    fit <- kappa_counts(counts, estimator = "fleiss_cuzick",
                        weight = "quadratic", values = values)
    estimate <- unname(stats::coef(fit)[["Fleiss"]])
    se <- sqrt(unname(stats::vcov(fit)["Fleiss", "Fleiss"]))
    focused <- .focused_intervals(
      x, estimate, se, bootstrap_reps = opts$bootstrap_reps,
      min_valid_fraction = min_valid_bootstrap_fraction,
      bootstrap_seed = bootstrap_seed,
      score_range = max(values) - min(values)
    )
    intervals <- focused$intervals
    data.frame(
      n = n, rep = rep, seed = seed, bootstrap_seed = bootstrap_seed,
      method = rownames(intervals), estimate = estimate, if_se = se,
      jackknife_se = focused$jackknife_se,
      lower = intervals[, "lower"], upper = intervals[, "upper"],
      truth = opts$target,
      fallback_used = unname(focused$fallback_used),
      disagreement_events = focused$disagreement_events,
      exact_agreement = focused$exact_agreement,
      leave_one_out_boundary = focused$leave_one_out_boundary,
      valid_bootstrap_t_fraction = focused$valid_bootstrap_t_fraction,
      bootstrap_boundary_fraction = focused$bootstrap_boundary_fraction,
      package_delta_estimate = focused$package_delta_estimate,
      package_delta_se = focused$package_delta_se,
      elapsed_seconds = proc.time()[["elapsed"]] - started,
      status = "ok", stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      n = n, rep = rep, seed = seed, bootstrap_seed = bootstrap_seed,
      method = names(method_labels), estimate = NA_real_, if_se = NA_real_,
      jackknife_se = NA_real_, lower = NA_real_, upper = NA_real_,
      truth = opts$target, fallback_used = FALSE,
      disagreement_events = NA_integer_, exact_agreement = NA,
      leave_one_out_boundary = NA, valid_bootstrap_t_fraction = NA_real_,
      bootstrap_boundary_fraction = NA_real_, package_delta_estimate = NA_real_,
      package_delta_se = NA_real_,
      elapsed_seconds = proc.time()[["elapsed"]] - started,
      status = conditionMessage(e), stringsAsFactors = FALSE
    )
  })
}

run_cell <- function(n) {
  if (opts$progress) {
    message(sprintf("n=%d, reps=%d, bootstrap B=%d, cores=%d",
                    n, opts$reps, opts$bootstrap_reps, opts$cores))
  }
  reps <- seq_len(opts$reps)
  if (opts$cores > 1L && .Platform$OS.type != "windows") {
    pieces <- parallel::mclapply(
      reps, function(rep) one_replicate(n, rep), mc.cores = opts$cores,
      mc.preschedule = TRUE
    )
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

summary_pieces <- split(
  replicates, interaction(replicates$n, replicates$method, drop = TRUE)
)
summary <- do.call(rbind, lapply(summary_pieces, function(d) {
  valid <- d$valid
  ok <- d$status == "ok"
  coverage <- mean(d$covered)
  data.frame(
    n = d$n[[1L]], method = d$method[[1L]],
    method_label = d$method_label[[1L]], reps = nrow(d),
    valid_reps = sum(valid), failures = sum(!valid),
    coverage = coverage,
    coverage_mcse = sqrt(coverage * (1 - coverage) / nrow(d)),
    coverage_valid = if (any(valid)) mean(d$covered[valid]) else NA_real_,
    miss_below = if (any(valid)) mean(d$miss_below[valid]) else NA_real_,
    miss_above = if (any(valid)) mean(d$miss_above[valid]) else NA_real_,
    mean_length = if (any(valid)) mean(d$length[valid]) else NA_real_,
    median_length = if (any(valid)) stats::median(d$length[valid]) else NA_real_,
    fallback_rate = mean(d$fallback_used[ok], na.rm = TRUE),
    bias = mean(d$estimate[ok], na.rm = TRUE) - opts$target,
    empirical_sd = stats::sd(d$estimate[ok], na.rm = TRUE),
    mean_if_se = mean(d$if_se[ok], na.rm = TRUE),
    mean_jackknife_se = mean(d$jackknife_se[ok], na.rm = TRUE),
    if_se_ratio = mean(d$if_se[ok], na.rm = TRUE) /
      stats::sd(d$estimate[ok], na.rm = TRUE),
    jackknife_se_ratio = mean(d$jackknife_se[ok], na.rm = TRUE) /
      stats::sd(d$estimate[ok], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
summary <- summary[order(summary$n, match(summary$method, names(method_labels))), ]
rownames(summary) <- NULL

diagnostics <- do.call(rbind, lapply(split(replicates, replicates$n), function(d) {
  one <- d[d$method == "wald_if_z" & d$status == "ok", ]
  data.frame(
    n = d$n[[1L]],
    exact_agreement_reps = sum(one$exact_agreement),
    exact_agreement_rate = mean(one$exact_agreement),
    leave_one_out_boundary_reps = sum(one$leave_one_out_boundary),
    leave_one_out_boundary_rate = mean(one$leave_one_out_boundary),
    one_or_fewer_disagreement_events = sum(one$disagreement_events <= 1L),
    mean_valid_bootstrap_t_fraction = mean(one$valid_bootstrap_t_fraction),
    min_valid_bootstrap_t_fraction = min(one$valid_bootstrap_t_fraction),
    mean_bootstrap_boundary_fraction = mean(one$bootstrap_boundary_fraction),
    max_abs_package_delta_estimate = max(abs(one$package_delta_estimate)),
    max_abs_package_delta_se = max(abs(one$package_delta_se)),
    median_seconds_per_replicate = stats::median(one$elapsed_seconds),
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
  key = c(
    "experiment", "mode", "reps", "bootstrap_reps", "n_grid", "raters",
    "categories", "target", "knowledge_probability", "seed_base", "cores",
    "min_valid_bootstrap_fraction", "methods", "elapsed_seconds", "R_version",
    "misskappa_version", "misskappa_library", "generated_at"
  ),
  value = c(
    "48-small-n-jackknife-boundary", if (opts$smoke) "smoke" else "focused",
    opts$reps, opts$bootstrap_reps, paste(opts$n_grid, collapse = ","),
    opts$raters, opts$categories, opts$target, knowledge_probability,
    opts$seed_base, opts$cores, min_valid_bootstrap_fraction,
    paste(names(method_labels), collapse = ","), elapsed_all,
    R.version.string, as.character(utils::packageVersion("misskappa")),
    find.package("misskappa"), format(Sys.time(), tz = "UTC", usetz = TRUE)
  ),
  stringsAsFactors = FALSE
)
write.csv(metadata, metadata_path, row.names = FALSE)

cat("Wrote:\n")
for (path in c(replicate_path, summary_path, diagnostic_path, metadata_path)) {
  cat("  ", path, "\n", sep = "")
}
