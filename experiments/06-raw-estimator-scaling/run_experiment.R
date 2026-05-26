#!/usr/bin/env Rscript
#
# Experiment 06: raw estimator scaling.
#
# Measures elapsed time for the raw estimators across sample sizes and compares
# them to counts-format estimators on the same simulated ratings. This is an
# engineering experiment for the raw O(n^2) chance-kernel cleanup.

usage <- function(status = 0L) {
  cat(paste0(
    "Usage: Rscript run_experiment.R [options]\n\n",
    "Options:\n",
    "  --help              Show this help and exit.\n",
    "  --smoke             Cheap mechanical check.\n",
    "  --overnight         Larger grid intended for unattended runs.\n",
    "  --n-grid CSV        Sample sizes. Default: 250,500,1000,2000.\n",
    "  --reps N            Replicates per sample size. Default: 3.\n",
    "  --seed-base N       Base seed for deterministic runs.\n",
    "  --out-dir PATH      Output directory. Default: script-local results/.\n",
    "  --progress          Print one line per timed fit.\n\n",
    "Examples:\n",
    "  Rscript run_experiment.R --smoke\n",
    "  Rscript run_experiment.R --overnight --progress\n"
  ))
  quit(save = "no", status = status)
}

parse_int_csv <- function(x) {
  out <- as.integer(strsplit(x, ",", fixed = TRUE)[[1L]])
  if (any(is.na(out)) || any(out < 2L)) stop("--n-grid must contain integers >= 2.", call. = FALSE)
  unique(out)
}

parse_args <- function(argv) {
  opts <- list(
    smoke = FALSE,
    overnight = FALSE,
    n_grid = c(250L, 500L, 1000L, 2000L),
    reps = 3L,
    seed_base = 606060L,
    out_dir = NULL,
    progress = FALSE
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
    if (arg == "--overnight") {
      opts$overnight <- TRUE
      i <- i + 1L
      next
    }
    if (arg == "--progress") {
      opts$progress <- TRUE
      i <- i + 1L
      next
    }
    needs_value <- c("--n-grid", "--reps", "--seed-base", "--out-dir")
    if (arg %in% needs_value) {
      if (i == length(argv)) stop(arg, " needs a value.", call. = FALSE)
      val <- argv[[i + 1L]]
      if (arg == "--n-grid") opts$n_grid <- parse_int_csv(val)
      if (arg == "--reps") opts$reps <- as.integer(val)
      if (arg == "--seed-base") opts$seed_base <- as.integer(val)
      if (arg == "--out-dir") opts$out_dir <- val
      i <- i + 2L
      next
    }
    stop("Unknown argument: ", arg, call. = FALSE)
  }
  if (opts$smoke) {
    opts$n_grid <- c(80L, 160L)
    opts$reps <- 1L
  }
  if (opts$overnight) {
    opts$n_grid <- c(250L, 500L, 1000L, 2000L, 4000L, 8000L)
    opts$reps <- max(opts$reps, 3L)
  }
  if (is.na(opts$reps) || opts$reps < 1L) stop("--reps must be >= 1.", call. = FALSE)
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
p <- colSums(as.matrix(dat.fleiss1971)) / sum(as.matrix(dat.fleiss1971))
pi_rater <- c(0.95, 0.90, 0.80, 0.65, 0.50, 0.35)

raw_methods <- c("available", "ipw", "gwet", "fiml")
counts_methods <- c("available", "fiml")

simulate_raw <- function(n) {
  x_star <- matrix(sample.int(C, n * R, replace = TRUE, prob = p), nrow = n, ncol = R)
  x <- x_star
  for (j in seq_len(R)) {
    observed <- stats::runif(n) < pi_rater[j]
    x[!observed, j] <- NA_integer_
  }
  x
}

to_counts <- function(x) {
  out <- matrix(0L, nrow = nrow(x), ncol = C)
  for (i in seq_len(nrow(x))) {
    vals <- x[i, ]
    vals <- vals[!is.na(vals)]
    if (length(vals) > 0L) {
      tab <- tabulate(vals, nbins = C)
      out[i, ] <- as.integer(tab)
    }
  }
  out
}

time_call <- function(expr) {
  gc(FALSE)
  start <- proc.time()
  result <- tryCatch(
    {
      value <- force(expr)
      list(ok = TRUE, value = value, error = "")
    },
    error = function(e) list(ok = FALSE, value = NULL, error = conditionMessage(e))
  )
  elapsed <- (proc.time() - start)[["elapsed"]]
  list(elapsed = as.numeric(elapsed), ok = result$ok, error = result$error)
}

fit_raw <- function(x, method) {
  args <- list(x = x, method = method, weight = "identity")
  if (method == "fiml") {
    args$em_options <- list(max_iter = 50000L, tol = 1e-7)
  }
  misskappa::kappa(x = args$x, method = args$method, weight = args$weight,
                   em_options = args$em_options %||% list())
}

`%||%` <- function(x, y) if (is.null(x)) y else x

fit_counts <- function(counts, method) {
  misskappa::kappa_counts(
    counts,
    method = method,
    weight = "identity",
    r_total = R,
    em_options = list(max_iter = 50000L, tol = 1e-7)
  )
}

rows <- list()
pos <- 1L
total_fits <- length(opts$n_grid) * opts$reps * (length(raw_methods) + length(counts_methods))
fit_id <- 0L

for (n in opts$n_grid) {
  for (rep in seq_len(opts$reps)) {
    seed <- opts$seed_base + n * 1000L + rep
    set.seed(seed)
    x <- simulate_raw(n)
    counts <- to_counts(x)

    for (method in raw_methods) {
      fit_id <- fit_id + 1L
      if (opts$progress) {
        message(sprintf("[%d/%d] n=%d rep=%d raw:%s", fit_id, total_fits, n, rep, method))
      }
      t <- time_call(fit_raw(x, method))
      rows[[pos]] <- data.frame(
        input = "raw", method = method, n = n, rep = rep, seed = seed,
        elapsed_sec = t$elapsed, ok = t$ok, error = t$error,
        stringsAsFactors = FALSE
      )
      pos <- pos + 1L
    }

    for (method in counts_methods) {
      fit_id <- fit_id + 1L
      if (opts$progress) {
        message(sprintf("[%d/%d] n=%d rep=%d counts:%s", fit_id, total_fits, n, rep, method))
      }
      t <- time_call(fit_counts(counts, method))
      rows[[pos]] <- data.frame(
        input = "counts", method = method, n = n, rep = rep, seed = seed,
        elapsed_sec = t$elapsed, ok = t$ok, error = t$error,
        stringsAsFactors = FALSE
      )
      pos <- pos + 1L
    }
  }
}

timings <- do.call(rbind, rows)

group_key <- interaction(timings$input, timings$method, timings$n, drop = TRUE, lex.order = TRUE)
summary <- do.call(rbind, lapply(split(timings, group_key), function(g) {
  ok <- g$ok
  elapsed <- g$elapsed_sec[ok]
  data.frame(
    input = g$input[[1L]],
    method = g$method[[1L]],
    n = g$n[[1L]],
    reps = nrow(g),
    n_ok = sum(ok),
    n_error = sum(!ok),
    median_elapsed_sec = if (length(elapsed) > 0L) stats::median(elapsed) else NA_real_,
    mean_elapsed_sec = if (length(elapsed) > 0L) mean(elapsed) else NA_real_,
    min_elapsed_sec = if (length(elapsed) > 0L) min(elapsed) else NA_real_,
    max_elapsed_sec = if (length(elapsed) > 0L) max(elapsed) else NA_real_,
    stringsAsFactors = FALSE
  )
}))
summary <- summary[order(summary$input, summary$method, summary$n), ]
rownames(summary) <- NULL

raw_available <- subset(summary, input == "raw" & method == "available")
if (nrow(raw_available) >= 2L) {
  fit <- stats::lm(log(median_elapsed_sec) ~ log(n), data = raw_available)
  scaling_exponent <- unname(stats::coef(fit)[["log(n)"]])
} else {
  scaling_exponent <- NA_real_
}

metadata <- data.frame(
  key = c(
    "timestamp", "smoke", "overnight", "n_grid", "reps", "seed_base",
    "C", "R", "raw_methods", "counts_methods", "misskappa_version",
    "r_version", "raw_available_loglog_exponent"
  ),
  value = c(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
    as.character(opts$smoke),
    as.character(opts$overnight),
    paste(opts$n_grid, collapse = ";"),
    as.character(opts$reps),
    as.character(opts$seed_base),
    as.character(C),
    as.character(R),
    paste(raw_methods, collapse = ";"),
    paste(counts_methods, collapse = ";"),
    as.character(utils::packageVersion("misskappa")),
    paste(R.version$major, R.version$minor, sep = "."),
    as.character(scaling_exponent)
  ),
  stringsAsFactors = FALSE
)

dir.create(opts$out_dir, recursive = TRUE, showWarnings = FALSE)
metadata_path <- file.path(opts$out_dir, "metadata.csv")
timings_path <- file.path(opts$out_dir, "timings.csv")
summary_path <- file.path(opts$out_dir, "summary.csv")

write.csv(metadata, metadata_path, row.names = FALSE)
write.csv(timings, timings_path, row.names = FALSE)
write.csv(summary, summary_path, row.names = FALSE)

cat("Wrote ", metadata_path, "\n", sep = "")
cat("Wrote ", timings_path, "\n", sep = "")
cat("Wrote ", summary_path, "\n", sep = "")
