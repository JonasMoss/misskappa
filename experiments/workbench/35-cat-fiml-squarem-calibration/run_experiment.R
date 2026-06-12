#!/usr/bin/env Rscript
#
# Experiment 35: compare ordinary strict-ML cat-FIML EM against the SQUAREM
# accelerated variant on study-29-shaped cells. The C++ helper is compiled
# twice against local copies of the instrumented estimator variants.

usage <- function(status = 0L) {
  cat(paste0(
    "Usage: Rscript run_experiment.R [mode/options]\n\n",
    "Modes:\n",
    "  --smoke          Cheap mechanical check: reps=2, n=40.\n",
    "  --calibration    Broader local run: reps=300, n=40,100,300.\n\n",
    "Options:\n",
    "  --reps N            Override replications per cell.\n",
    "  --n-grid CSV        Override sample sizes.\n",
    "  --truth-n N         Complete-data Monte Carlo truth size.\n",
    "  --dgps CSV          DGP ids: balanced6,highagree6,sparse6,sparse5r7.\n",
    "  --mechanisms CSV    Mechanisms: complete,mcar30,mar_anchor30,designed_random2,designed_bib3.\n",
    "  --weights CSV       Weights: identity,quadratic.\n",
    "  --seed-base N       Base seed. Default: 353500.\n",
    "  --out-dir PATH      Output directory. Default: script-local results/.\n",
    "  --progress          Print one line per design cell.\n",
    "  --help, -h          Show this help and exit.\n\n",
    "Examples:\n",
    "  Rscript run_experiment.R --smoke --progress\n",
    "  Rscript run_experiment.R --calibration --reps 500 --progress\n"
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

mode_defaults <- function(mode) {
  if (mode == "smoke") {
    return(list(
      mode = mode,
      reps = 2L,
      n_grid = 40L,
      truth_n = 30000L,
      dgps = c("balanced6", "sparse6"),
      mechanisms = c("mcar30", "designed_random2"),
      weights = c("identity", "quadratic")
    ))
  }
  if (mode == "calibration") {
    return(list(
      mode = mode,
      reps = 300L,
      n_grid = c(40L, 100L, 300L),
      truth_n = 200000L,
      dgps = c("balanced6", "highagree6", "sparse6"),
      mechanisms = c("mcar30", "mar_anchor30", "designed_random2", "designed_bib3"),
      weights = c("identity", "quadratic")
    ))
  }
  stop("Unknown mode: ", mode, call. = FALSE)
}

parse_args <- function(argv) {
  mode <- "smoke"
  if ("--smoke" %in% argv) mode <- "smoke"
  if ("--calibration" %in% argv) mode <- "calibration"
  opts <- mode_defaults(mode)
  opts$seed_base <- 353500L
  opts$out_dir <- file.path(script_dir, "results")
  opts$progress <- FALSE

  i <- 1L
  while (i <= length(argv)) {
    arg <- argv[[i]]
    if (arg %in% c("--help", "-h")) usage(0L)
    if (arg %in% c("--smoke", "--calibration")) {
      i <- i + 1L
      next
    }
    if (arg == "--progress") {
      opts$progress <- TRUE
      i <- i + 1L
      next
    }
    needs_value <- c("--reps", "--n-grid", "--truth-n", "--dgps",
                     "--mechanisms", "--weights", "--seed-base", "--out-dir")
    if (arg %in% needs_value) {
      if (i == length(argv)) stop(arg, " needs a value.", call. = FALSE)
      val <- argv[[i + 1L]]
      if (arg == "--reps") opts$reps <- as.integer(val)
      else if (arg == "--n-grid") opts$n_grid <- parse_csv_int(val, arg)
      else if (arg == "--truth-n") opts$truth_n <- as.integer(val)
      else if (arg == "--dgps") opts$dgps <- parse_csv_chr(val)
      else if (arg == "--mechanisms") opts$mechanisms <- parse_csv_chr(val)
      else if (arg == "--weights") opts$weights <- parse_csv_chr(val)
      else if (arg == "--seed-base") opts$seed_base <- as.integer(val)
      else if (arg == "--out-dir") opts$out_dir <- val
      i <- i + 2L
      next
    }
    stop("Unknown argument: ", arg, call. = FALSE)
  }

  if (is.na(opts$reps) || opts$reps < 1L) stop("--reps must be >= 1.", call. = FALSE)
  if (any(is.na(opts$n_grid)) || any(opts$n_grid < 5L)) {
    stop("--n-grid entries must be integers >= 5.", call. = FALSE)
  }
  if (is.na(opts$truth_n) || opts$truth_n < 1000L) {
    stop("--truth-n must be an integer >= 1000.", call. = FALSE)
  }
  if (is.na(opts$seed_base)) stop("--seed-base must be an integer.", call. = FALSE)
  opts
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))

valid_dgps <- c("balanced6", "highagree6", "sparse6", "sparse5r7")
valid_mechanisms <- c("complete", "mcar30", "mar_anchor30",
                      "designed_random2", "designed_bib3")
valid_weights <- c("identity", "nominal", "quadratic")
if (!all(opts$dgps %in% valid_dgps)) {
  stop("--dgps must contain only: ", paste(valid_dgps, collapse = ","), call. = FALSE)
}
if (!all(opts$mechanisms %in% valid_mechanisms)) {
  stop("--mechanisms must contain only: ", paste(valid_mechanisms, collapse = ","), call. = FALSE)
}
if (!all(opts$weights %in% valid_weights)) {
  stop("--weights must contain only: ", paste(valid_weights, collapse = ","), call. = FALSE)
}
opts$weights <- ifelse(opts$weights == "nominal", "identity", opts$weights)

write_csv_atomic <- function(x, path) {
  tmp <- paste0(path, ".tmp")
  write.csv(x, tmp, row.names = FALSE)
  if (file.exists(path)) unlink(path)
  ok <- file.rename(tmp, path)
  if (!ok) stop("Failed to move temporary file into place: ", path, call. = FALSE)
}

find_eigen <- function() {
  candidates <- unique(c(
    Sys.getenv("EIGEN", unset = NA_character_),
    file.path(script_dir, "../../../build/dev/_deps/eigen3-src"),
    file.path(script_dir, "../../../build/opt/_deps/eigen3-src"),
    "/usr/include/eigen3",
    "/usr/local/include/eigen3",
    "/home/jonas/Files/research/magmaan/build/fast/_deps/eigen3-src"
  ))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  hits <- candidates[file.exists(file.path(candidates, "Eigen", "Core"))]
  if (!length(hits)) {
    stop("Could not find Eigen. Set EIGEN=/path/to/eigen3-src.", call. = FALSE)
  }
  normalizePath(hits[[1L]], mustWork = TRUE)
}

compile_variant <- function(variant, source, eigen_dir) {
  bin_dir <- file.path(script_dir, "bin")
  dir.create(bin_dir, recursive = TRUE, showWarnings = FALSE)
  exe <- file.path(bin_dir, paste0("calibration_", variant))
  cxx <- Sys.getenv("CXX", unset = "g++")
  args <- c(
    "-std=c++17", "-O3", "-g", "-DNDEBUG", "-DEIGEN_NO_EXCEPTIONS",
    paste0("-I", file.path(script_dir, "../../../include")),
    paste0("-I", file.path(script_dir, "../../../src")),
    paste0("-I", script_dir),
    paste0("-I", eigen_dir),
    file.path(script_dir, "calibration.cpp"),
    file.path(script_dir, source),
    file.path(script_dir, "../../../src/loss.cpp"),
    file.path(script_dir, "../../../src/estimate_raw.cpp"),
    file.path(script_dir, "../../../src/detail_inverse_weights.cpp"),
    "-o", exe
  )
  message("Compiling ", variant, " helper.")
  status <- system2(cxx, args = args)
  if (!identical(status, 0L)) stop("Compilation failed for variant: ", variant, call. = FALSE)
  normalizePath(exe, mustWork = TRUE)
}

run_variant <- function(variant, exe) {
  out <- file.path(opts$out_dir, paste0("replicates-", variant, ".csv"))
  args <- c(
    "--variant", variant,
    "--out", out,
    "--reps", as.character(opts$reps),
    "--truth-n", as.character(opts$truth_n),
    "--seed-base", as.character(opts$seed_base),
    "--n-grid", paste(opts$n_grid, collapse = ","),
    "--dgps", paste(opts$dgps, collapse = ","),
    "--mechanisms", paste(opts$mechanisms, collapse = ","),
    "--weights", paste(opts$weights, collapse = ",")
  )
  if (opts$progress) args <- c(args, "--progress")
  message("Running ", variant, " helper.")
  status <- system2(exe, args = args)
  if (!identical(status, 0L)) stop("Run failed for variant: ", variant, call. = FALSE)
  out
}

clip_unit <- function(x) pmin(pmax(x, -1 + 1e-10), 1 - 1e-10)

add_intervals <- function(df) {
  z <- stats::qnorm(0.975)
  tcrit <- stats::qt(0.975, df = pmax(1, df$n_eff - 1L))
  ok <- is.finite(df$estimate) & is.finite(df$se) & df$se >= 0 &
    is.finite(df$n_eff) & df$n_eff >= 2

  df$wald_z_lower <- df$estimate - z * df$se
  df$wald_z_upper <- df$estimate + z * df$se
  df$wald_t_lower <- df$estimate - tcrit * df$se
  df$wald_t_upper <- df$estimate + tcrit * df$se

  est_c <- clip_unit(df$estimate)
  se_f <- df$se / (1 - est_c^2)
  g_f <- atanh(est_c)
  df$fisher_t_lower <- tanh(g_f - tcrit * se_f)
  df$fisher_t_upper <- tanh(g_f + tcrit * se_f)

  se_a <- df$se / sqrt(1 - est_c^2)
  g_a <- asin(est_c)
  df$asin_t_lower <- sin(g_a - tcrit * se_a)
  df$asin_t_upper <- sin(g_a + tcrit * se_a)

  interval_cols <- grep("_(lower|upper)$", names(df), value = TRUE)
  df[!ok, interval_cols] <- NA_real_
  df
}

split_keys <- function(data, keys) {
  interaction(data[, keys], drop = TRUE, lex.order = TRUE)
}

summarize_replicates <- function(df) {
  keys <- c("variant", "method", "estimator", "dgp", "dgp_label", "dgp_family",
            "C", "R", "mechanism", "mechanism_family", "n", "weight_label",
            "coefficient")
  interval_names <- c("wald_z", "wald_t", "fisher_t", "asin_t")
  pieces <- lapply(split(df, split_keys(df, keys)), function(g) {
    ok <- is.finite(g$estimate)
    se_ok <- ok & is.finite(g$se) & g$se > 0
    err <- g$estimate - g$truth
    base <- data.frame(
      variant = g$variant[[1L]],
      method = g$method[[1L]],
      estimator = g$estimator[[1L]],
      dgp = g$dgp[[1L]],
      dgp_label = g$dgp_label[[1L]],
      dgp_family = g$dgp_family[[1L]],
      C = g$C[[1L]],
      R = g$R[[1L]],
      mechanism = g$mechanism[[1L]],
      mechanism_family = g$mechanism_family[[1L]],
      n = g$n[[1L]],
      weight_label = g$weight_label[[1L]],
      coefficient = g$coefficient[[1L]],
      reps = length(unique(g$rep)),
      n_valid = sum(ok),
      failures = sum(!ok),
      truth = g$truth[[1L]],
      mean_estimate = if (any(ok)) mean(g$estimate[ok]) else NA_real_,
      bias = if (any(ok)) mean(err[ok]) else NA_real_,
      sd_estimate = if (sum(ok) > 1L) stats::sd(g$estimate[ok]) else NA_real_,
      rmse = if (any(ok)) sqrt(mean(err[ok]^2)) else NA_real_,
      mean_se = if (any(se_ok)) mean(g$se[se_ok]) else NA_real_,
      se_over_sd = if (sum(se_ok) > 1L && stats::sd(g$estimate[se_ok]) > 1e-12) {
        mean(g$se[se_ok]) / stats::sd(g$estimate[se_ok])
      } else NA_real_,
      mean_null_frac = mean(g$null_frac, na.rm = TRUE),
      mean_elapsed_ms = mean(g$elapsed_ms, na.rm = TRUE),
      median_elapsed_ms = stats::median(g$elapsed_ms, na.rm = TRUE),
      mean_em_iters = mean(g$em_iters, na.rm = TRUE),
      median_em_iters = stats::median(g$em_iters, na.rm = TRUE),
      mean_support = mean(g$support, na.rm = TRUE),
      mean_active = mean(g$active, na.rm = TRUE),
      mean_observed_fraction = mean(g$observed_fraction, na.rm = TRUE),
      mean_min_pair_count = mean(g$min_pair_count, na.rm = TRUE),
      stringsAsFactors = FALSE
    )

    do.call(rbind, lapply(interval_names, function(interval) {
      lo <- g[[paste0(interval, "_lower")]]
      hi <- g[[paste0(interval, "_upper")]]
      ci_ok <- se_ok & is.finite(lo) & is.finite(hi)
      covered <- ci_ok & lo <= g$truth & g$truth <= hi
      cbind(
        base,
        data.frame(
          interval = interval,
          interval_n = sum(ci_ok),
          coverage95 = if (any(ci_ok)) mean(covered[ci_ok]) else NA_real_,
          mean_ci_length = if (any(ci_ok)) mean(hi[ci_ok] - lo[ci_ok]) else NA_real_,
          stringsAsFactors = FALSE
        )
      )
    }))
  })
  ans <- do.call(rbind, pieces)
  ans <- ans[order(ans$dgp, ans$mechanism, ans$n, ans$weight_label,
                   ans$coefficient, ans$variant, ans$interval), ]
  rownames(ans) <- NULL
  ans
}

summarize_comparison <- function(df) {
  em <- df[df$variant == "em" & is.finite(df$estimate), ]
  sq <- df[df$variant == "squarem" & is.finite(df$estimate), ]
  by <- c("dgp", "mechanism", "n", "rep", "weight_label", "coefficient")
  m <- merge(em, sq, by = by, suffixes = c("_em", "_squarem"))
  if (!nrow(m)) return(data.frame())
  m$estimate_diff <- m$estimate_squarem - m$estimate_em
  m$abs_estimate_diff <- abs(m$estimate_diff)
  m$abs_estimate_diff_over_em_se <- m$abs_estimate_diff / m$se_em
  m$se_ratio <- m$se_squarem / m$se_em
  m$elapsed_speedup <- m$elapsed_ms_em / m$elapsed_ms_squarem
  m$iter_ratio <- m$em_iters_squarem / m$em_iters_em
  m$support_delta <- m$support_squarem - m$support_em
  m$wald_t_cover_em <- m$wald_t_lower_em <= m$truth_em & m$truth_em <= m$wald_t_upper_em
  m$wald_t_cover_squarem <- m$wald_t_lower_squarem <= m$truth_em & m$truth_em <= m$wald_t_upper_squarem
  m$fisher_t_cover_em <- m$fisher_t_lower_em <= m$truth_em & m$truth_em <= m$fisher_t_upper_em
  m$fisher_t_cover_squarem <- m$fisher_t_lower_squarem <= m$truth_em & m$truth_em <= m$fisher_t_upper_squarem

  pieces <- lapply(split(m, split_keys(m, c("dgp", "mechanism", "n",
                                            "weight_label", "coefficient"))), function(g) {
    ok_se <- is.finite(g$se_ratio) & g$se_ratio > 0
    data.frame(
      dgp = g$dgp[[1L]],
      dgp_label = g$dgp_label_em[[1L]],
      dgp_family = g$dgp_family_em[[1L]],
      C = g$C_em[[1L]],
      R = g$R_em[[1L]],
      mechanism = g$mechanism[[1L]],
      mechanism_family = g$mechanism_family_em[[1L]],
      n = g$n[[1L]],
      weight_label = g$weight_label[[1L]],
      coefficient = g$coefficient[[1L]],
      paired_reps = nrow(g),
      max_abs_estimate_diff = max(g$abs_estimate_diff, na.rm = TRUE),
      median_abs_estimate_diff = stats::median(g$abs_estimate_diff, na.rm = TRUE),
      p95_abs_estimate_diff = unname(stats::quantile(g$abs_estimate_diff, 0.95, na.rm = TRUE)),
      median_abs_estimate_diff_over_em_se = stats::median(g$abs_estimate_diff_over_em_se, na.rm = TRUE),
      median_se_ratio = stats::median(g$se_ratio[ok_se], na.rm = TRUE),
      p10_se_ratio = unname(stats::quantile(g$se_ratio[ok_se], 0.10, na.rm = TRUE)),
      p90_se_ratio = unname(stats::quantile(g$se_ratio[ok_se], 0.90, na.rm = TRUE)),
      mean_support_delta = mean(g$support_delta, na.rm = TRUE),
      mean_speedup = mean(g$elapsed_speedup, na.rm = TRUE),
      median_speedup = stats::median(g$elapsed_speedup, na.rm = TRUE),
      mean_iter_ratio = mean(g$iter_ratio, na.rm = TRUE),
      coverage_wald_t_em = mean(g$wald_t_cover_em, na.rm = TRUE),
      coverage_wald_t_squarem = mean(g$wald_t_cover_squarem, na.rm = TRUE),
      coverage_wald_t_delta = mean(g$wald_t_cover_squarem, na.rm = TRUE) -
        mean(g$wald_t_cover_em, na.rm = TRUE),
      coverage_fisher_t_em = mean(g$fisher_t_cover_em, na.rm = TRUE),
      coverage_fisher_t_squarem = mean(g$fisher_t_cover_squarem, na.rm = TRUE),
      coverage_fisher_t_delta = mean(g$fisher_t_cover_squarem, na.rm = TRUE) -
        mean(g$fisher_t_cover_em, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  ans <- do.call(rbind, pieces)
  ans <- ans[order(ans$dgp, ans$mechanism, ans$n, ans$weight_label, ans$coefficient), ]
  rownames(ans) <- NULL
  ans
}

t0 <- Sys.time()
dir.create(opts$out_dir, recursive = TRUE, showWarnings = FALSE)
eigen_dir <- find_eigen()
executables <- c(
  em = compile_variant("em", "estimate_fiml_em.cpp", eigen_dir),
  squarem = compile_variant("squarem", "estimate_fiml_squarem.cpp", eigen_dir)
)
raw_paths <- mapply(run_variant, names(executables), executables, SIMPLIFY = TRUE)

replicates <- do.call(rbind, lapply(raw_paths, function(path) {
  read.csv(path, stringsAsFactors = FALSE, na.strings = "NA")
}))
replicates <- add_intervals(replicates)
summary <- summarize_replicates(replicates)
comparison <- summarize_comparison(replicates)
truth <- unique(replicates[, c("dgp", "dgp_label", "dgp_family", "C", "R",
                               "weight_label", "coefficient", "truth")])
truth <- truth[order(truth$dgp, truth$weight_label, truth$coefficient), ]
cell_plan <- unique(replicates[, c("dgp", "mechanism", "n")])
cell_plan <- cell_plan[order(cell_plan$dgp, cell_plan$mechanism, cell_plan$n), ]

write_csv_atomic(replicates, file.path(opts$out_dir, "replicates.csv"))
write_csv_atomic(summary, file.path(opts$out_dir, "summary.csv"))
write_csv_atomic(comparison, file.path(opts$out_dir, "comparison.csv"))
write_csv_atomic(truth, file.path(opts$out_dir, "truth.csv"))
write_csv_atomic(cell_plan, file.path(opts$out_dir, "cell_plan.csv"))

pkg_version <- if (requireNamespace("misskappa", quietly = TRUE)) {
  as.character(utils::packageVersion("misskappa"))
} else {
  NA_character_
}
metadata <- data.frame(
  key = c("generated_at", "script", "command", "mode", "reps", "n_grid",
          "truth_n", "dgps", "mechanisms", "weights", "seed_base",
          "eigen_dir", "elapsed_seconds", "misskappa_version", "r_version"),
  value = c(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
    normalizePath(file.path(script_dir, "run_experiment.R"), mustWork = FALSE),
    paste(commandArgs(FALSE), collapse = " "),
    opts$mode,
    as.character(opts$reps),
    paste(opts$n_grid, collapse = ","),
    as.character(opts$truth_n),
    paste(opts$dgps, collapse = ","),
    paste(opts$mechanisms, collapse = ","),
    paste(opts$weights, collapse = ","),
    as.character(opts$seed_base),
    eigen_dir,
    as.character(round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 3)),
    pkg_version,
    R.version.string
  ),
  stringsAsFactors = FALSE
)
write_csv_atomic(metadata, file.path(opts$out_dir, "metadata.csv"))

cat("Wrote:\n")
for (name in c("replicates.csv", "summary.csv", "comparison.csv", "truth.csv",
               "cell_plan.csv", "metadata.csv")) {
  cat("  ", normalizePath(file.path(opts$out_dir, name), mustWork = FALSE), "\n", sep = "")
}

headline <- comparison[comparison$coefficient == "Conger" &
                         comparison$weight_label == "quadratic",
                       c("dgp", "mechanism", "n", "paired_reps",
                         "median_abs_estimate_diff", "median_se_ratio",
                         "coverage_wald_t_em", "coverage_wald_t_squarem",
                         "median_speedup", "mean_iter_ratio")]
headline <- headline[order(headline$dgp, headline$mechanism, headline$n), ]
if (nrow(headline)) {
  cat("\nConger quadratic paired comparison:\n")
  print(utils::head(headline, 24L), row.names = FALSE, digits = 3)
}
