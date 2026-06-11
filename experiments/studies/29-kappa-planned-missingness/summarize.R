# Shared aggregation helpers for Experiment 29, sourced by both
# run_experiment.R (single-process end-of-run summary) and combine.R (merge of
# per-cell checkpoints across replicate shards). Keeping one definition here
# guarantees the sharded merge and the single-process run summarise identically.

split_keys <- function(data, keys) {
  interaction(data[, keys], drop = TRUE, lex.order = TRUE)
}

summarize_replicates <- function(df) {
  keys <- c("dgp", "dgp_label", "dgp_family", "C", "R", "mechanism",
            "mechanism_family", "n", "method", "estimator", "weight_label",
            "coefficient")
  interval_names <- c("wald_z", "wald_t", "fisher_t", "asin_t")
  pieces <- lapply(split(df, split_keys(df, keys)), function(g) {
    ok <- is.finite(g$estimate)
    se_ok <- ok & is.finite(g$se) & g$se > 0
    err <- g$estimate - g$truth
    base <- data.frame(
      dgp = g$dgp[[1L]],
      dgp_label = g$dgp_label[[1L]],
      dgp_family = g$dgp_family[[1L]],
      C = g$C[[1L]],
      R = g$R[[1L]],
      mechanism = g$mechanism[[1L]],
      mechanism_family = g$mechanism_family[[1L]],
      n = g$n[[1L]],
      method = g$method[[1L]],
      estimator = g$estimator[[1L]],
      weight_label = g$weight_label[[1L]],
      coefficient = g$coefficient[[1L]],
      reps = length(unique(g$rep)),
      n_valid = sum(ok),
      failures = sum(!ok),
      truth = g$truth[[1L]],
      mean_estimate = if (any(ok)) mean(g$estimate[ok]) else NA_real_,
      bias = if (any(ok)) mean(err[ok]) else NA_real_,
      sd_estimate = if (sum(ok) > 1L) stats::sd(g$estimate[ok]) else NA_real_,
      mc_se_bias = if (sum(ok) > 1L) stats::sd(g$estimate[ok]) / sqrt(sum(ok)) else NA_real_,
      mse = if (any(ok)) mean(err[ok]^2) else NA_real_,
      rmse = if (any(ok)) sqrt(mean(err[ok]^2)) else NA_real_,
      mean_se = if (any(se_ok)) mean(g$se[se_ok]) else NA_real_,
      se_over_sd = if (sum(se_ok) > 1L && stats::sd(g$estimate[se_ok]) > 1e-12) {
        mean(g$se[se_ok]) / stats::sd(g$estimate[se_ok])
      } else NA_real_,
      mean_elapsed_ms = mean(g$elapsed_ms, na.rm = TRUE),
      median_elapsed_ms = stats::median(g$elapsed_ms, na.rm = TRUE),
      mean_observed_fraction = mean(g$observed_fraction, na.rm = TRUE),
      mean_subjects_used = mean(g$subjects_used, na.rm = TRUE),
      mean_empty_rows = mean(g$empty_rows, na.rm = TRUE),
      mean_complete_rows = mean(g$complete_rows, na.rm = TRUE),
      mean_n_eff = mean(g$n_eff, na.rm = TRUE),
      min_pair_count_min = min(g$min_pair_count, na.rm = TRUE),
      mean_observed_patterns = mean(g$observed_patterns, na.rm = TRUE),
      stringsAsFactors = FALSE
    )

    do.call(rbind, lapply(interval_names, function(interval) {
      lo <- g[[paste0(interval, "_lower")]]
      hi <- g[[paste0(interval, "_upper")]]
      ci_ok <- se_ok & is.finite(lo) & is.finite(hi)
      covered <- ci_ok & lo <= g$truth & g$truth <= hi
      below <- ci_ok & hi < g$truth
      above <- ci_ok & lo > g$truth
      cbind(
        base,
        data.frame(
          interval = interval,
          interval_n = sum(ci_ok),
          coverage95 = if (any(ci_ok)) mean(covered[ci_ok]) else NA_real_,
          miss_below = if (any(ci_ok)) mean(below[ci_ok]) else NA_real_,
          miss_above = if (any(ci_ok)) mean(above[ci_ok]) else NA_real_,
          mean_ci_length = if (any(ci_ok)) mean(hi[ci_ok] - lo[ci_ok]) else NA_real_,
          stringsAsFactors = FALSE
        )
      )
    }))
  })
  ans <- do.call(rbind, pieces)
  ans <- ans[order(ans$dgp, ans$mechanism, ans$n, ans$weight_label,
                   ans$estimator, ans$coefficient, ans$interval), ]
  rownames(ans) <- NULL
  ans
}

write_csv_atomic <- function(x, path) {
  tmp <- paste0(path, ".tmp")
  write.csv(x, tmp, row.names = FALSE)
  if (file.exists(path)) unlink(path)
  ok <- file.rename(tmp, path)
  if (!ok) stop("Failed to move temporary file into place: ", path, call. = FALSE)
}
