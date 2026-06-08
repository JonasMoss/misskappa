# Aggregate per-replicate records into one summary row per (cell, estimator).
# Bias and Monte Carlo SD are over the finite estimates; the calibration ratio
# mean(SE)/MC-SD separates "is the SE formula right" from "is alpha-hat biased",
# which Wald coverage alone conflates.

summarise_cell <- function(reps_df, truth) {
  by <- split(reps_df, reps_df$estimator)
  rows <- lapply(names(by), function(est) {
    d <- by[[est]]
    fin <- is.finite(d$estimate)
    e <- d$estimate[fin]; se <- d$se[fin]
    covered <- is.finite(d$lwr) & is.finite(d$upr) & d$lwr <= truth & truth <= d$upr
    data.frame(
      estimator      = est,
      truth_alpha    = truth,
      n_rep          = nrow(d),
      n_valid        = sum(fin),
      pct_undefined  = mean(d$undefined %in% TRUE),
      mean_est       = if (length(e)) mean(e) else NA_real_,
      bias           = if (length(e)) mean(e) - truth else NA_real_,
      mc_sd          = if (length(e) > 1) stats::sd(e) else NA_real_,
      mean_se        = if (any(is.finite(se))) mean(se[is.finite(se)]) else NA_real_,
      mean_se_boot   = if (any(is.finite(d$se_boot))) mean(d$se_boot[is.finite(d$se_boot)]) else NA_real_,
      se_sd_ratio    = NA_real_,  # filled below
      coverage       = if (any(is.finite(d$lwr))) mean(covered, na.rm = TRUE) else NA_real_,
      mean_ci_width  = if (any(is.finite(d$upr - d$lwr))) mean((d$upr - d$lwr)[is.finite(d$upr - d$lwr)]) else NA_real_,
      pct_npd        = if (any(is.finite(d$npd))) mean(d$npd[is.finite(d$npd)]) else NA_real_,
      pct_nonfinite_se = mean(!is.finite(d$se)),
      mean_time      = mean(d$time, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$se_sd_ratio <- out$mean_se / out$mc_sd
  out
}
