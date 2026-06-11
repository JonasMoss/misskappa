# Reduce study-32 raw results to one row per (dgp, mech, n, coef).
# Per-dataset spreads across starts, then cell-level aggregates.

summarize_starts <- function(res, ref_start = 0.1) {
  ok <- res[!is.na(res$est), ]
  key <- c("dgp", "mech", "n", "rep")

  # Start-dependent convergence: datasets where some starts error and others
  # succeed (all-fail datasets are design-guard failures, not start effects).
  # Failure rows are one per start; success rows are one per coefficient, so
  # collapse to one row per dataset x start first.
  per_start <- aggregate(list(failed = is.na(res$est)),
                         by = res[, c(key, "start_alpha")], FUN = any)
  conv <- aggregate(list(fail_share = per_start$failed),
                    by = per_start[, key], FUN = mean)
  conv$start_dep_conv <- conv$fail_share > 0 & conv$fail_share < 1
  conv$all_fail <- conv$fail_share == 1
  conv_cell <- aggregate(conv[, c("start_dep_conv", "all_fail")],
                         by = conv[, c("dgp", "mech", "n")], FUN = mean)

  # Spread metrics use datasets where every start converged.
  n_starts <- length(unique(res$start_alpha))
  per <- do.call(rbind, by(ok, ok[, c(key, "coef")], function(d) {
    if (nrow(d) < n_starts) return(NULL)
    ref <- d[d$start_alpha == ref_start, ]
    data.frame(
      dgp = d$dgp[1], mech = d$mech[1], n = d$n[1], rep = d$rep[1],
      coef = d$coef[1],
      spread_est = max(d$est) - min(d$est),
      spread_se_rel = max(d$se) / min(d$se) - 1,
      spread_nf = max(d$null_frac) - min(d$null_frac),
      se_ref = ref$se[1], nf_ref = ref$null_frac[1],
      stringsAsFactors = FALSE
    )
  }))

  per$spread_rel <- per$spread_est / per$se_ref
  per$est_moves <- per$spread_rel > 0.1
  per$se_moves <- per$spread_se_rel > 0.10

  cell <- do.call(rbind, by(per, per[, c("dgp", "mech", "n", "coef")], function(d) {
    movers <- d$est_moves | d$se_moves
    data.frame(
      dgp = d$dgp[1], mech = d$mech[1], n = d$n[1], coef = d$coef[1],
      n_complete = nrow(d),
      spread_rel_med = median(d$spread_rel),
      spread_rel_p90 = unname(quantile(d$spread_rel, 0.9)),
      spread_rel_max = max(d$spread_rel),
      frac_est_moves = mean(d$est_moves),
      frac_se_moves = mean(d$se_moves),
      # Does the null_frac diagnostic flag the datasets the start can move?
      flag_rate_movers = if (any(movers)) mean(d$nf_ref[movers] > 0.01) else NA_real_,
      flag_rate_all = mean(d$nf_ref > 0.01),
      stringsAsFactors = FALSE
    )
  }))

  out <- merge(cell, conv_cell, by = c("dgp", "mech", "n"), all.x = TRUE)
  out[order(out$dgp, out$mech, out$n, out$coef), ]
}
