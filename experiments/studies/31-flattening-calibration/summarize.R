# Reduce study-31 raw results to one row per (dgp, mech, n, flatten, coef).
# Shared by run_calibration.R and ad-hoc reanalysis.

summarize_calibration <- function(res) {
  res$failed <- is.na(res$est)
  key <- c("dgp", "mech", "n", "flatten")

  fail <- aggregate(failed ~ dgp + mech + n + flatten, data = res,
                    FUN = function(z) mean(z))
  # Failures have coef = NA; per-coefficient stats use successful fits only.
  ok <- res[!res$failed, ]
  stats <- do.call(rbind, by(ok, ok[, c(key, "coef")], function(d) {
    err <- d$est - d$truth
    covered <- abs(err) <= qnorm(0.975) * d$se
    data.frame(
      dgp = d$dgp[1], mech = d$mech[1], n = d$n[1], flatten = d$flatten[1],
      coef = d$coef[1], n_ok = nrow(d),
      bias = mean(err), sd = sd(d$est), rmse = sqrt(mean(err^2)),
      se_over_sd = mean(d$se, na.rm = TRUE) / sd(d$est),
      cover95 = mean(covered, na.rm = TRUE),
      null_frac_mean = mean(d$null_frac, na.rm = TRUE),
      null_frac_p90 = unname(quantile(d$null_frac, 0.9, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  }))
  out <- merge(stats, fail, by = key, all.x = TRUE)
  names(out)[names(out) == "failed"] <- "fail_rate"
  out[order(out$dgp, out$mech, out$n, out$coef, out$flatten), ]
}
