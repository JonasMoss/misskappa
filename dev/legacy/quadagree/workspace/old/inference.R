#' @keywords internal
quadagree_internal <- function(calc,
                               transform,
                               conf_level,
                               alternative = c("two.sided", "greater", "less"),
                               bootstrap,
                               n_reps,
                               call,
                               fun,
                               ...) {
  alternative <- match.arg(alternative)
  transformer <- get_transformer(transform)
  quants <- limits(alternative, conf_level)
  est_var <- fun(calc)
  est <- est_var$est
  sd <- sqrt(max(est_var$var, 0))

  out <- if (!bootstrap) {
    ci_asymptotic(est, sd, calc$n, transformer, quants)
  } else {
    ci_boot(calc, fun, transformer, quants, n_reps)
  }

  out[2] <- min(out[2], 1)
  names(out) <- quants
  structure(out,
    conf_level = conf_level,
    alternative = alternative,
    type = calc$type,
    n = calc$n,
    transform = transform,
    bootstrap = bootstrap,
    n_reps = n_reps,
    estimate = est,
    sd = sd,
    call = call,
    class = "quadagree"
  )
}

#' @keywords internal
ci_asymptotic <- function(est, sd, n, transformer, quants) {
  est_t <- transformer$est(est)
  sd_t <- transformer$sd(est, sd)
  multiplier <- stats::qt(quants, n - 1) / sqrt(n - 1)
  sort(transformer$inv(est_t + multiplier * sd_t))
}

#' @keywords internal
ci_boot <- function(calc, fun, transformer, quants, n_reps) {
  est_var <- fun(calc)
  est_t <- transformer$est(est_var$est)
  sd_t <- transformer$sd(est_var$est, sqrt(est_var$var))
  boots <- bootstrapper(calc, est_t, fun, transformer, n_reps)
  multiplier <- stats::quantile(boots, quants, na.rm = TRUE)
  sort(transformer$inv(est_t + multiplier * sd_t))
}

#' @keywords internal
bootstrapper <- function(calc, est_t, fun, transformer, n_reps) {
  calc_new <- calc
  trans_est <- transformer$est
  trans_sd <- transformer$sd
  future.apply::future_replicate(n_reps,
    {
      indices_star <- sample.int(calc$n, replace = TRUE)
      calc_new$xx <- calc$xx[indices_star, ]
      est_var <- fun(calc_new)
      est_star <- est_var$est
      sd_star <- sqrt(est_var$var)
      (trans_est(est_star) - est_t) / trans_sd(est_star, sd_star)
    },
    future.seed = TRUE
  )
}
