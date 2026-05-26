out <- kappa_raw(x_missing, transform = "none", bootstrap = FALSE)
out$fleiss$conf.low  # should now match the R version
out$fleiss$sd        # plug this into ci_from_est_sd

ci_from_est_sd <- function(estimate, sd, conf_level = 0.95) {
  z <- qnorm(1 - (1 - conf_level)/2)
  c(estimate - z * sd, estimate + z * sd)
}


ci_from_est_sd(out$fleiss$est, out$fleiss$sd)

out

x <- seq(0.0001, 0.999, length.out = 100)
abs(sapply(x, TestQNormCpp) - qnorm(x))
