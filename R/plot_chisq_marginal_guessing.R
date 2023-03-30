### ===========================================================================
###
### Testing the chisq test for the van Oest model.
###
### ===========================================================================
source("R/functions.R")

set.seed(313)
results_mgdm <- replicate(10000, {
  n <- 100
  s <- 0.7
  true_dist <- c(0.2, 0.1, 0.1, 0.3, 0.1)
  x <- agreeable::simulate_jsm(n, c(s, s), model = "fleiss", true_dist = true_dist)
  tab <- Rfast::Table(x = x[, 1], y = x[, 2], names = FALSE)

  mgdm_test_2(tab)
})

set.seed(313)
results_plm <- replicate(10000, {
  n <- 100
  s <- 0.7
  true_dist <- c(0.2, 0.1, 0.1, 0.3, 0.1)
  x <- agreeable::simulate_jsm(n, c(s, s), model = "bp", true_dist = true_dist)
  tab <- Rfast::Table(x = x[, 1], y = x[, 2], names = FALSE)

  perreault_leigh_test(tab)
})


xx <- seq(0, 100)
n_cat <- 5

par(mfrow = c(1, 2))
hist(results_plm[1, ],
  freq = FALSE, breaks = 100, main = "Perreault-Leigh",
  xlab = "Statistic", ylab = "Density"
)
lines(xx, dchisq(xx, n_cat^2 - n_cat - 2))

hist(results_mgdm[1, ],
  freq = FALSE, breaks = 100, main = "Marginal guessing distribution",
  xlab = "Statistic", ylab = "Density"
)
lines(xx, dchisq(xx, n_cat^2 - n_cat - 2))
