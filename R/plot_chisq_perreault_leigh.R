### ===========================================================================
###
### Testing and experimentation with Chisquare test for the
###   Perreault--Leigh model.
###
### ===========================================================================
source("R/functions.R")

set.seed(313)
results <- replicate(10000, {
  n <- 100
  s <- 0.7
  true_dist <- c(0.2, 0.1, 0.1, 0.3, 0.1)
  x <- agreeable::simulate_jsm(n, c(s, s), model = "bp", true_dist = true_dist)
  tab <- Rfast::Table(x = x[, 1], y = x[, 2], names = FALSE)

  perreault_leigh_test(tab)
})

hist(results[1, ],
  freq = FALSE, breaks = 100, main = "Perreault-Leigh model",
  xlab = "Statistic", ylab = "Density"
)
xx <- seq(0, 50)
n_cat <- 5
lines(xx, dchisq(xx, n_cat^2 - n_cat - 2))
