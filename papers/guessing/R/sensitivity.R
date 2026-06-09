#### ===========================================================================
#### Code for the sensitivity analysis of Section 4.
####
#### - We run 8 simulations. 6 with dependence, 2 with independence.
#### - We export the data to a .csv file.
#### - The values from the .csv file were manually copied into the paper.
#### ===========================================================================
source("functions.R")

set.seed(313)
vars <- expand.grid(c(0, 10, 0.5), c(0, 10, 0.5))
n_reps <- 1000

sensitivies_ <- list(
  apply(vars, 1, sens_all_gcb, n_reps, "uniform", rho = 0.9),
  apply(vars, 1, sens_all_gcb, n_reps, "marginal", rho = 0.9),
  apply(vars, 1, sens_all_gcb, n_reps, "uniform", rho = 0.5),
  apply(vars, 1, sens_all_gcb, n_reps, "marginal", rho = 0.5),
  apply(vars, 1, sens_all_gcb, n_reps, "uniform", rho = 0.2),
  apply(vars, 1, sens_all_gcb, n_reps, "marginal", rho = 0.2),
  apply(vars, 1, sens_all_jsm, n_reps, "uniform"),
  apply(vars, 1, sens_all_jsm, n_reps, "marginal")
)

sensitivies <- lapply(sensitivies_, formater)
write.csv(do.call(rbind, sensitivies), file = "sensitivites.csv")
