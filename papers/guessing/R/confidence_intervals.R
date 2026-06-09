#### ===========================================================================
#### Code for confidence intervals in Section 5.
####
#### - The values were manually copied into the paper.
#### ===========================================================================
source("functions.R")

# We first simulate the CIs in the papers.
set.seed(313)
vars <- expand.grid(c(0, 10, 0.5), c(0, 10, 0.5), c(20, 100))
colnames(vars) <- c("guessing", "true", "n")
n_reps <- 10000
uniforms_ <- apply(vars, 1, sim_helper, "uniform", n_reps)
marginals_ <- apply(vars, 1, sim_helper, "marginal", n_reps)

# Then transform the simulated data into a suitable shape.
v <- c(1, 6, 2, 7, 3, 8, 4, 9, 5, 10)

uniforms <- rbind(
  round(uniforms_[v, 1:9], 2),
  round(uniforms_[v, 10:18], 2)
)

marginals <- rbind(
  round(marginals_[v, 1:9], 2),
  round(marginals_[v, 10:18], 2)
)

write.csv(marginals, file = "ci_marginals.csv")
write.csv(uniforms, file = "ci_uniforms.csv")
