parameters = expand.grid(
  n = c(10, 40, 100),
  j = c(2, 5, 20),
  weight = c("abs", "quadratic")
)

rho = 0.9
n_reps = 1000

current_time = proc.time() # Use this to find out how long the simulation takes.
cohen_normal_sim = sapply(seq(nrow(parameters)), function(i) {
  sim_ci_normal_model(rho = rho, j = parameters[i, 2], n = parameters[i, 1], type = "cohen",
                      weight = as.character(parameters[i, 3]),n_reps = n_reps)

})
elapsed_time = proc.time() - current_time  # Number of seconds the sim took.

current_time = proc.time() # Use this to find out how long the simulation takes.
fleiss_normal_sim = sapply(seq(nrow(parameters)), function(i) {
  sim_ci_normal_model(rho = rho, j = parameters[i, 2], n = parameters[i, 1], type = "fleiss",
                      weight = as.character(parameters[i, 3]),n_reps = n_reps)

})
elapsed_time = proc.time() - current_time  # Number of seconds the sim took.



# Copied into the document, Cohen's kappa
cohen_normal <- round(rbind(
  cohen_normal_sim[c(1, 4, 2, 5, 3, 6),10:18],
  cohen_normal_sim[c(1, 4, 2, 5, 3, 6),1:9]),2)

# Copied into the document, Fleiss's kappa
fleiss_normal <- round(rbind(
  fleiss_normal_sim[c(1, 4, 2, 5, 3, 6),10:18],
  fleiss_normal_sim[c(1, 4, 2, 5, 3, 6),1:9]),2)

