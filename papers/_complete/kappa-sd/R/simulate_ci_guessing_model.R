### We will use `sim_ci_guessing_model` in this file, hence we have to
### source the function `simulate_ci_functions.R`.
source("R/simulate_ci_functions.R")

### The guessing models involve 3 ** 3 = 27 parameters, neatly tucked inside
### this data frame.

parameters = expand.grid(
  n = c(10, 40, 100),
  j = c(2, 5, 20),
  weight = c("unit", "abs", "quadratic")
)

### This is the Cohen guessing model CI simulation.

n_reps = 10000
current_time = proc.time() # Use this to find out how long the simulation takes.
cohen_guessing_sim = apply(parameters, 1, function(row) {

  n = row[1]
  j = row[2]
  weight = row[3]

  sim_ci_guessing_model(
    s = rep(sqrt(0.8), j),
    true_dist = rep(1,5) / 5,
    model = "bp",
    n = n,
    type = "cohen",
    weight = weight,
    n_reps = n_reps)

})
elapsed_time = proc.time() - current_time  # Number of seconds the sim took.

### This is the Fleiss guessing model CI simulation.
current_time = proc.time() # Use this to find out how long the simulation takes.
fleiss_guessing_sim = apply(parameters, 1, function(row) {

  n = row[1]
  j = row[2]
  weight = row[3]

  sim_ci_guessing_model(
    s = rep(sqrt(0.8), j),
    true_dist = rep(1,5) / 5,
    model = "bp",
    n = n,
    type = "fleiss",
    weight = weight,
    n_reps = n_reps)

})
elapsed_time = proc.time() - current_time  # Number of seconds the sim took.

### These transformations make the simulation data easier to understand.
t(cbind(parameters, t(cohen_guessing_sim)))
t(cbind(parameters, t(fleiss_guessing_sim)))


# Copied into the document, Cohen's kappa
coehn_guessing_sim <- round(rbind(
  cohen_guessing_sim[c(1, 4, 2, 5, 3, 6),1:9],
  cohen_guessing_sim[c(1, 4, 2, 5, 3, 6),19:27],
  cohen_guessing_sim[c(1, 4, 2, 5, 3, 6),10:18]),2)

# Copied into the document, Fleiss's kappa
fleiss_guessing_sim <- round(rbind(
  fleiss_guessing_sim[c(1, 4, 2, 5, 3, 6),1:9],
  fleiss_guessing_sim[c(1, 4, 2, 5, 3, 6),19:27],
  fleiss_guessing_sim[c(1, 4, 2, 5, 3, 6),10:18]),2)