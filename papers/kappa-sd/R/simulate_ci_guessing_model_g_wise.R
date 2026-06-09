### We will use `sim_ci_guessing_model_g_wise` in this file, hence we have to
### source the function `simulate_ci_functions_g_wise.R`.
### Running this file took a little over an hour on my 2 kernel 2 ghz computer.

source("R/simulate_ci_functions_g_wise.R")
source("R/simulate_jsm_true_value.R")

### Parameters
r <- 5
true_dist <- rep(1, r) / r
s <- rep(sqrt(0.8), r)
model <- "bp"

library("future.apply")
plan(multisession, workers = availableCores() - 2)

### Helper function
sim_g_wise <- function(parameters, n_reps, r, type = "cohen") {
  future.apply::future_apply(parameters, 1, function(row) {
    n <- as.numeric(row[1])
    g <- as.numeric(row[2])
    disagreement <- row[3]
    true_value <- as.numeric(row[4])

    sim_ci_guessing_model_g_wise(
      s = rep(sqrt(0.8), r),
      true_dist = true_dist,
      model = model,
      n = n,
      type = type,
      disagreement = disagreement,
      g = g,
      true_value = true_value,
      n_reps = n_reps
    )
  }, future.seed = NULL)
}

### The guessing models involve 3 ** 2 * 4 = 36 parameters, neatly tucked inside
### a data frame.
if (!exists("sim_jsm")) {
  sim_jsm <- memoise::memoise(sim_jsm_true_value)
}

n <- 1000000
parameters <- expand.grid(
  n = c(10, 40, 100),
  g = c(3, 4, 5),
  disagreement = c("nominal", "absolute", "quadratic", "hubert")
)

trues <- apply(parameters, 1, function(x) {
  sim_jsm(n, r, true_dist, s, model, as.numeric(x[2]), x[3])
})

parameters <- cbind(parameters, trues)

### This is the Cohen guessing model CI simulation.
set.seed(313)
n_reps <- 1000
current_time <- proc.time() # Use this to find out how long the simulation takes.
cohen_sim_g_wise <- sim_g_wise(parameters, n_reps, r, type = "cohen")
elapsed_time <- proc.time() - current_time # Number of seconds the sim took.

cohens = cbind(parameters, t(cohen_sim_g_wise))

cohens_g_wise_out <- round(rbind(
  cohen_sim_g_wise[c(1, 3, 2, 4),1:9],
  cohen_sim_g_wise[c(1, 3, 2, 4),10:18],
  cohen_sim_g_wise[c(1, 3, 2, 4),19:27],
  cohen_sim_g_wise[c(1, 3, 2, 4),28:36]),2)

### This is the Fleiss guessing model CI simulation.
set.seed(313)
n_reps = 1000
current_time = proc.time() # Use this to find out how long the simulation takes.
fleiss_sim_g_wise = sim_g_wise(parameters, n_reps, r, type = "fleiss")
elapsed_time = proc.time() - current_time  # Number of seconds the sim took.

cbind(parameters, t(fleiss_sim_g_wise))

fleiss_g_wise_out <- round(rbind(
  fleiss_sim_g_wise[c(1, 3, 2, 4),1:9],
  fleiss_sim_g_wise[c(1, 3, 2, 4),10:18],
  fleiss_sim_g_wise[c(1, 3, 2, 4),19:27],
  fleiss_sim_g_wise[c(1, 3, 2, 4),28:36]),2)

rbind(
  fleiss_sim_g_wise[c(1, 3, 2, 4),1:9],
  fleiss_sim_g_wise[c(1, 3, 2, 4),10:18],
  fleiss_sim_g_wise[c(1, 3, 2, 4),19:27],
  fleiss_sim_g_wise[c(1, 3, 2, 4),28:36])
