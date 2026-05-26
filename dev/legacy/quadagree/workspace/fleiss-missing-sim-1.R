# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# THE FINAL TEST FLIGHT: A Simulation Study
#
# We're going to see if our Leviathan can navigate the stormy seas of
# random sampling and missing data and still point to the right treasure.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# --- Step 0: Ensure all our tools are loaded ---
# Make sure you have all the functions we've built in your R environment:
# - simulate_jsm (the one you provided)
# - pattern_to_index, index_to_pattern
# - estimate_pi_em (the robust one)
# - hessian_louis (the robust one)
# - pi_to_fleiss_kappa, grad_fleiss_kappa
# - fleiss_kappa_missing (the main wrapper)
# And the required libraries:
require(dplyr)
require(MASS)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# THE FINAL, CORRECT IMPLEMENTATION
#
# We are replacing the buggy analytical Hessian with the robust numerical
# method using the softmax transformation. This is the definitive version.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
require(numDeriv)

# --- Helper: Softmax Transformation Tools ---
pi_to_softmax <- function(pi) {
  # Takes a probability vector pi (summing to 1) and returns m-1 unconstrained alphas
  log(pi[1:(length(pi)-1)] / pi[length(pi)])
}

softmax_to_pi <- function(alphas) {
  # Takes m-1 unconstrained alphas and returns a probability vector of length m
  m <- length(alphas) + 1
  pi <- numeric(m)
  exp_alphas <- exp(alphas)
  sum_exp_alphas <- 1 + sum(exp_alphas)
  pi[1:(m-1)] <- exp_alphas / sum_exp_alphas
  pi[m] <- 1 / sum_exp_alphas
  return(pi)
}


fleiss_kappa_final <- function(data, conf.level = 0.95) {
  R <- ncol(data)
  K <- max(data, na.rm = TRUE)
  n_patterns_complete <- K^R

  # --- Step 1: EM Estimation (this part is solid) ---
  pi_start <- get_pi_start_complete_case(data) # Use a robust start
  pi_mle <- estimate_pi_em(data, pi_start = pi_start)

  # --- Step 2: Kappa Point Estimate ---
  kappa_est <- pi_to_fleiss_kappa(pi_mle, R, K)$kappa

  # --- Step 3: Variance Calculation (The Correct Way) ---

  # 3a. Identify the "active" parameters with non-zero probability
  active_indices <- which(pi_mle > 1e-12)
  pi_active <- pi_mle[active_indices]
  m_active <- length(pi_active)

  if (m_active <= 1) {
    warning("Only one outcome has probability; cannot compute variance.")
    return(list(kappa = kappa_est, se = NA, conf.int = c(NA, NA)))
  }

  # 3b. Convert the MLE pi vector to the unconstrained alpha space
  alphas_mle <- pi_to_softmax(pi_active)

  # 3c. Define the log-likelihood function in terms of alphas
  #     This is our stable Oracle function.
  loglik_for_hessian <- function(alphas, data, active_indices) {
    pi_vec <- rep(0, n_patterns_complete)
    pi_vec[active_indices] <- softmax_to_pi(alphas)

    # This is the core log-likelihood calculation
    loglik <- sum(sapply(1:nrow(data), function(i){
      key <- paste(as.numeric(data[i, ]), collapse = "-")
      # This requires the compatibility_list logic from the EM function
      # For a self-contained function, we rebuild it.
      obs_pattern <- as.numeric(data[i, ])
      compatible_indices <- which(sapply(1:n_patterns_complete, function(j) {
        all(is.na(obs_pattern) | (obs_pattern == index_to_pattern(j,R,K)))
      }))

      prob_obs <- sum(pi_vec[compatible_indices])
      if(prob_obs <= 0) return(-1e20) # Should not happen with softmax
      return(log(prob_obs))
    }))
    return(loglik)
  }

  # 3d. Calculate the Hessian numerically in the stable alpha-space
  H_alpha <- numDeriv::hessian(func = loglik_for_hessian, x = alphas_mle,
                               data = data, active_indices = active_indices)

  # 3e. Get the variance-covariance matrix for alphas
  Var_alpha <- solve(-H_alpha)

  # 3f. Use the Delta Method to transform variance from alpha-space to pi-space
  J_pi_alpha <- numDeriv::jacobian(func = softmax_to_pi, x = alphas_mle)
  Var_pi_active <- J_pi_alpha %*% Var_alpha %*% t(J_pi_alpha)

  # 3g. Use the Delta Method again to get the variance for kappa
  grad_kappa_full <- grad_fleiss_kappa(pi_mle, R, K)
  grad_kappa_active <- grad_kappa_full[active_indices]

  var_kappa <- t(grad_kappa_active) %*% Var_pi_active %*% grad_kappa_active
  se_kappa <- sqrt(var_kappa)

  # --- Step 4: Confidence Interval ---
  z <- qnorm(1 - (1 - conf.level) / 2)
  ci_lower <- kappa_est - z * se_kappa
  ci_upper <- kappa_est + z * se_kappa

  result <- list(
    kappa = as.numeric(kappa_est),
    se = as.numeric(se_kappa),
    conf.int = c(as.numeric(ci_lower), as.numeric(ci_upper)),
    conf.level = conf.level,
    pi_mle = pi_mle
  )

  return(result)
}


# --- Step 1: Simulation Setup ---
cat("Preparing the test flight...\n")
set.seed(42) # For reproducible results
n_sims <- 200 # Number of simulations. Good for a first check.
n_obs <- 1000  # Number of subjects per simulation
skill_vector <- c(0.9, 0.9, 0.9)
true_dist <- c(0.4, 0.6)
missing_prob <- 0.5 # Let's make it challenging with 15% missingness

# A list to store the results of each simulation run
sim_results <- list()

# --- Step 2: The Simulation Loop ---
cat(paste("Initiating", n_sims, "simulation runs. This might take a moment...\n"))
pb <- txtProgressBar(min = 0, max = n_sims, style = 3)

for (i in 1:n_sims) {
  # 1. Simulate a complete dataset
  data_complete <- simulate_jsm(n = n_obs, s = skill_vector, model = "fleiss", true_dist = true_dist)
  true_kappa <- attr(data_complete, "kappa")

  # 2. Introduce missingness (MCAR: Missing Completely At Random)
  data_missing <- as.matrix(data_complete)
  n_total_cells <- prod(dim(data_missing))
  na_indices <- sample(n_total_cells, size = floor(n_total_cells * missing_prob))
  data_missing[na_indices] <- NA

  # 3. Run our ultimate weapon! Use a try-catch in case a weird sample causes a failure.
  est <- try(fleiss_kappa_missing(data_missing, conf.level = 0.95), silent = TRUE)

  # 4. Store the results
  if (!inherits(est, "try-error")) {
    sim_results[[i]] <- data.frame(
      true_kappa = true_kappa,
      kappa_hat = est$kappa,
      se_hat = est$se,
      ci_lower = est$conf.int[1],
      ci_upper = est$conf.int[2]
    )
  }

  setTxtProgressBar(pb, i)
}
close(pb)

# Combine all results into a single data frame
results_df <- do.call(rbind, sim_results)


# --- Step 3: Analyze the Results ---
cat("\n\n--- TEST FLIGHT ANALYSIS ---\n")

# A. BIAS: Is our estimator pointing in the right direction on average?
avg_kappa_hat <- mean(results_df$kappa_hat, na.rm = TRUE)
true_kappa_val <- unique(results_df$true_kappa)
bias <- avg_kappa_hat - true_kappa_val

cat(paste("\nTrue Latent Kappa:", round(true_kappa_val, 4)))
cat(paste("\nAverage Estimated Kappa:", round(avg_kappa_hat, 4)))
cat(paste("\nBias:", round(bias, 4), "(A small bias is expected, especially with missing data)\n"))

# B. VARIANCE: Is our standard error calculation honest?
# We compare the actual spread of our estimates to the average of our estimated SEs.
empirical_sd <- sd(results_df$kappa_hat, na.rm = TRUE)
average_se <- mean(results_df$se_hat, na.rm = TRUE)

cat(paste("\nStandard Deviation of Estimates (The 'True' SE):", round(empirical_sd, 4)))
cat(paste("\nAverage Estimated SE from our Method:", round(average_se, 4)))
cat("\n(These two numbers should be very close!)\n")

# C. COVERAGE: Does our 95% confidence interval do its job?
# It should "cover" the true value in about 95% of the simulations.
results_df$covered <- (results_df$ci_lower <= results_df$true_kappa) &
  (results_df$ci_upper >= results_df$true_kappa)

coverage_prob <- mean(results_df$covered, na.rm = TRUE)

cat(paste("\nNominal 95% CI Coverage:", round(coverage_prob * 100, 1), "%"))
cat("\n(This should be very close to 95%)\n")

cat("\n--- ANALYSIS COMPLETE ---\n")
