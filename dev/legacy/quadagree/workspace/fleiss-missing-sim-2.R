# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The Ultimate Diagnostic Toolkit
#
# Part 1: An Alternative Starting Point (Complete Case Analysis)
# Part 2: The Ground Truth Log-Likelihood Function
# Part 3: The Verification Script: Putting our EM to the Test
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# You will definitely need this library
require(numDeriv)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PART 1: A NEW STARTING POINT - COMPLETE CASE
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

get_pi_start_complete_case <- function(data) {
  R <- ncol(data)
  K <- max(data, na.rm = TRUE)
  n_patterns_complete <- K^R

  # Filter for complete cases only
  data_complete <- data[complete.cases(data), ]
  N_complete <- nrow(data_complete)

  if (N_complete == 0) {
    # If no complete cases, start with a uniform distribution
    warning("No complete cases found. Starting EM with a uniform distribution.")
    return(rep(1 / n_patterns_complete, n_patterns_complete))
  }

  # Calculate pi from the frequencies of complete cases
  pattern_counts <- rep(0, n_patterns_complete)
  for (i in 1:N_complete) {
    idx <- pattern_to_index(data_complete[i, ], K)
    pattern_counts[idx] <- pattern_counts[idx] + 1
  }

  # Return the proportions
  return(pattern_counts / N_complete)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PART 2: THE ORACLE - THE OBSERVED LOG-LIKELIHOOD
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This function calculates the log-likelihood directly.
# IMPORTANT: For numerical differentiation, we can't have the sum(pi)=1 constraint.
# So we pass in K^R - 1 parameters, and calculate the last one by subtraction.
# This reparameterization trick makes the Hessian invertible.

observed_log_likelihood <- function(params, data) {
  R <- ncol(data)
  K <- max(data, na.rm = TRUE)
  N <- nrow(data)
  n_patterns_complete <- K^R

  # --- Reparameterization Trick ---
  # Construct the full pi vector from the K^R-1 params
  pi <- numeric(n_patterns_complete)
  pi[1:(n_patterns_complete - 1)] <- params
  pi[n_patterns_complete] <- 1 - sum(params)

  # If any pi are negative or the last one is, this is an invalid
  # parameter set. Return a very large negative number.
  if (any(pi < 0)) return(-1e20)

  # Pre-compute compatibility map for speed (same as in the EM)
  unique_obs_patterns_df <- data %>% as.data.frame() %>% distinct()
  compatibility_list <- list()
  all_complete_patterns_list <- lapply(1:n_patterns_complete, index_to_pattern, R = R, K = K)
  for (i in 1:nrow(unique_obs_patterns_df)) {
    obs_pattern <- as.numeric(unique_obs_patterns_df[i, ])
    key <- paste(obs_pattern, collapse = "-")
    compatible_indices <- which(sapply(all_complete_patterns_list, function(comp_pattern) {
      all(is.na(obs_pattern) | (obs_pattern == comp_pattern))
    }))
    compatibility_list[[key]] <- compatible_indices
  }

  # --- Calculate Log-Likelihood ---
  # Formula: Sum_i log( Sum_j pi_j ) for compatible j
  loglik <- 0
  for (i in 1:N) {
    key <- paste(as.numeric(data[i, ]), collapse = "-")
    compatible_indices <- compatibility_list[[key]]

    # Sum of probabilities of all complete patterns consistent with this observation
    prob_obs <- sum(pi[compatible_indices])

    # If prob is zero, log-lik is -Inf. Return a huge negative number.
    if (prob_obs <= 0) return(-1e20)

    loglik <- loglik + log(prob_obs)
  }

  return(loglik)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PART 3: THE VERIFICATION SCRIPT
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Let's use the sample data from before
set.seed(123)
sample_data <- matrix(c(
  2, 2, 2, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 2,
  2, 1, 1, 1, 2, 2, 2, 1, 2, 1, NA, 1, 2, 2, NA
), ncol = 3, byrow = TRUE)

# --- Step 1: Run our EM to get the MLE ---
# Let's try the new 'complete case' starting point
pi_start_cc <- get_pi_start_complete_case(sample_data)
pi_mle <- estimate_pi_em(sample_data, pi_start = pi_start_cc)

# --- Step 2: Prepare for the Oracle ---
# We need to give numDeriv the K^R-1 vector of parameters
params_mle <- pi_mle[1:(length(pi_mle) - 1)]

# --- Step 3: Call the Oracle (numDeriv) ---
cat("\n--- VERIFICATION PROTOCOL INITIATED ---\n")

# Check 1: Is the gradient at our MLE basically zero?
# If not, our EM didn't truly find the maximum.
cat("\n[CHECK 1] Gradient at EM solution:\n")
grad_numeric <- grad(func = observed_log_likelihood, x = params_mle, data = sample_data)
print(grad_numeric)
cat("...these values should be very close to zero.\n")

# Check 2: Does our Louis's Method Hessian match the numerical Hessian?
# This is the ultimate test of our standard error calculation.
cat("\n[CHECK 2] Hessian Comparison:\n")
# Get our Hessian (for the active, non-zero parameters)
hessian_info <- hessian_louis(pi_mle, sample_data)
I_obs_louis <- hessian_info$I_observed
# Note: The information is the NEGATIVE Hessian
H_louis <- -I_obs_louis

# Get the numerical Hessian (for the K^R-1 subspace)
H_numeric <- hessian(func = observed_log_likelihood, x = params_mle, data = sample_data)

cat("--> Louis's Method Hessian (for active params):\n")
print(H_louis)

cat("\n--> Numerical Hessian from numDeriv (for first K^R-1 params):\n")
print(H_numeric)
cat("...these matrices should be very similar (up to which params are included).\n")

cat("\n--- VERIFICATION PROTOCOL COMPLETE ---\n")
