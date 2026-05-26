# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The Core Engine: A Robust EM for the MLE of Pi
#
# This function's only job is to find the Maximum Likelihood Estimate
# of the K^R pattern probabilities. It is the solid foundation for
# everything else.
#
# It incorporates our best practices:
# - A robust default starting point (complete case analysis).
# - A pre-computed compatibility map for speed.
# - Clean, debugged E-step and M-step logic.
# - Returns useful diagnostics like convergence status and iterations.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# --- Prerequisite Helper Functions ---
# (Including them here so the code block is self-contained)

get_pi_start_complete_case <- function(data) {
  R <- ncol(data)
  K <- max(data, na.rm = TRUE)
  n_patterns_complete <- K^R
  data_complete <- data[complete.cases(data), ]
  N_complete <- nrow(data_complete)
  if (N_complete == 0) {
    warning("No complete cases found. Starting EM with a uniform distribution.")
    return(rep(1 / n_patterns_complete, n_patterns_complete))
  }
  pattern_counts <- rep(0, n_patterns_complete)
  for (i in 1:N_complete) {
    idx <- pattern_to_index(data_complete[i, ], K)
    pattern_counts[idx] <- pattern_counts[idx] + 1
  }
  return(pattern_counts / N_complete)
}

pattern_to_index <- function(pattern, K) {
  R <- length(pattern)
  idx <- 0
  for (r in 1:R) {
    idx <- idx + (pattern[r] - 1) * (K^(r - 1))
  }
  return(idx + 1)
}

index_to_pattern <- function(idx, R, K) {
  pattern <- numeric(R)
  temp_idx <- idx - 1
  for (r in 1:R) {
    pattern[r] <- (temp_idx %% K) + 1
    temp_idx <- floor(temp_idx / K)
  }
  return(pattern)
}

# --- THE MAIN FUNCTION ---

#' Estimate Multinomial Probabilities via EM for Missing Data
#'
#' @description
#' This function calculates the Maximum Likelihood Estimate (MLE) of the probabilities
#' for all K^R possible rating patterns from a data matrix with missing values (NAs).
#' It assumes the data are Missing at Random (MAR) or Missing Completely at Random (MCAR).
#'
#' @param data An N x R matrix or data frame of categorical ratings.
#' @param pi_start Optional starting vector for the probabilities. If NULL, a
#'   complete-case analysis is used to generate a smart starting point.
#' @param max_iter Maximum number of EM iterations.
#' @param tol Convergence tolerance. The algorithm stops when the maximum
#'   absolute change in any pi parameter is less than this value.
#'
#' @return A list containing:
#' - `pi_mle`: The final MLE vector of probabilities (length K^R).
#' - `iterations`: The number of iterations performed.
#' - `converged`: A logical indicating if the algorithm converged.
#' - `loglik`: The final log-likelihood value at the MLE.
#'
get_pi_mle <- function(data, pi_start = NULL, max_iter = 1000, tol = 1e-9) {
  N <- nrow(data)
  R <- ncol(data)
  K <- max(data, na.rm = TRUE)
  n_patterns_complete <- K^R

  # --- Step 1: Get Starting Values ---
  if (is.null(pi_start)) {
    pi_current <- get_pi_start_complete_case(data)
  } else {
    pi_current <- pi_start
  }

  # --- Step 2: Pre-compute Compatibility Map for speed ---
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

  # --- Step 3: The Main EM Loop ---
  converged <- FALSE
  for (iter in 1:max_iter) {
    # --- E-Step: Calculate expected counts `n_hat` ---
    n_hat <- rep(0, n_patterns_complete)
    for (i in 1:N) {
      key <- paste(as.numeric(data[i, ]), collapse = "-")
      compatible_indices <- compatibility_list[[key]]
      if (length(compatible_indices) == 0) next

      pi_compatible <- pi_current[compatible_indices]
      prob_obs <- sum(pi_compatible)

      if (prob_obs > 1e-12) {
        n_hat[compatible_indices] <- n_hat[compatible_indices] + (pi_compatible / prob_obs)
      }
    }

    # --- M-Step: Update pi ---
    pi_new <- n_hat / N

    # --- Check Convergence ---
    if (max(abs(pi_new - pi_current)) < tol) {
      converged <- TRUE
      pi_current <- pi_new
      break
    }
    pi_current <- pi_new
  }

  # --- Step 4: Finalize and Calculate Log-Likelihood ---
  pi_mle <- pi_current

  final_loglik <- sum(sapply(1:N, function(i) {
    key <- paste(as.numeric(data[i, ]), collapse = "-")
    compatible_indices <- compatibility_list[[key]]
    prob_obs <- sum(pi_mle[compatible_indices])
    return(log(prob_obs))
  }))

  if (!converged) {
    warning(paste("EM did not converge within", max_iter, "iterations."))
  }

  return(list(
    pi_mle = pi_mle,
    iterations = iter,
    converged = converged,
    loglik = final_loglik
  ))
}

# --- EXAMPLE USAGE ---
set.seed(42)
sample_data <- matrix(c(
  2, 2, 2, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 2,
  2, 1, 1, 1, 2, 2, 2, 1, 2, 1, NA, 1, 2, 2, NA
), ncol = 3, byrow = TRUE)

# Run the function
mle_results <- get_pi_mle(sample_data)

# See the beautiful results
print(paste("Converged:", mle_results$converged, "in", mle_results$iterations, "iterations."))
print("Final MLE for pi:")
print(round(mle_results$pi_mle, 4))
