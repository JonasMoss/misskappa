#' Calculate Observed Data Log-Likelihood
#'
#' @param x The incomplete data matrix.
#' @param phi The current phi parameter vector.
#' @param constants Pre-computed model constants.
#' @return The total log-likelihood of the observed data.
calculate_log_likelihood <- function(x, phi, constants) {
  # This function's logic is very similar to the E-step's inner loop
  n <- nrow(x)
  R <- ncol(x)
  C <- nrow(constants$K_mat)

  # This setup is identical to the E-step and could be pre-computed once
  C_k <- apply(constants$K_mat, 2, function(k) exp(lfactorial(R) - sum(lfactorial(k))))
  theta_k <- phi / C_k
  theta_k[C_k == 0] <- 0

  all_full_vectors <- t(gtools::permutations(R, C, 1:C, repeats.allowed = TRUE))
  full_vec_counts <- apply(all_full_vectors, 2, function(v) tabulate(v, nbins = C))

  all_patterns_list <- as.list(as.data.frame(constants$K_mat))
  map_to_k_idx <- vapply(as.list(as.data.frame(full_vec_counts)), function(p_counts) {
    which(vapply(all_patterns_list, function(k) all(k == p_counts), logical(1)))
  }, integer(1))

  pi_z_full <- theta_k[map_to_k_idx]

  total_log_lik <- 0
  for (i in 1:n) {
    obs_data_i <- x[i, ]
    is_missing <- is.na(obs_data_i)

    consistent_z_mask <- apply(all_full_vectors, 2, function(z) {
      all(z[!is_missing] == obs_data_i[!is_missing])
    })

    lik_i <- sum(pi_z_full[consistent_z_mask])
    if (lik_i > 0) {
      total_log_lik <- total_log_lik + log(lik_i)
    }
  }
  return(total_log_lik)
}

#' Helper to calculate Kappa from Phi (to avoid code duplication)
calculate_kappa_from_phi <- function(phi, R, C, loss_matrix, constants) {
  K_mat <- constants$K_mat
  d_vec <- constants$d_vec

  p_estimate <- (1 / R) * (K_mat %*% phi)
  observed_disagreement <- sum(d_vec * phi)
  chance_disagreement <- t(p_estimate) %*% loss_matrix %*% p_estimate

  kappa <- 1 - observed_disagreement / as.numeric(chance_disagreement)

  return(list(
    kappa = kappa,
    observed_agreement = 1 - observed_disagreement,
    chance_agreement = 1 - as.numeric(chance_disagreement)
  ))
}

#' Estimate Model-Based Kappa from Incomplete Count Data via EM
#'
#' Estimates weighted Fleiss' Kappa from count data that may be incomplete
#' (i.e., contain missing ratings). Missingness is inferred when a row's sum
#' is less than the total number of raters, R.
#'
#' @param x A data matrix where rows are subjects and columns are categories.
#'        Each cell x_ic contains the number of raters who assigned category c
#'        to subject i.
#' @param R Integer. The total number of raters in the study. If NULL, it is
#'        inferred from the row with the maximum sum.
#' @param loss_matrix A C x C disagreement/loss matrix. Defaults to unweighted.
#' @param alpha A small numeric value for additive smoothing.
#' @param tol Convergence tolerance for the EM algorithm.
#' @param max_iter Maximum number of iterations for the EM algorithm.
#' @return A list containing the kappa estimate and other model results.
estimate_kappa_em <- function(x, R = NULL, loss_matrix = NULL, alpha = 1e-9,
                              tol = 1e-8, max_iter = 500) {

  x <- as.matrix(x)
  # --- 1. Setup and Initialization ---
  n <- nrow(x)
  C <- ncol(x)

  if (is.null(R)) {
    R <- max(rowSums(x))
    message(paste("R not provided. Inferred R =", R, "from data."))
  }

  if (any(rowSums(x) > R)) {
    stop("A row in the data has more ratings than the total number of raters R.")
  }

  if (is.null(loss_matrix)) {
    loss_matrix <- matrix(1, nrow = C, ncol = C) - diag(C)
  }

  constants <- precompute_kappa_constants(R, C, loss_matrix)
  K_mat <- constants$K_mat
  K_patterns <- ncol(K_mat)

  # --- 2. Initial Parameter Guess (phi_t) ---
  phi_t <- rep(1 / K_patterns, K_patterns)

  # --- 3. EM Iteration Loop ---
  log_likelihood_old <- -Inf

  for (iter in 1:max_iter) {
    # --- E-Step: Calculate Expected Complete Counts ---
    # This is now much simpler.

    # Pre-calculate probabilities of all full count patterns, C(k)*theta_k
    C_k <- apply(K_mat, 2, function(k) exp(lfactorial(R) - sum(lfactorial(k))))

    # Check for phi_t being exactly zero, which can happen with alpha=0
    # In that case, theta_k should also be zero.
    theta_k <- ifelse(C_k > 0, phi_t / C_k, 0)

    expected_n_k <- numeric(K_patterns)

    # Loop over each subject (row in the count data)
    for (i in 1:n) {
      observed_counts_i <- x[i, ]
      n_observed_i <- sum(observed_counts_i)
      n_missing_i <- R - n_observed_i

      if (n_missing_i == 0) {
        # Subject has complete data, no imputation needed.
        # Find which pattern k this corresponds to.
        k_idx <- which(apply(K_mat, 2, function(k) all(k == observed_counts_i)))
        expected_n_k[k_idx] <- expected_n_k[k_idx] + 1
      } else {
        # Subject has missing data, perform imputation.

        # Find all complete patterns k that are 'consistent' with the observed counts.
        # A pattern k is consistent if k_c >= x_ic for all categories c.
        is_consistent <- apply(K_mat, 2, function(k) all(k >= observed_counts_i))
        consistent_k_indices <- which(is_consistent)

        # Calculate the 'imputation weights' for this subject.
        # The weight is proportional to the probability of seeing the missing part,
        # given the full pattern k. This is a multinomial probability.

        # Prob(observed part | full pattern k) is proportional to C(k - obs)
        # Prob(full pattern k) is phi_t[k]
        # Weight w_k is proportional to phi_t[k] * Prob(obs | k) -> which is complex
        # A simpler path: Prob(k | obs) is proportional to Prob(k) * Prob(obs | k)
        # Using the generative model: P(k | obs) proportional to theta_k * C(k - obs)

        imputation_weights <- sapply(consistent_k_indices, function(k_idx) {
          k_vec <- K_mat[, k_idx]
          # This is theta_k * C(k|obs)
          lfactorial(n_missing_i) - sum(lfactorial(k_vec - observed_counts_i)) + log(theta_k[k_idx])
        })

        # Normalize weights to sum to 1 (after converting from log-scale)
        max_log_weight <- max(imputation_weights)
        norm_weights <- exp(imputation_weights - max_log_weight)
        norm_weights <- norm_weights / sum(norm_weights)

        # Add this subject's 'soft' counts to the total expected counts
        expected_n_k[consistent_k_indices] <- expected_n_k[consistent_k_indices] + norm_weights
      }
    }

    # --- M-Step: Update Phi ---
    phi_t1 <- (expected_n_k + alpha) / (n + K_patterns * alpha)

    # --- Check for Convergence ---
    # We can monitor the change in the phi vector itself for simplicity
    change <- sqrt(sum((phi_t1 - phi_t)^2))
    if (change < tol) {
      phi_t <- phi_t1
      # cat(sprintf("EM converged in %d iterations.\n", iter))
      break
    }

    phi_t <- phi_t1

    if (iter == max_iter) {
      warning("EM algorithm did not converge within max_iter.")
    }
  }

  # --- 4. Final Kappa Calculation ---
  kappa_results <- calculate_kappa_from_phi(phi_t, R, C, loss_matrix, constants)

  return(list(
    kappa = kappa_results$kappa,
    phi_estimate = phi_t,
    observed_agreement = kappa_results$observed_agreement,
    chance_agreement = kappa_results$chance_agreement,
    params = list(n=n, R=R, C=C, alpha=alpha)
  ))
}

# Helper functions precompute_kappa_constants and calculate_kappa_from_phi
# are assumed to be loaded from our previous correct versions.
# --- Verification Script for EM with Count Data ---

# Assume dat.fleiss1971 is loaded
# It's a 30x5 matrix of counts. R=6.
complete_data <- as.matrix(dat.fleiss1971)

# 1. Create incomplete data from the counts
# Let's remove 1 rating from the first subject and 2 from the second.
incomplete_data <- complete_data
# Subject 1: was (0,0,0,6,0). Let's say one 'neurosis' rating is missing.
incomplete_data[1, ] <- c(0, 0, 0, 5, 0) # Now sums to 5
# Subject 2: was (0,3,0,0,3). Let's say one 'pers dis' and one 'other' are missing.
incomplete_data[2, ] <- c(0, 2, 0, 0, 2) # Now sums to 4

# 2. Run our estimators
# Use the original complete-data function (kappa_model) on the original data
kappa_complete <- kappa_model(complete_data)$kappa

# Use our new EM function on the incomplete data
kappa_em <- estimate_kappa_em(incomplete_data, R = 6)$kappa


# 3. Compare results
cat("--- EM Verification with Incomplete Count Data ---\n\n")
cat(sprintf("Kappa from original complete data:  %.4f\n", kappa_complete))
cat(sprintf("Kappa from incomplete data via EM: %.4f\n\n", kappa_em))

cat("Explanation: The EM estimate should be close to the original, but not identical.\n")
cat("It represents the most likely kappa value given the information remaining\n")
cat("in the incomplete dataset under the MAR assumption.\n")
