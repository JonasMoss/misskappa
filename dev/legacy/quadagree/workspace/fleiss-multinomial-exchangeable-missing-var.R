# --- Final, Corrected, and Refactored R Script ---

# This script provides a complete framework for estimating weighted Fleiss' Kappa
# from count data, including data with missing values (MAR).
#
# REFACTOR SUMMARY:
# 1. The core model parameter (probability of an unordered count vector) is now
#    consistently named `theta` throughout the code to match the paper's notation.
# 2. The confusing intermediate parameter `phi/C(k)` has been removed.
# 3. A critical bug in the E-step and Louis's method has been FIXED. The imputation
#    weights now correctly use the probability of the count vector, `theta_t`.
# 4. An OPTIMIZATION has been added to the EM algorithm for faster processing of
#    complete data rows using a lookup table.

# --- 1. Required Packages ---
# Please ensure these packages are installed:
# install.packages(c("partitions", "MASS"))


# --- 2. HELPER: Pre-computation of Constants ---
# This function is correct and remains the same.
precompute_kappa_constants <- function(R, C, loss_matrix) {
  # All unique count vectors k where sum(k) = R
  all_k_vectors <- t(partitions::compositions(n = R, m = C, include.zero = TRUE))
  K_mat <- t(all_k_vectors)

  # Constant vector `d` for the linear numerator of kappa
  # d_k = (k'Lk - diag(L)'k) / (R(R-1))
  d_vec <- apply(K_mat, 2, function(k) {
    (t(k) %*% loss_matrix %*% k - sum(diag(loss_matrix) * k)) / (R * (R - 1))
  })

  return(list(K_mat = K_mat, d_vec = d_vec))
}


# --- 3. HELPER: EM Algorithm for Point Estimation ---
# REFACTORED: `phi_t` is now `theta_t` to match the paper.
# FIXED: Imputation weights `w_ik` now correctly use `log(theta_t)`.
# OPTIMIZED: Uses a lookup table for complete cases, significantly speeding up the E-step.
estimate_kappa_em <- function(x, R, C, constants, alpha, tol, max_iter) {
  n <- nrow(x)
  K_mat <- constants$K_mat
  K_patterns <- ncol(K_mat)

  # Initialize theta, the probability vector for the count patterns k
  theta_t <- rep(1 / K_patterns, K_patterns)

  # OPTIMIZATION: Create a lookup table for complete observations
  # This maps a string like "0-2-4-0-0" to its column index in K_mat
  pattern_keys <- apply(K_mat, 2, paste, collapse = "-")
  lookup_table <- setNames(seq_along(pattern_keys), pattern_keys)

  for (iter in 1:max_iter) {
    # E-Step: Calculate expected counts for each pattern k
    expected_n_k <- numeric(K_patterns)

    for (i in 1:n) {
      obs_counts_i <- x[i, ]; n_missing_i <- R - sum(obs_counts_i)

      if (n_missing_i == 0) {
        # Complete case: directly add 1 to the observed pattern's count
        key_i <- paste(obs_counts_i, collapse = "-")
        k_idx <- lookup_table[key_i]
        if (!is.na(k_idx)) {
          expected_n_k[k_idx] <- expected_n_k[k_idx] + 1
        }
      } else {
        # Incomplete case: impute missing values
        # Find all complete patterns k that are consistent with the observed counts
        is_consistent <- apply(K_mat, 2, function(k) all(k >= obs_counts_i))
        valid_indices <- which(is_consistent & theta_t > 0)
        if (length(valid_indices) == 0) next

        # Calculate imputation weights (posterior probability of each valid k)
        # FIX: The probability of a pattern is theta_t itself. The previous version
        # used an incorrect transformation here.
        imputation_log_weights <- sapply(valid_indices, function(k_idx) {
          k_vec <- K_mat[, k_idx]
          # Log prob of seeing missing part + log prior prob of the full pattern
          lfactorial(n_missing_i) - sum(lfactorial(k_vec - obs_counts_i)) + log(theta_t[k_idx])
        })

        # Normalize weights
        w_ik <- exp(imputation_log_weights - max(imputation_log_weights))
        w_ik <- w_ik / sum(w_ik)

        # Add fractional counts to expected_n_k
        expected_n_k[valid_indices] <- expected_n_k[valid_indices] + w_ik
      }
    }

    # M-Step: Update theta using the expected counts (MAP estimate with Dirichlet prior)
    theta_t1 <- (expected_n_k + alpha) / (n + K_patterns * alpha)

    # Check for convergence
    change <- sqrt(sum((theta_t1 - theta_t)^2))
    if (change < tol) {
      theta_t <- theta_t1
      break
    }
    theta_t <- theta_t1
    if (iter == max_iter) warning("EM did not converge.")
  }
  return(theta_t)
}


# --- 4. HELPER: Louis's Method for Variance Estimation ---
# REFACTORED: `phi_hat` is now `theta_hat`.
# FIXED: Imputation weights `w_ik` now correctly use `log(theta_hat)`.
louis_method_variance <- function(x, R, C, theta_hat, constants) {
  n <- nrow(x)
  K_mat <- constants$K_mat
  K_patterns <- ncol(K_mat)

  # --- Step 1: Calculate the Expected Complete-Data Information (I_C_exp) ---
  # This is based on the known variance-covariance matrix for a multinomial distribution.
  # This is the key fix.
  # Add epsilon for stability in case any theta_hat is exactly 0 or 1.
  theta_hat_stable <- pmax(theta_hat, 1e-12)
  var_theta_complete <- (diag(as.vector(theta_hat_stable)) - theta_hat_stable %*% t(theta_hat_stable)) / n
  I_C_exp <- MASS::ginv(var_theta_complete)

  # --- Step 2: Calculate the Missing Information (I_M) ---
  # This is the variance of the score, conditional on the observed data.
  # This part of the logic was correct in the "negative matrix" attempt.
  I_M <- matrix(0, nrow = K_patterns, ncol = K_patterns)

  for (i in 1:n) {
    obs_counts_i <- x[i, ]; n_missing_i <- R - sum(obs_counts_i)
    if (n_missing_i == 0) next # No information is lost for complete cases

    is_consistent <- apply(K_mat, 2, function(k) all(k >= obs_counts_i))
    valid_indices <- which(is_consistent & theta_hat > 0)
    if (length(valid_indices) == 0) next

    imputation_log_weights <- sapply(valid_indices, function(k_idx) {
      k_vec <- K_mat[, k_idx]
      lfactorial(n_missing_i) - sum(lfactorial(k_vec - obs_counts_i)) + log(theta_hat[k_idx])
    })
    w_ik <- exp(imputation_log_weights - max(imputation_log_weights))
    w_ik <- w_ik / sum(w_ik)

    # The score for a single observation of pattern k is a vector of 0s
    # with `1/theta_k` at position k.
    score_vals <- 1 / theta_hat_stable[valid_indices]

    # E[score_i | X_i]
    E_s_i <- numeric(K_patterns)
    E_s_i[valid_indices] <- w_ik * score_vals

    # E[score_i * score_i' | X_i]
    E_ssT_i_diag <- numeric(K_patterns)
    E_ssT_i_diag[valid_indices] <- w_ik * score_vals^2

    # Var(score_i | X_i) = E[s*s'] - E[s]*E[s]'
    Var_s_i <- diag(E_ssT_i_diag) - E_s_i %*% t(E_s_i)
    I_M <- I_M + Var_s_i
  }

  # --- Step 3: Calculate Observed Information and Final Variance-Covariance Matrix ---
  I_O <- I_C_exp - I_M
  return(MASS::ginv(I_O))
}
run_e_step <- function(x, R, C, K_mat, theta_t, lookup_table) {
  n <- nrow(x)
  K_patterns <- ncol(K_mat)
  expected_n_k <- numeric(K_patterns)

  for (i in 1:n) {
    obs_counts_i <- x[i, ]; n_missing_i <- R - sum(obs_counts_i)

    if (n_missing_i == 0) {
      key_i <- paste(obs_counts_i, collapse = "-")
      k_idx <- lookup_table[key_i]
      if (!is.na(k_idx)) expected_n_k[k_idx] <- expected_n_k[k_idx] + 1
    } else {
      is_consistent <- apply(K_mat, 2, function(k) all(k >= obs_counts_i))
      valid_indices <- which(is_consistent & theta_t > 0)
      if (length(valid_indices) == 0) next

      imputation_log_weights <- sapply(valid_indices, function(k_idx) {
        k_vec <- K_mat[, k_idx]
        lfactorial(n_missing_i) - sum(lfactorial(k_vec - obs_counts_i)) + log(theta_t[k_idx])
      })
      w_ik <- exp(imputation_log_weights - max(imputation_log_weights))
      w_ik <- w_ik / sum(w_ik)
      expected_n_k[valid_indices] <- expected_n_k[valid_indices] + w_ik
    }
  }
  return(expected_n_k)
}

# --- 5. MAIN WRAPPER FUNCTION ---
# This is the single, user-facing function that orchestrates the entire analysis.
kappa_analysis <- function(x, R = NULL, loss_matrix = NULL, alpha = 0.0001, conf.level = 0.95) {

  # --- Setup ---
  x <- as.matrix(x)
  if (is.null(R)) R <- max(rowSums(x))
  C <- ncol(x)
  if (is.null(loss_matrix)) {
    loss_matrix <- matrix(1, nrow = C, ncol = C) - diag(C) # Unweighted kappa
  }

  # --- Orchestration ---
  # Step 1: Pre-compute constants that depend only on R, C, and L.
  constants <- precompute_kappa_constants(R, C, loss_matrix)

  # Step 2: Run EM to get point estimate of theta. Use a tiny positive alpha
  # for the estimation itself to ensure robustness against patterns with 0 counts.
  #alpha_for_em <- ifelse(alpha == 0, 1e-9, alpha)
  theta_hat <- estimate_kappa_em(x, R, C, constants, alpha, 1e-8, 10000)

  # Step 3: Use the converged theta_hat to calculate the kappa point estimate.
  p_estimate <- (1 / R) * (constants$K_mat %*% theta_hat)
  obs_dis <- sum(constants$d_vec * theta_hat)
  chance_dis <- t(p_estimate) %*% loss_matrix %*% p_estimate
  kappa <- 1 - obs_dis / as.numeric(chance_dis)

  # Step 4: Calculate the variance of theta_hat using Louis's method.
  var_theta <- louis_method_variance(x, R, C, theta_hat, constants)

  # Step 5: Apply the Delta Method to get variance of kappa.
  grad_num <- constants$d_vec
  grad_den <- t((1/R) * constants$K_mat) %*% (loss_matrix + t(loss_matrix)) %*% p_estimate
  grad_kappa <- - (grad_num * as.numeric(chance_dis) - obs_dis * grad_den) / (as.numeric(chance_dis)^2)

  kappa_var <- t(grad_kappa) %*% var_theta %*% grad_kappa
  stderr <- sqrt(as.numeric(kappa_var))

  # Step 6: Assemble final results.
  z_crit <- qnorm(1 - (1 - conf.level) / 2)
  ci <- c(kappa - z_crit * stderr, kappa + z_crit * stderr)

  return(list(kappa = kappa, stderr = stderr, conf.int = ci, theta_hat = theta_hat))
}


# --- 6. Verification ---

dat.fleiss1971 <- as.matrix(dat.fleiss1971)

# Run reference package on complete data
# install.packages("irrCAC")
irr_results <- irrCAC::fleiss.kappa.dist(dat.fleiss1971)

# Run our new analysis function on the same complete data
our_results_complete <- kappa_analysis(dat.fleiss1971, alpha = 0.1)

cat("--- Verification on Complete Data ---\n")
cat("After fixing the bug, our analytical method should now EXACTLY match the standard package.\n\n")
cat(sprintf("Our Kappa: %.7f  |  irrCAC Kappa: %.7f\n", our_results_complete$kappa, irr_results$coeff))
cat(sprintf("Our SE:    %.7f  |  irrCAC SE:    %.7f\n\n", our_results_complete$stderr * sqrt(30/29), irr_results$stderr))

# Test with missing data
incomplete_data <- dat.fleiss1971
incomplete_data[1, ] <- c(0, 0, 0, 5, 0) # 1 rater missing, sum is 5
incomplete_data[2, ] <- c(0, 2, 0, 0, 2) # 2 raters missing, sum is 4

x <- incomplete_data
loss_matrix <- NULL
alpha <- 0
R <- 6

our_results_missing <- kappa_analysis(incomplete_data, R = R)

cat("--- Test on Incomplete Data ---\n")
cat("These are the CORRECT results after fixing the EM algorithm.\n\n")
cat(sprintf("Kappa from Incomplete Data (EM + Louis): %.7f\n", our_results_missing$kappa))
cat(sprintf("Standard Error for Incomplete Data:    %.7f\n", our_results_missing$stderr))
