# --- The Verified Main Function ---
# Logic is taken directly from the successful identity test.
kappa_acov_verified <- function(x) {

  # --- 1. INITIAL SETUP ---
  n <- nrow(x)
  r <- ncol(x)
  M <- !is.na(x)

  mu_hat <- colMeans(x, na.rm = TRUE)
  Y_hat <- sweep(x, 2, mu_hat, "-")

  # Isolate the effective data used for Z*
  complete_rows_idx <- which(rowSums(M) > 0)
  Y_hat_eff <- Y_hat[complete_rows_idx, , drop = FALSE]
  M_eff <- M[complete_rows_idx, , drop = FALSE]
  n_eff <- nrow(Y_hat_eff)

  # Probabilities calculated from the effective sample
  p1 <- colSums(M_eff) / n_eff
  p2 <- crossprod(M_eff) / n_eff


  # --- 2. CLEVER Z METHOD (VERIFIED) ---
  Z_star_matrix_eff <- matrix(NA, nrow = n_eff, ncol = 3)
  mu_diff <- mu_hat - mean(mu_hat)
  inv_p1 <- 1 / p1
  inv_p2 <- 1 / p2

  for (i in 1:n_eff) {
    Mi <- M_eff[i,]
    Yi <- Y_hat_eff[i,]
    obs_raters <- which(Mi)

    z1 <- sum(mu_diff[Mi] * Yi[Mi] * inv_p1[Mi])
    z3 <- sum(Yi[Mi]^2 * inv_p1[Mi])

    z2_diag <- z3
    z2_offdiag <- 0
    if (length(obs_raters) > 1) {
      pairs <- combn(obs_raters, 2)
      for (k in 1:ncol(pairs)) {
        j1 <- pairs[1, k]; j2 <- pairs[2, k]
        z2_offdiag <- z2_offdiag + 2 * Yi[j1] * Yi[j2] * inv_p2[j1, j2]
      }
    }
    z2 <- z2_diag + z2_offdiag
    Z_star_matrix_eff[i, ] <- c(z1, z2, z3)
  }

  # 'Meat' of the sandwich estimator, using n_eff denominator
  Xi_star <- cov(Z_star_matrix_eff) * (n_eff - 1) / n_eff


  # --- 3. KAPPA ESTIMATES AND GRADIENTS ---
  sigma_hat <- cov(x, use = "pairwise.complete.obs")
  tr_sigma <- sum(diag(sigma_hat))
  d_mu <- sum((mu_hat - mean(mu_hat))^2)

  n_c <- sum(sigma_hat) - tr_sigma
  n_f <- n_c - d_mu
  d_c <- (r - 1) * tr_sigma + r * d_mu
  d_f <- (r-1) * tr_sigma + (r-1) * d_mu

  kappa_c <- n_c / d_c
  kappa_f <- n_f / d_f

  # Gradients
  cohen_grad <- c(2 * r * kappa_c, -1, 1 + (r - 1) * kappa_c) / d_c
  fleiss_grad <- c(2 * (1 + (r - 1) * kappa_f), -1, 1 + (r - 1) * kappa_f) / d_f
  grad_matrix <- cbind(fleiss_grad, cohen_grad)

  # This is the asymptotic covariance of sqrt(n_eff) * (kappa - E[kappa])
  scaled_acov <- t(grad_matrix) %*% Xi_star %*% grad_matrix

  list(
    kappa_fleiss = kappa_f,
    kappa_cohen = kappa_c,
    scaled_acov = scaled_acov,
    n_eff = n_eff
  )
}


# --- Simulation Study with the Verified Function ---

set.seed(789)
n_sim <- 5000
n_reps <- 200
r_sim <- 4
missing_prob <- 0.2 # The challenging case

# True population parameters for data generation
true_mu <- c(1, 1.5, 1.2, 2)
true_sigma <- matrix(0.7, r_sim, r_sim)
diag(true_sigma) <- 1
L <- chol(true_sigma)

# Storage for simulation results
kappa_estimates <- matrix(NA, nrow = n_reps, ncol = 2)
scaled_acov_estimates <- array(NA, dim = c(2, 2, n_reps))
n_eff_vec <- numeric(n_reps)

cat("--- Starting Final Monte Carlo Verification ---\n")
pb <- txtProgressBar(min = 0, max = n_reps, style = 3)
for (i in 1:n_reps) {
  # Generate new complete data
  complete_data <- matrix(rnorm(n_sim * r_sim), n_sim, r_sim) %*% L
  complete_data <- sweep(complete_data, 2, true_mu, "+")

  # Introduce missingness
  x_missing <- complete_data
  x_missing[sample(1:(n_sim * r_sim), size = n_sim * r_sim * missing_prob)] <- NA

  # Run the VERIFIED function
  res <- kappa_acov_verified(x_missing)

  # Store results
  kappa_estimates[i,] <- c(res$kappa_fleiss, res$kappa_cohen)
  scaled_acov_estimates[,,i] <- res$scaled_acov
  n_eff_vec[i] <- res$n_eff

  setTxtProgressBar(pb, i)
}
close(pb)

# --- Analyze and Compare Results ---

mean_n_eff <- mean(n_eff_vec)
empirical_scaled_acov <- cov(kappa_estimates*sqrt(n_eff_vec))
mean_estimated_scaled_acov <- apply(scaled_acov_estimates, c(1, 2), mean)

cat("\n\n--- Final Monte Carlo Results (N_sim =", n_sim, ", Missing Prob =", missing_prob, ") ---\n")
cat("--- Using Fully Verified Code ---\n")

cat("\n'True' Scaled ACov (from empirical variance * n_eff):\n")
print(empirical_scaled_acov)

cat("\nMean Estimated Scaled ACov (from Verified 'Clever Z' function):\n")
print(mean_estimated_scaled_acov)
