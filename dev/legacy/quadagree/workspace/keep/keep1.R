# Helper function to get initial estimates and probabilities
get_initials <- function(x) {
  n <- nrow(x)
  r <- ncol(x)

  # Indicator matrix for non-missing values
  M <- !is.na(x)

  # --- Pairwise estimates for mu, sigma, and higher moments ---
  mu_hat <- colMeans(x, na.rm = TRUE)
  sigma_hat <- cov(x, use = "pairwise.complete.obs")

  # For C and Gamma, we need centered data
  Y_hat <- sweep(x, 2, mu_hat, "-")

  # --- Observation probabilities ---
  p1 <- colSums(M) / n # Vector p_j
  p2 <- crossprod(M) / n # Matrix p_jk

  # We need p_jkl and p_jklm for the formal method. These are expensive.
  # Let's compute them here for use in the formal function.
  p3 <- array(NA, dim = c(r, r, r))
  p4 <- array(NA, dim = c(r, r, r, r))

  for (j in 1:r) {
    for (k in j:r) {
      for (l in k:r) {
        p_jkl <- mean(M[,j] * M[,k] * M[,l])
        p3[j,k,l] <- p3[j,l,k] <- p3[k,j,l] <- p3[k,l,j] <- p3[l,j,k] <- p3[l,k,j] <- p_jkl
        for (m in l:r) {
          p_jklm <- mean(M[,j] * M[,k] * M[,l] * M[,m])
          # This is just a sketch; a full implementation would handle all permutations
          p4[j,k,l,m] <- p_jklm
        }
      }
    }
  }

  list(
    n = n, r = r, x = x, Y_hat = Y_hat, M = M,
    mu_hat = mu_hat, sigma_hat = sigma_hat,
    p1 = p1, p2 = p2, p3 = p3, p4 = p4
  )
}

kappa_acov_formal <- function(x) {

  # --- Get initial quantities ---
  d <- get_initials(x)
  n <- d$n; r <- d$r; Y_hat <- d$Y_hat

  # --- Build the full Asymptotic Covariance Matrix Omega_pw ---
  # Total parameters: r means + r*r covariances
  n_params <- r + r^2
  Omega_pw <- matrix(NA, nrow = n_params, ncol = n_params)

  # Define indices for mu and sigma in the large vector
  mu_idx <- 1:r
  sig_idx <- (r + 1):n_params

  # Block 1: Omega_mumu (Covariance of means)
  Pi_Sigma <- d$p2 / (d$p1 %o% d$p1)
  Omega_mumu <- Pi_Sigma * d$sigma_hat
  Omega_pw[mu_idx, mu_idx] <- Omega_mumu

  # Blocks 2 & 3: Omega_musigma (Cross-covariance)
  # This requires 3rd order moments
  C_mat <- matrix(NA, nrow = r, ncol = r^2)
  Pi_C <- matrix(NA, nrow = r, ncol = r^2)

  for (j in 1:r) {
    for (k in 1:r) {
      for (l in 1:r) {
        # Using listwise deletion for each specific moment
        c_j_kl <- mean(Y_hat[,j] * Y_hat[,k] * Y_hat[,l], na.rm = TRUE)
        # Get column index in the vectorized r^2 matrix
        col_idx <- (k - 1) * r + l
        C_mat[j, col_idx] <- c_j_kl
        Pi_C[j, col_idx] <- d$p3[j,k,l] / (d$p1[j] * d$p2[k,l])
      }
    }
  }
  Omega_musigma <- Pi_C * C_mat
  Omega_pw[mu_idx, sig_idx] <- Omega_musigma
  Omega_pw[sig_idx, mu_idx] <- t(Omega_musigma)

  # Block 4: Omega_sigmasigma (Gamma matrix)
  # This requires 4th order moments
  Gamma_mat <- matrix(NA, nrow = r^2, ncol = r^2)
  Pi_Gamma <- matrix(NA, nrow = r^2, ncol = r^2)

  for (i in 1:r) { for (j in 1:r) {
    for (k in 1:r) { for (l in 1:r) {
      row_idx <- (i - 1) * r + j
      col_idx <- (k - 1) * r + l

      mu_ijkl <- mean(Y_hat[,i]*Y_hat[,j]*Y_hat[,k]*Y_hat[,l], na.rm = TRUE)
      Gamma_mat[row_idx, col_idx] <- mu_ijkl - d$sigma_hat[i,j] * d$sigma_hat[k,l]

      # Getting the right p4 value is tricky, this is a simplification
      p_ijkl <- d$p4[min(i,j,k,l), sort(c(i,j,k,l))[2], sort(c(i,j,k,l))[3], max(i,j,k,l)]
      Pi_Gamma[row_idx, col_idx] <- p_ijkl / (d$p2[i,j] * d$p2[k,l])
    }}
  }}
  Omega_sigmasigma <- Pi_Gamma * Gamma_mat
  Omega_pw[sig_idx, sig_idx] <- Omega_sigmasigma

  # --- Delta Method Gradient ---
  # The gradient of (kappa_f, kappa_c) w.r.t (mu, sigma_vec)
  # This is a complex derivation, we'll represent it as a placeholder
  # In a real scenario, this would be a large function.
  # For now, we know the final result only depends on the Z-vector gradients.
  # Let's just return Omega_pw for comparison.
  # A full implementation would compute the full gradient and the final 2x2 matrix.

  return(Omega_pw) # For debugging/comparison
}

kappa_acov_clever_Z <- function(x) {

  # --- Get initial quantities (don't need p3, p4) ---
  d <- get_initials(x)
  n <- d$n; r <- d$r; Y_hat <- d$Y_hat; M <- d$M
  mu_hat <- d$mu_hat; p1 <- d$p1; p2 <- d$p2

  # --- Construct the n x 3 Clever Z matrix ---
  Z_star_matrix <- matrix(NA, nrow = n, ncol = 3)

  mu_bar <- mean(mu_hat)
  mu_diff <- mu_hat - mu_bar

  # Pre-calculate inverse probabilities to avoid division in loop
  inv_p1 <- 1 / p1
  inv_p2 <- 1 / p2

  for (i in 1:n) {
    Mi <- M[i,]
    if (sum(Mi) == 0) next

    Yi <- Y_hat[i,]

    # Z_i1*: Mean cross-product term
    z1 <- sum(mu_diff[Mi] * Yi[Mi] * inv_p1[Mi])

    # Z_i3*: Trace term
    z3 <- sum(Yi[Mi]^2 * inv_p1[Mi])

    # Z_i2*: Total sum of squares term
    z2_diag <- z3 # Diagonal part is the same as the trace

    z2_offdiag <- 0
    obs_raters <- which(Mi)
    if (length(obs_raters) > 1) {
      # Loop over unique pairs of observed raters
      pairs <- combn(obs_raters, 2)
      for (k in 1:ncol(pairs)) {
        j1 <- pairs[1, k]
        j2 <- pairs[2, k]
        # Factor of 2 for symmetry
        z2_offdiag <- z2_offdiag + 2 * Yi[j1] * Yi[j2] * inv_p2[j1, j2]
      }
    }
    z2 <- z2_diag + z2_offdiag

    Z_star_matrix[i, ] <- c(z1, z2, z3)
  }

  # Remove rows corresponding to subjects with no ratings
  Z_star_matrix <- Z_star_matrix[complete.cases(Z_star_matrix), ]

  # The covariance matrix of the Z* vectors
  Xi_star <- cov(Z_star_matrix)

  # --- Compute Final Kappa ACov ---
  # First, get the kappa estimates
  tr_sigma <- sum(diag(d$sigma_hat))
  d_mu <- sum((mu_hat - mu_bar)^2)

  n_c <- sum(d$sigma_hat) - tr_sigma
  n_f <- n_c - d_mu
  d_c <- (r - 1) * tr_sigma + r * d_mu
  d_f <- (r-1) * tr_sigma + (r-1) * d_mu

  kappa_c <- n_c / d_c
  kappa_f <- n_f / d_f

  # Gradients w.r.t W = [W1, W2, W3]
  # Where W1=(mu-mu_bar)'mu, W2=1'S1, W3=tr(S)
  # This is the gradient of kappa w.r.t the EXPECTATIONS of Z*
  # xi_F and xi_C from your original paper are exactly these gradients.

  # Wait, the original gradients were w.r.t. Z = [Z1, Z2, Z3] not E[Z]
  # Z1 = (mu-mu_bar)'Y, Z2 = (1'Y)^2, Z3 = Y'Y
  # E[Z1]=0, E[Z2]=1'S1, E[Z3]=tr(S)
  # This is correct. The delta method is applied to the function of Z.

  cohen_grad <- c(2 * r * kappa_c, -1, 1 + (r - 1) * kappa_c) / d_c
  fleiss_grad <- c(2 * (1 + (r - 1) * kappa_f), -1, 1 + (r - 1) * kappa_f) / d_f

  grad_matrix <- cbind(fleiss_grad, cohen_grad) # Note order F, C

  # The final 2x2 Asymptotic Covariance Matrix
  # The n denominator comes from the sqrt(n) in the CLT
  acov <- t(grad_matrix) %*% Xi_star %*% grad_matrix

  return(list(
    kappa_fleiss = kappa_f,
    kappa_cohen = kappa_c,
    acov = acov, # The final variance is acov/n
    Xi_star = Xi_star # For debugging
  ))
}

set.seed(123)
n_sim <- 1000
r_sim <- 4

# Generate complete data from a multivariate normal distribution
true_mu <- c(1, 1.5, 1.2, 2)
true_sigma <- matrix(0.7, r_sim, r_sim)
diag(true_sigma) <- 1
L <- chol(true_sigma)
complete_data <- matrix(rnorm(n_sim * r_sim), n_sim, r_sim) %*% L
complete_data <- sweep(complete_data, 2, true_mu, "+")

# Introduce missingness (MCAR)
missing_prob <- 0.2
x_missing <- complete_data
x_missing[sample(1:(n_sim * r_sim), size = n_sim * r_sim * missing_prob)] <- NA

# --- Run the Clever Z function ---
clever_results <- kappa_acov_clever_Z(x_missing)

cat("--- Clever Z Method Results ---\n")
print(clever_results$acov)

# --- The problem with the Formal Method ---
# As written, the formal method is too complex to verify quickly.
# The crucial insight is that the Z-method *is* the delta method.
# The covariance matrix Xi_star *is* the result of G' * Omega_pw * G
# where G is the gradient of Z w.r.t theta = (mu, sigma).
#
# So, instead of comparing the final 2x2 matrices, let's verify that
# Xi_star is indeed capturing the full covariance structure. We can
# compare elements of Xi_star to manually computed elements using the
# formal theory.

# Let's check Cov(Z2*, Z3*)
# Cov(Z2*, Z3*) = Cov(IPW(1'S1), IPW(tr(S)))
# From the delta method, this depends on all 3rd and 4th order moments.
# The 'kappa_acov_formal' function is already a sketch of this.
# A full verification would require completing it and debugging the indices,
# which is a significant task.

# The core of the verification lies in the fact that:
# E[Z_i*] = [0, 1'S1, tr(S)]
# And by the properties of IPW estimators, the covariance matrix
# of the Z* vectors will converge to the correct asymptotic
# covariance matrix required by the delta method.

# --- Simulation to check correctness ---
# We can run a Monte Carlo simulation.
# 1. Repeatedly generate data with missingness.
# 2. In each replication, calculate kappa_f and kappa_c.
# 3. The empirical covariance of these estimates across replications should match
#    the average of the estimated ACov matrices from our function.

n_reps <- 2000
kappa_estimates <- matrix(NA, nrow = n_reps, ncol = 2)
acov_estimates <- array(NA, dim = c(2, 2, n_reps))

cat("\n--- Starting Monte Carlo Verification ---\n")
for (i in 1:n_reps) {
  # Generate new data for each rep
  sim_data <- matrix(rnorm(n_sim * r_sim), n_sim, r_sim) %*% L
  sim_data <- sweep(sim_data, 2, true_mu, "+")
  sim_data[sample(1:(n_sim * r_sim), size = n_sim * r_sim * missing_prob)] <- NA

  n_eff <- sum(apply(sim_data, 1, \(x) sum(is.na(x)) != R))
  res <- kappa_acov_clever_Z(sim_data)
  kappa_estimates[i,] <- c(res$kappa_fleiss, res$kappa_cohen)
  acov_estimates[,,i] <- res$acov
  if(i %% 20 == 0) cat("Rep:", i, "/", n_reps, "\n")
}

# Empirical variance from the simulation
empirical_acov <- cov(kappa_estimates)

# Average estimated variance from our function
mean_estimated_acov <- apply(acov_estimates, c(1, 2), mean)

cat("\n--- Monte Carlo Results ---\n")
cat("Empirical (True) ACov from Simulation:\n")
print(empirical_acov*n_sim)

cat("\nAverage Estimated ACov from 'Clever Z' function:\n")
print(mean_estimated_acov)
