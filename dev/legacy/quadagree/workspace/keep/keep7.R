# --- The "Aggregate Statistics" Method ---
kappa_acov_aggregate <- function(x) {

  # --- 1. INITIAL SETUP (Same as formal method) ---
  n <- nrow(x); r <- ncol(x); M <- !is.na(x)
  mu_hat <- colMeans(x, na.rm=TRUE)
  sigma_hat <- cov(x, use="pairwise.complete.obs")
  Y_hat <- sweep(x, 2, mu_hat, "-")
  p1 <- colMeans(M); p2 <- crossprod(M) / n

  # Pre-calculate expensive moments and probabilities
  mu3 <- array(0, dim=c(r,r,r)); mu4 <- array(0, dim=c(r,r,r,r))
  p3 <- array(0, dim=c(r,r,r)); p4 <- array(0, dim=c(r,r,r,r))
  # This is the bottleneck, but necessary
  for(i in 1:r) for(j in 1:r) for(k in 1:r) {
    mu3[i,j,k] <- mean(Y_hat[,i]*Y_hat[,j]*Y_hat[,k], na.rm=TRUE)
    p3[i,j,k] <- mean(M[,i]*M[,j]*M[,k])
    for(l in 1:r) {
      mu4[i,j,k,l] <- mean(Y_hat[,i]*Y_hat[,j]*Y_hat[,k]*Y_hat[,l], na.rm=TRUE)
      p4[i,j,k,l] <- mean(M[,i]*M[,j]*M[,k]*M[,l])
    }
  }

  # --- 2. 3x3 Asymptotic Covariance Matrix of the Aggregates ---
  # Let theta_hat = [hat(1'S1), hat(tr(S)), hat(D_mu)]
  Cov_theta <- matrix(NA, 3, 3)

  # Helper to compute Cov(sum(v1_ij*s_ij), sum(v2_kl*s_kl))
  # where v1, v2 are selection vectors/matrices (e.g., 1 for total, c for trace)
  get_cov_of_sums <- function(v1, v2, type) {
    total_cov <- 0
    if (type == "sigma-sigma") {
      for(i in 1:r) for(j in 1:r) for(k in 1:r) for(l in 1:r) {
        gamma <- mu4[i,j,k,l] - sigma_hat[i,j]*sigma_hat[k,l]
        pi_gamma <- p4[i,j,k,l] / (p2[i,j]*p2[k,l])
        total_cov <- total_cov + v1[i,j] * v2[k,l] * pi_gamma * gamma
      }
    } else if (type == "mu-mu") {
      for(i in 1:r) for(j in 1:r) {
        omega_ij <- (p2[i,j]/(p1[i]*p1[j])) * sigma_hat[i,j]
        total_cov <- total_cov + v1[i] * v2[j] * omega_ij
      }
    } else if (type == "mu-sigma") { # v1 for mu, v2 for sigma
      for(i in 1:r) for(j in 1:r) for(k in 1:r) {
        omega_ijk <- (p3[i,j,k]/(p1[i]*p2[j,k])) * mu3[i,j,k]
        total_cov <- total_cov + v1[i] * v2[j,k] * omega_ijk
      }
    }
    return(total_cov)
  }

  # Define selection vectors for the sums
  v_total <- matrix(1, r, r)
  v_trace <- diag(1, r, r)
  v_dmu <- 2 * (mu_hat - mean(mu_hat)) # This is the gradient of D_mu wrt mu

  # Var(hat(1'S1)) = 1' * Omega_ss * 1
  Cov_theta[1,1] <- get_cov_of_sums(v_total, v_total, "sigma-sigma")
  # Var(hat(tr(S))) = c' * Omega_ss * c
  Cov_theta[2,2] <- get_cov_of_sums(v_trace, v_trace, "sigma-sigma")
  # Var(hat(D_mu)) = (dD/dmu)' * Omega_mm * (dD/dmu)
  Cov_theta[3,3] <- get_cov_of_sums(v_dmu, v_dmu, "mu-mu")

  # Cov(hat(1'S1), hat(tr(S)))
  Cov_theta[1,2] <- Cov_theta[2,1] <- get_cov_of_sums(v_total, v_trace, "sigma-sigma")
  # Cov(hat(1'S1), hat(D_mu))
  Cov_theta[1,3] <- Cov_theta[3,1] <- get_cov_of_sums(v_dmu, v_total, "mu-sigma")
  # Cov(hat(tr(S)), hat(D_mu))
  Cov_theta[2,3] <- Cov_theta[3,2] <- get_cov_of_sums(v_dmu, v_trace, "mu-sigma")


  # --- 3. GRADIENTS of kappas w.r.t the Aggregates ---
  # theta = [theta1, theta2, theta3] = [1'S1, tr(S), D_mu]
  # kF = (theta1 - theta2 - theta3) / ((r-1)*theta2 + (r-1)*theta3)
  # kC = (theta1 - theta2) / ((r-1)*theta2 + r*theta3)

  tr_s <- sum(diag(sigma_hat));
  d_mu <- sum((mu_hat - mean(mu_hat))^2)
  theta1_h <- sum(sigma_hat);
  theta2_h <- tr_s; theta3_h <- d_mu

  # Fleiss
  NF <- theta1_h - theta2_h - theta3_h
  DF <- (r-1)*theta2_h + (r-1)*theta3_h
  grad_F <- c( 1/DF,                                  # d/d_theta1
               (-DF - NF*(r-1)) / DF^2,              # d/d_theta2
               (-DF - NF*(r-1)) / DF^2 )             # d/d_theta3

  # Conger
  NC <- theta1_h - theta2_h
  DC <- (r-1)*theta2_h + r*theta3_h
  grad_C <- c( 1/DC,                                  # d/d_theta1
               (-DC - NC*(r-1)) / DC^2,              # d/d_theta2
               (-NC*r) / DC^2 )                      # d/d_theta3

  grad_matrix <- cbind(grad_F, grad_C)

  # --- 4. FINAL DELTA METHOD ---
  scaled_acov <- t(grad_matrix) %*% Cov_theta %*% grad_matrix

  list(scaled_acov = scaled_acov)
}


# --- Verification Run ---
# We run this new function against our trusted formal one.
set.seed(135)
n_sim <- 5000; r_sim <- 4; missing_prob <- 0.2
true_mu <- c(1, 1.5, 1.2, 2); true_sigma <- matrix(0.7, r_sim, r_sim); diag(true_sigma) <- 1
L <- chol(true_sigma)
complete_data <- matrix(rnorm(n_sim * r_sim), n_sim, r_sim) %*% L
complete_data <- sweep(complete_data, 2, true_mu, "+")
x_missing <- complete_data
x_missing[sample(1:(n_sim * r_sim), size = n_sim * r_sim * missing_prob)] <- NA

cat("\n--- Aggregate Statistics Method ---\n")
print(kappa_acov_aggregate(x_missing))

cat("\n--- Correct Formal Method (for comparison) ---\n")
print(kappa_acov_formal_CORRECTED(x_missing)$scaled_acov)
