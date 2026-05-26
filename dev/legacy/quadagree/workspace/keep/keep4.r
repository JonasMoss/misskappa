# --- The Final, Formal Delta Method Function ---
kappa_acov_formal_full <- function(x) {

  # --- 1. INITIAL SETUP ---
  n <- nrow(x)
  r <- ncol(x)
  M <- !is.na(x)

  mu_hat <- colMeans(x, na.rm = TRUE)
  sigma_hat <- cov(x, use = "pairwise.complete.obs")
  Y_hat <- sweep(x, 2, mu_hat, "-")

  # --- Needed Probabilities ---
  p1 <- colSums(M) / n
  p2 <- crossprod(M) / n

  # --- 2. THE GRADIENT (d(k_f, k_c) / d(mu, sigma_vec)) ---
  n_params <- r + r^2
  grad_mat <- matrix(NA, nrow = n_params, ncol = 2) # cols for F, C

  # Kappa point estimates
  tr_s <- sum(diag(sigma_hat))
  mu_bar <- mean(mu_hat)
  d_mu <- sum((mu_hat - mu_bar)^2)

  NF <- sum(sigma_hat) - tr_s - d_mu
  DF <- (r-1)*tr_s + (r-1)*d_mu
  kF <- NF/DF

  NC <- sum(sigma_hat) - tr_s
  DC <- (r-1)*tr_s + r*d_mu
  kC <- NC/DC

  # Gradient w.r.t mu
  dkF_dmu <- - (NF + DF) / (DF^2) * 2 * (mu_hat - mu_bar)
  dkC_dmu <- - (NC + DC) / (DC^2) * 2 * r/r * (mu_hat - mu_bar) # Note: d(r*d_mu)/dmu = 2r(mu-mubar)
  # Correction for dC's d_mu term
  dkC_dmu <- (-(DC * r * 2 * (mu_hat - mu_bar)) - (NC * r * 2 * (mu_hat - mu_bar))) / DC^2
  dkC_dmu <- -((NC+DC)*r*2*(mu_hat - mu_bar)) / DC^2

  grad_mat[1:r, 1] <- dkF_dmu
  grad_mat[1:r, 2] <- dkC_dmu

  # Gradient w.r.t sigma_kl
  for(k in 1:r) { for(l in 1:r) {
    idx <- r + (k-1)*r + l
    is_diag <- as.numeric(k == l)

    dkF_ds <- (DF*(1-is_diag) - NF*((r-1)*is_diag + (r-1))) / DF^2
    dkC_ds <- (DC*(1-is_diag) - NC*((r-1)*is_diag + r      )) / DC^2
    grad_mat[idx, 1] <- dkF_ds
    grad_mat[idx, 2] <- dkC_ds
  }}

  # --- 3. THE OMEGA_PW MATRIX ---
  Omega_pw <- matrix(NA, nrow = n_params, ncol = n_params)

  # Block 1: Omega_mumu
  Pi_S <- p2 / (p1 %o% p1)
  Omega_pw[1:r, 1:r] <- Pi_S * sigma_hat

  # Blocks 2 & 3: Omega_musigma
  C_mat <- matrix(NA, nrow = r, ncol = r^2)
  p3_mat <- array(NA, dim = c(r, r, r))
  for(j in 1:r) for(k in 1:r) for(l in 1:r) p3_mat[j,k,l] <- mean(M[,j]*M[,k]*M[,l])

  for(j in 1:r) { for(k in 1:r) { for(l in 1:r) {
    idx <- (k-1)*r + l
    C_mat[j, idx] <- mean(Y_hat[,j]*Y_hat[,k]*Y_hat[,l], na.rm = TRUE)
    Pi_C_val <- p3_mat[j,k,l] / (p1[j] * p2[k,l])
    Omega_pw[j, r + idx] <- Omega_pw[r + idx, j] <- Pi_C_val * C_mat[j, idx]
  }}}

  # Block 4: Omega_sigmasigma
  p4_mat <- array(NA, dim=c(r,r,r,r))
  # This is the slowest part
  for(i in 1:r) for(j in 1:r) for(k in 1:r) for(l in 1:r) p4_mat[i,j,k,l]<-mean(M[,i]*M[,j]*M[,k]*M[,l])

  for(i in 1:r) { for(j in 1:r) {
    for(k in 1:r) { for(l in 1:r) {
      row_idx <- r + (i-1)*r + j
      col_idx <- r + (k-1)*r + l
      mu_ijkl <- mean(Y_hat[,i]*Y_hat[,j]*Y_hat[,k]*Y_hat[,l], na.rm=TRUE)
      gamma_val <- mu_ijkl - sigma_hat[i,j] * sigma_hat[k,l]
      pi_gamma_val <- p4_mat[i,j,k,l] / (p2[i,j] * p2[k,l])
      Omega_pw[row_idx, col_idx] <- pi_gamma_val * gamma_val
    }}
  }}

  # --- 4. FINAL CALCULATION ---
  scaled_acov <- t(grad_mat) %*% Omega_pw %*% grad_mat

  n_eff <- sum(rowSums(M) > 0)

  list(
    kappa_fleiss = kF,
    kappa_cohen = kC,
    scaled_acov = scaled_acov,
    n_eff = n_eff
  )
}

# --- Simulation Study with the Formal Function ---
# NOTE: This will be VERY slow due to the O(R^4) loops for p4 and Omega.
# We will use a smaller n_reps.

set.seed(135)
n_sim <- 5000
n_reps <- 200 # Reduced for speed
r_sim <- 4
missing_prob <- 0.2

true_mu <- c(1, 1.5, 1.2, 2)
true_sigma <- matrix(0.7, r_sim, r_sim); diag(true_sigma) <- 1
L <- chol(true_sigma)

kappa_estimates <- matrix(NA, nrow = n_reps, ncol = 2)
scaled_acov_estimates <- array(NA, dim = c(2, 2, n_reps))
n_eff_vec <- numeric(n_reps)

cat("--- Starting Final Formal Monte Carlo ---\n")
pb <- txtProgressBar(min = 0, max = n_reps, style = 3)
for (i in 1:n_reps) {
  complete_data <- matrix(rnorm(n_sim * r_sim), n_sim, r_sim) %*% L
  complete_data <- sweep(complete_data, 2, true_mu, "+")
  x_missing <- complete_data
  x_missing[sample(1:(n_sim * r_sim), size = n_sim * r_sim * missing_prob)] <- NA

  res <- kappa_acov_formal_full(x_missing)

  kappa_estimates[i,] <- c(res$kappa_fleiss, res$kappa_cohen)
  scaled_acov_estimates[,,i] <- res$scaled_acov
  n_eff_vec[i] <- res$n_eff

  setTxtProgressBar(pb, i)
}
close(pb)

mean_n_eff <- mean(n_eff_vec)
empirical_scaled_acov <- cov(kappa_estimates) * mean_n_eff
mean_estimated_scaled_acov <- apply(scaled_acov_estimates, c(1, 2), mean)

cat("\n\n--- Full Formal Delta Method Results (N_sim =", n_sim, ", Missing Prob =", missing_prob, ") ---\n")
cat("\n'True' Scaled ACov (from empirical variance * n_eff):\n")
print(empirical_scaled_acov)

cat("\nMean Estimated Scaled ACov (from Formal Delta Method function):\n")
print(mean_estimated_scaled_acov)
