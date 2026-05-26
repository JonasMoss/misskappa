# --- The Fixed Clever Z with Cross-Fitting ---
kappa_acov_Z_fixed <- function(x) {

  n <- nrow(x)
  r <- ncol(x)

  # 1. Split the data
  fold_ids <- sample(1:2, size = n, replace = TRUE)

  # --- Helper function to calculate Z* for one fold using params from the other ---
  calculate_Z_for_fold <- function(data_fold, M_fold, params_other_fold) {
    n_fold <- nrow(data_fold)
    Z_star_fold <- matrix(NA, nrow = n_fold, ncol = 3)

    # Unpack params from the other fold
    mu_hat <- params_other_fold$mu
    mu_bar <- mean(mu_hat)
    mu_diff <- mu_hat - mu_bar
    inv_p1 <- 1 / params_other_fold$p1
    inv_p2 <- 1 / params_other_fold$p2

    Y_hat_fold <- sweep(data_fold, 2, mu_hat, "-")

    for (i in 1:n_fold) {
      Mi <- M_fold[i,]
      if(sum(Mi) == 0) next
      Yi <- Y_hat_fold[i,]
      obs_raters <- which(Mi)

      z1 <- sum(mu_diff[Mi] * Yi[Mi] * inv_p1[Mi])
      z3 <- sum(Yi[Mi]^2 * inv_p1[Mi])

      z2_diag <- z3
      z2_offdiag <- 0
      if (length(obs_raters) > 1) {
        pairs <- combn(obs_raters, 2)
        for(k in 1:ncol(pairs)) {
          j1 <- pairs[1,k]; j2 <- pairs[2,k]
          z2_offdiag <- z2_offdiag + 2 * Yi[j1] * Yi[j2] * inv_p2[j1,j2]
        }
      }
      z2 <- z2_diag + z2_offdiag
      Z_star_fold[i, ] <- c(z1, z2, z3)
    }
    return(Z_star_fold)
  }

  # --- Process Fold 1 ---
  x1 <- x[fold_ids == 1, ]; M1 <- !is.na(x1)
  params1 <- list(mu = colMeans(x1, na.rm=TRUE), p1 = colMeans(M1), p2 = crossprod(M1)/nrow(M1))

  # --- Process Fold 2 ---
  x2 <- x[fold_ids == 2, ]; M2 <- !is.na(x2)
  params2 <- list(mu = colMeans(x2, na.rm=TRUE), p1 = colMeans(M2), p2 = crossprod(M2)/nrow(M2))

  # 2. & 3. Cross-fitting
  Z_star_for_fold1 <- calculate_Z_for_fold(x1, M1, params2)
  Z_star_for_fold2 <- calculate_Z_for_fold(x2, M2, params1)

  # 4. Combine
  Z_star_full <- rbind(Z_star_for_fold1, Z_star_for_fold2)

  # Filter out any rows that had no observations (should be few)
  complete_rows_idx <- which(!is.na(Z_star_full[,1]))
  Z_star_eff <- Z_star_full[complete_rows_idx, ]
  n_eff <- nrow(Z_star_eff)

  # --- Calculate final results using the full dataset's parameters ---
  Xi_star_fixed <- cov(Z_star_eff) * (n_eff - 1) / n_eff

  # Use full-sample estimates for the point estimate and gradient
  mu_hat_full <- colMeans(x, na.rm=TRUE)
  sigma_hat_full <- cov(x, use="pairwise.complete.obs")
  # ... (kappa and gradient calculations as before, using full-sample estimates)
  tr_s <- sum(diag(sigma_hat_full)); mu_bar <- mean(mu_hat_full)
  d_mu <- sum((mu_hat_full - mu_bar)^2)
  NF <- sum(sigma_hat_full)-tr_s-d_mu; DF <- (r-1)*tr_s+(r-1)*d_mu; kF <- NF/DF
  NC <- sum(sigma_hat_full)-tr_s; DC <- (r-1)*tr_s+r*d_mu; kC <- NC/DC
  grad_c <- c(2*r*kC,-1,1+(r-1)*kC)/DC; grad_f <- c(2*(1+(r-1)*kF),-1,1+(r-1)*kF)/DF
  grad_matrix <- cbind(grad_f, grad_c)

  scaled_acov <- t(grad_matrix) %*% Xi_star_fixed %*% grad_matrix

  list(
    scaled_acov = scaled_acov,
    kappa_fleiss = kF,
    kappa_cohen = kC,
    n_eff = n_eff
  )
}

kappa_acov_formal_CORRECTED <- function(x) {

  # --- 1. INITIAL SETUP ---
  n <- nrow(x)
  r <- ncol(x)
  M <- !is.na(x)

  mu_hat <- colMeans(x, na.rm = TRUE)
  sigma_hat <- cov(x, use = "pairwise.complete.obs")
  Y_hat <- sweep(x, 2, mu_hat, "-")

  p1 <- colSums(M) / n
  p2 <- crossprod(M) / n

  # --- 2. THE CORRECTED GRADIENT ---
  n_params <- r + r^2
  grad_mat <- matrix(NA, nrow = n_params, ncol = 2)

  # Point estimates
  tr_s <- sum(diag(sigma_hat))
  mu_bar <- mean(mu_hat)
  d_mu <- sum((mu_hat - mu_bar)^2)

  NF <- sum(sigma_hat) - tr_s - d_mu
  DF <- (r-1)*tr_s + (r-1)*d_mu
  kF <- NF/DF

  NC <- sum(sigma_hat) - tr_s
  DC <- (r-1)*tr_s + r*d_mu
  kC <- NC/DC

  # Gradient w.r.t mu (CORRECTED)
  mu_diff_term <- 2 * (mu_hat - mu_bar)
  dkF_dmu <- -mu_diff_term * (DF + NF * (r-1)) / (DF^2)
  dkC_dmu <- -mu_diff_term * (DC * 0 + NC * r) / (DC^2) # d(NC)/dmu = 0

  grad_mat[1:r, 1] <- dkF_dmu
  grad_mat[1:r, 2] <- dkC_dmu

  # Gradient w.r.t sigma_kl (CORRECTED)
  for(k in 1:r) { for(l in 1:r) {
    idx <- r + (k-1)*r + l
    is_diag <- as.numeric(k == l)

    dkF_ds <- (DF * (1 - is_diag) - NF * (r - 1) * is_diag) / DF^2
    dkC_ds <- (DC * (1 - is_diag) - NC * (r - 1) * is_diag) / DC^2

    grad_mat[idx, 1] <- dkF_ds
    grad_mat[idx, 2] <- dkC_ds
  }}

  # --- 3. THE OMEGA_PW MATRIX (Same as before, it was likely correct) ---
  Omega_pw <- matrix(NA, nrow = n_params, ncol = n_params)

  Pi_S <- p2 / (p1 %o% p1)
  Omega_pw[1:r, 1:r] <- Pi_S * sigma_hat

  p3_mat <- array(NA, dim = c(r, r, r))
  for(j in 1:r) for(k in 1:r) for(l in 1:r) p3_mat[j,k,l] <- mean(M[,j]*M[,k]*M[,l])

  for(j in 1:r) { for(k in 1:r) { for(l in 1:r) {
    idx <- (k-1)*r + l
    C_val <- mean(Y_hat[,j]*Y_hat[,k]*Y_hat[,l], na.rm = TRUE)
    Pi_C_val <- p3_mat[j,k,l] / (p1[j] * p2[k,l])
    Omega_pw[j, r + idx] <- Omega_pw[r + idx, j] <- Pi_C_val * C_val
  }}}

  p4_mat <- array(NA, dim=c(r,r,r,r))
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
    scaled_acov = scaled_acov,
    n_eff = n_eff,
    kappa_fleiss = kF,
    kappa_cohen = kC
  )
}



set.seed(135)
n_sim <- 5000; r_sim <- 4; missing_prob <- 0.2
true_mu <- c(1, 1.5, 1.2, 2); true_sigma <- matrix(0.7, r_sim, r_sim); diag(true_sigma) <- 1
L <- chol(true_sigma)
complete_data <- matrix(rnorm(n_sim * r_sim), n_sim, r_sim) %*% L
complete_data <- sweep(complete_data, 2, true_mu, "+")
x_missing <- complete_data
x_missing[sample(1:(n_sim * r_sim), size = n_sim * r_sim * missing_prob)] <- NA
kappa_acov_Z_fixed(x_missing)
kappa_acov_formal_CORRECTED(x_missing)
