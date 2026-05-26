# Your original, elegant functions (for reference)
fleiss_aggr_fun <- function(calc) {
  means <- colMeans(calc$xx)
  # Corrected covariance calculation for clarity
  theta <- cov(calc$xx) * (nrow(calc$xx) - 1) / nrow(calc$xx)

  # Your kappa formula re-expressed in terms of the means
  # D_denom = means[3] - (1/calc$r)*means[1]^2
  # N_num = means[2] - means[3] - (1 - 1/calc$r)*means[1]^2
  # est_alt = (1/(calc$r-1)) * (N_num/D_denom)

  # Let's use your original, more compact form for the gradient
  k <- calc$r / (means[3] * calc$r - means[1]^2)
  km <- k * (means[2] - means[1]^2)
  calc_r_inv <- 1 / (calc$r - 1)

  # Your original gradient from the prompt
  # Note: This gradient is w.r.t [phi1, phi2, phi3], not [phi1, phi1^2, phi3]
  # We need to be careful here. Let's use the verified one from the aggregate method.

  # Let's recalculate the gradient wrt means = [mean(s), mean(s^2), mean(q)]
  phi1_h <- means[1]; phi2_h <- means[2]; phi3_h <- means[3]
  R <- calc$r

  D_denom <- phi3_h - (1/R) * phi1_h^2 # This isn't right, phi2 is mean(s^2), not (mean(s))^2
  # The kappa formula expects E[S], E[S^2], E[Q]
  # Your 'calc$xx' provides s_i, s_i^2, q_i
  # So means[1]=E[s_i], means[2]=E[s_i^2], means[3]=E[q_i]
  # This is correct.

  D_denom <- means[3] - (1/R)*means[1]^2 # This assumes E[S^2] = (E[S])^2 which is only true if Var(S)=0.

  # Let's use the explicit kappa formula, which is clearer
  # phi1 = E[S], phi2 = E[S^2], phi3 = E[Q]
  phi1 <- means[1]; phi2 <- means[2]; phi3 <- means[3]
  R_total <- calc$r

  # The kappa formula is in terms of Var(S) and E[S]
  # Var(S) = E[S^2] - (E[S])^2 = phi2 - phi1^2
  # 1'S1 = Var(S) + (E[S])^2 = phi2

  # kappa_F = ( (1'S1) - tr(S) - D_mu ) / ( (r-1)tr(S) + (r-1)D_mu )
  # This is the raw data formulation. Let's stick to your distribution form one.
  # k_F = 1/(R-1) * ( (phi2 - phi1^2) / (phi3 - phi1^2/R) - 1 )

  N_term <- phi2 - phi1^2
  D_term <- phi3 - phi1^2 / R_total
  est <- (1/(R_total-1)) * (N_term / D_term - 1)

  # Gradient of est w.r.t [phi1, phi2, phi3]
  dk_dphi1 <- (1/(R_total-1)) * ( -2*phi1/D_term + (N_term * (2*phi1/R_total))/(D_term^2) )
  dk_dphi2 <- (1/(R_total-1)) / D_term
  dk_dphi3 <- (1/(R_total-1)) * (-N_term / (D_term^2))

  grad <- c(dk_dphi1, dk_dphi2, dk_dphi3)

  var <- c(t(grad) %*% theta %*% grad)

  list(est = est, var = var / nrow(calc$xx))
}

# --- The NEW preparation function for missing data ---
fleiss_aggr_prepare_missing <- function(x, values = seq_len(ncol(x))) {
  # Assumes x is a n_subjects x R_raters matrix with NAs
  y <- as.matrix(x)
  n <- nrow(y)
  R_total <- ncol(y)

  # Calculate per-subject weights under exchangeability assumption
  r_i <- rowSums(!is.na(y))
  weights <- ifelse(r_i > 0, R_total / r_i, 0)

  # Use 0 for NAs in matrix multiplication
  y_zeroed <- y
  y_zeroed[is.na(y)] <- 0

  # Calculate OBSERVED sums (S_obs,i and Q_obs,i)
  s_obs <- c(y_zeroed %*% values)
  q_obs <- c(y_zeroed %*% (values^2))

  # Apply weights to get the reweighted statistics
  s_i <- weights * s_obs
  q_i <- weights * q_obs

  list(
    # The output is a matrix with columns [s_i, s_i^2, q_i]
    # This matches the structure needed by fleiss_aggr_fun
    xx = cbind(s_i, s_i^2, q_i),
    r = R_total
  )
}

# --- Example Usage ---
set.seed(123)
# Assume values are categories 1, 2, 3, 4
values_vec <- 1:4
# Generate raw subject-by-rater data with NAs
x_missing <- matrix(sample(values_vec, 50*4, replace=TRUE), nrow=50, ncol=4)
x_missing[sample(1:200, 40)] <- NA


# 1. Prepare the reweighted statistics
prepared_data <- fleiss_aggr_prepare_missing(x_missing, values = values_vec)

# 2. Run the analysis function on the prepared data
results <- fleiss_aggr_fun(prepared_data)

print(results)

kappa_dist_rcpp(x_missing, values_vec, 4)
