# Your R analysis function (with the corrections from before)
fleiss_aggr_fun_corrected <- function(calc) {
  is_effective_row <- rowSums(abs(calc$xx)) > 1e-9
  xx_eff <- calc$xx[is_effective_row, , drop = FALSE]
  n_eff <- nrow(xx_eff)
  if (n_eff < 2) {
    return(list(est = NA, var = NA))
  }

  means <- colMeans(xx_eff)
  theta <- cov(xx_eff) * (n_eff - 1) / n_eff # Use method="ML" for biased cov

  phi1 <- means[1]
  phi2 <- means[2]
  phi3 <- means[3]
  R_total <- calc$r

  N_term <- phi2 - phi1^2
  D_term <- phi3 - phi1^2 / R_total
  if (abs(D_term) < 1e-9) {
    return(list(est = NA, var = NA))
  }

  est <- (1 / (R_total - 1)) * (N_term / D_term - 1)

  dk_dphi1 <- (1 / (R_total - 1)) * (-2 * phi1 / D_term + (N_term * (2 * phi1 / R_total)) / (D_term^2))
  dk_dphi2 <- (1 / (R_total - 1)) / D_term
  dk_dphi3 <- (1 / (R_total - 1)) * (-N_term / (D_term^2))
  grad <- c(dk_dphi1, dk_dphi2, dk_dphi3)

  scaled_var <- c(t(grad) %*% theta %*% grad)
  final_var <- scaled_var / n_eff

  list(est = est, var = final_var)
}

# Your R preparation function for DISTRIBUTION (count) data
fleiss_aggr_prepare_dist_missing <- function(x_counts, R, values = NULL) {
  y <- as.matrix(x_counts)
  k_cat <- ncol(y)
  if (is.null(values)) values <- 1:k_cat

  r_i <- rowSums(y)
  weights <- ifelse(r_i > 0, R / r_i, 0)

  s_obs <- c(y %*% values)
  q_obs <- c(y %*% (values^2))

  s_i <- weights * s_obs
  q_i <- weights * q_obs

  list(
    xx = cbind(s_i, s_i^2, q_i),
    r = R
  )
}

# --- The Definitive Test ---
set.seed(123)
values_vec <- 1:4
R_total <- 4

# Generate raw subject-by-rater data with NAs
raw_data <- matrix(sample(values_vec, 50 * R_total, replace = TRUE), nrow = 50, ncol = R_total)
raw_data[sample(1:(50 * R_total), 40)] <- NA

# --- THIS IS THE CRITICAL STEP ---
# Convert the raw data into the correct n x k count matrix format
dist_data <- t(apply(raw_data, 1, function(row) {
  table(factor(row, levels = values_vec))
}))


# --- Test the R Prototype ---
# 1. Prepare data using the correct function for count data
prepared_data_r <- fleiss_aggr_prepare_dist_missing(dist_data, R = R_total, values = values_vec)
# 2. Run the analysis function
results_r <- fleiss_aggr_fun_corrected(prepared_data_r)

# --- Test the Rcpp Function ---
# 3. Run the Rcpp function on the SAME count data
#results_rcpp <- fleissci_aggr_rcpp(dist_data, values_vec, R_total)
#results_rcpp <- results_rcpp$sd^2

# --- Compare ---
cat("--- R Prototype Result (Distribution Form) ---\n")
print(results_r)

cat("\n--- Rcpp Result (Distribution Form) ---\n")
#print(results_rcpp)

kappa_aggr(dist_data, values_vec, R_total)
