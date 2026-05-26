

# We need to pass R as an argument AND have it default to the maximal number of rated items.
  # Note: This function now requires the raw (or semi-raw) subject-by-category matrix
  # to calculate r_i for each subject.
fleiss_aggr_missing <- function(x, values = seq_len(ncol(x))) {

  # --- Setup ---
  # x is now the n x k matrix of counts for each subject
  y <- as.matrix(x)
  R <- sum(x[1, ]) # Total possible raters (assumes first subject is representative)
  # A more robust way is to pass R as an argument.
  # Let's assume R is the true total number of raters/columns.
  R <- ncol(values) # This is safer if values represents all possible raters.
  # Let's assume x is now a subject-by-rater matrix with NAs.

  # --- This is the key change ---
  # We now assume x is a n_subjects x R_raters matrix with NAs
  n <- nrow(y)
  R_total <- ncol(y) # The total number of possible raters

  # Calculate per-subject weights
  r_i <- rowSums(!is.na(y))
  # Avoid division by zero for subjects with no ratings
  weights <- ifelse(r_i > 0, R_total / r_i, 0)

  # Set NA to 0 for matrix multiplication, since they won't be counted anyway
  y_zeroed <- y
  y_zeroed[is.na(y)] <- 0

  # --- Calculate OBSERVED sums using vectorized operations ---
  # These are S_obs,i and Q_obs,i
  s_obs <- y_zeroed %*% values
  q_obs <- y_zeroed %*% (values^2)

  # --- Apply weights to get the reweighted statistics s_i and q_i ---
  s_i <- weights * s_obs
  q_i <- weights * q_obs

  # --- Estimate the phi parameters using means of reweighted stats ---
  phi1_hat <- mean(s_i)
  phi2_hat <- mean(s_i^2)
  phi3_hat <- mean(q_i)

  # --- Final kappa calculation using the same formula as before ---
  # The denominator D from your original derivation
  D_denom <- phi3_hat - (1/R_total) * phi1_hat^2

  # The numerator N
  N_num <- phi2_hat - phi3_hat - (1 - 1/R_total) * phi1_hat^2

  kappa_est <- (1 / (R_total - 1)) * (N_num / D_denom)

  return(as.numeric(kappa_est))
}




fleiss_aggr <- \(x, values = seq_len(ncol(x))) {
  r <- sum(x[1, ])
  stopifnot(ncol(x) == length(values))

  y <- as.matrix(x)
  xtx <- tcrossprod(values^2, y)
  xt1 <- tcrossprod(values, y)

  extx <- mean(xtx)
  ext1 <- mean(xt1)
  ext2 <- mean(xt1^2)

  1 / (r - 1) * ((ext2 - ext1^2) / (extx - ext1^2 / r) - 1)
}


fleiss_aggr_prepare <- \(x, values) {
  y <- as.matrix(x)
  r <- sum(y[1, ])
  xtx <- c(tcrossprod(values^2, y))
  xt1 <- c(tcrossprod(values, y))
  xt12 <- xt1^2

  list(
    xx = cbind(xt1, xt12, xtx),
    n = nrow(y),
    r = r
  )
}

fleiss_aggr_fun <- \(calc) {
  means <- colMeans(calc$xx)
  theta <- tcrossprod(t(calc$xx) - means) / calc$n
  k <- calc$r / (means[3] * calc$r - means[1]^2)
  km <- k * (means[2] - means[1]^2)
  calc_r_inv <- 1 / (calc$r - 1)

  grad <- k * calc_r_inv * c(2 * means[1] * (km / calc$r - 1), 1, -km)
  est <- calc_r_inv * (km - 1)
  var <- c(crossprod(grad, theta %*% grad))

  list(est = est, var = var)
}
