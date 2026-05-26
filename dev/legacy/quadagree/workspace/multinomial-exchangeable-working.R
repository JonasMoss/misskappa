#' Pre-compute Constant Matrices for Kappa Calculation
#'
#' This function generates the constant matrices used in the model-based kappa
#' calculation, which depend only on the dimensions of the problem.
#'
#' @param R Integer. The number of raters.
#' @param C Integer. The number of categories.
#' @param loss_matrix A C x C matrix where L_jl is the disagreement/loss
#'        for a pair of ratings (j, l). Defaults to unweighted (0/1) loss.
#' @return A list containing the constant matrices: K_mat (the count matrix),
#'         d_vec (the linear disagreement vector), and Q_mat (the quadratic
#'         chance disagreement matrix).
precompute_kappa_constants <- function(R, C, loss_matrix = NULL) {
  if (is.null(loss_matrix)) {
    loss_matrix <- matrix(1, nrow = C, ncol = C) - diag(C)
  }
  all_k_vectors <- t(partitions::compositions(n = R, m = C, include.zero = TRUE))
  K_mat <- t(all_k_vectors)
  Q_mat <- (1 / R^2) * t(K_mat) %*% loss_matrix %*% K_mat
  d_vec <- apply(K_mat, 2, function(k) {
    quad_term <- t(k) %*% loss_matrix %*% k
    lin_term <- sum(diag(loss_matrix) * k)
    return((quad_term - lin_term) / (R * (R - 1)))
  })
  return(list(K_mat = K_mat, d_vec = d_vec, Q_mat = Q_mat))
}

#' Estimate Model-Based Kappa for Exchangeable Raters
#'
#' Calculates weighted Fleiss' Kappa using a model-based approach, which can
#' be extended to handle missing data. This implementation uses additive
#' smoothing to prevent boundary estimates.
#'
#' @param x A data matrix where rows are subjects and columns are categories.
#'        Each cell x_ic contains the number of raters who assigned category c
#'        to subject i.
#' @param loss_matrix A C x C disagreement/loss matrix. Defaults to unweighted.
#' @param alpha A small numeric value for additive (Lidstone) smoothing to
#'        prevent boundary estimates. Defaults to a near-zero value.
#' @return A list containing the kappa estimate, the estimated phi parameters,
#'         observed and expected agreement, and model constants.
estimate_kappa <- function(x, loss_matrix = NULL, alpha = 1e-9) {
  # --- 1. Input Validation and Setup ---
  x <- as.matrix(x)
  if (any(rowSums(x) != sum(x[1, ]))) {
    stop("All subjects must be rated by the same number of raters.")
  }
  n <- nrow(x)
  C <- ncol(x)
  R <- sum(x[1, ])

  # --- 2. Pre-compute Constants ---
  constants <- precompute_kappa_constants(R, C, loss_matrix)
  K_mat <- constants$K_mat
  d_vec <- constants$d_vec
  Q_mat <- constants$Q_mat
  K_patterns <- ncol(K_mat)

  # --- 3. Count Observed Patterns (Robust String Matching Method) ---
  # This section is completely rewritten for simplicity and robustness.

  # Create unique string identifiers for all possible patterns
  all_pattern_strings <- apply(t(K_mat), 1, paste, collapse = "-")

  # Create unique string identifiers for each row in the observed data
  observed_pattern_strings <- apply(x, 1, paste, collapse = "-")

  # Use table() to count the occurrences of each observed pattern string
  observed_counts_table <- table(observed_pattern_strings)

  # Create the final n_k vector of length K_patterns, ensuring correct order
  n_k <- numeric(K_patterns)
  names(n_k) <- all_pattern_strings

  # Fill the n_k vector with the counts from our table
  n_k[names(observed_counts_table)] <- observed_counts_table

  # --- 4. Estimate Phi using Smoothing ---
  phi_estimate <- (n_k + alpha) / (n + K_patterns * alpha)

  # --- 5. Calculate Kappa ---
  numerator <- sum(d_vec * phi_estimate)
  denominator <- t(phi_estimate) %*% Q_mat %*% phi_estimate

  observed_disagreement <- numerator
  chance_disagreement <- as.numeric(denominator)
  kappa <- 1 - (observed_disagreement / chance_disagreement)

  # --- 6. Return Results ---
  results <- list(
    kappa = kappa,
    phi_estimate = phi_estimate,
    observed_disagreement = observed_disagreement,
    chance_disagreement = as.numeric(denominator),
    observed_agreement = 1 - observed_disagreement,
    chance_agreement = 1 - chance_disagreement,
    params = list(n = n, R = R, C = C, K_patterns = K_patterns, alpha = alpha),
    model_constants = constants
  )
  return(results)
}

#' Compute Kappa and its Variance from Model Parameters
#'
#' This function takes the estimated model parameters (phi) and pre-computed
#' constants to calculate the kappa statistic, its gradient, the variance of
#' the parameters, and the final variance of the kappa estimate.
#'
#' @param phi_estimate A numeric vector of the estimated probabilities for each
#'        count pattern. Must be in the interior of the parameter space (no zeros).
#' @param n Integer. The number of subjects.
#' @param constants A list containing the pre-computed model constants from
#'        `precompute_kappa_constants()`. Must include K_mat, d_vec.
#' @param loss_matrix A C x C disagreement/loss matrix.
#' @return A list containing the kappa estimate, its gradient, the variance-
#'         covariance matrix of phi, and the variance of kappa.
compute_inference <- function(phi_estimate, n, constants, loss_matrix) {
  # --- 1. Extract Constants and Dimensions ---
  K_mat <- constants$K_mat
  d_vec <- constants$d_vec

  R <- sum(K_mat[, 1]) # R can be derived from the first count vector
  C <- nrow(K_mat)
  K_patterns <- ncol(K_mat)

  # --- 2. Calculate Kappa and its Components ---
  # Using the efficient "on-the-fly" method for the denominator
  p_estimate <- (1 / R) * (K_mat %*% phi_estimate)

  observed_disagreement <- sum(d_vec * phi_estimate)
  chance_disagreement <- t(p_estimate) %*% loss_matrix %*% p_estimate

  kappa <- 1 - observed_disagreement / as.numeric(chance_disagreement)

  # --- 3. Compute Gradient of Kappa wrt Phi ---
  # Gradient of numerator (d' * phi) is the constant vector d
  grad_num <- d_vec

  # Gradient of denominator (p' * L * p) using the chain rule
  A_map <- (1/R) * K_mat # Mapping from phi to p
  grad_den <- t(A_map) %*% (loss_matrix + t(loss_matrix)) %*% p_estimate

  # Quotient rule for the gradient of kappa
  num_val <- observed_disagreement
  den_val <- as.numeric(chance_disagreement)

  grad_kappa <- - (grad_num * den_val - num_val * grad_den) / (den_val^2)

  # --- 4. Compute Variance-Covariance Matrix of Phi ---
  # This is the Hessian-based part. The variance of a smoothed multinomial
  # parameter vector phi is given by this formula.
  # The denominator `n` here assumes `phi` was estimated from `n` subjects.
  # If phi came from a model with smoothing parameter alpha, n should be n + K*alpha.
  # We will assume `n` is the effective sample size for simplicity here.
  var_phi <- (diag(as.vector(phi_estimate)) - phi_estimate %*% t(phi_estimate)) / n

  # Note: The inverse of the Hessian of the multinomial log-likelihood is
  # approximately equal to var_phi. This is a standard result.
  # So, var_phi serves as our proxy for the inverse Hessian.

  # --- 5. Compute Kappa Variance via Delta Method ---
  kappa_var <- t(grad_kappa) %*% var_phi %*% grad_kappa
  kappa_var <- as.numeric(kappa_var) # Ensure it's a scalar

  # --- 6. Return All Components ---
  results <- list(
    kappa = kappa,
    grad_kappa = grad_kappa,
    var_phi = var_phi, # This is our proxy for the inverse Hessian
    kappa_var = kappa_var
  )
  return(results)
}

kappa_model <- function(x, loss_matrix = NULL, alpha = 1e-9,
                        conf.level = 0.95, variance = TRUE) {

  x <- as.matrix(x)
  n <- nrow(x)
  C <- ncol(x)
  R <- sum(x[1, ])

  if (is.null(loss_matrix)) {
    loss_matrix <- matrix(1, nrow = C, ncol = C) - diag(C)
  }

  constants <- precompute_kappa_constants(R, C, loss_matrix)
  K_mat <- constants$K_mat
  K_patterns <- ncol(K_mat)

  # Count observed patterns
  all_pattern_strings <- apply(t(K_mat), 1, paste, collapse = "-")
  observed_pattern_strings <- apply(x, 1, paste, collapse = "-")
  observed_counts_table <- table(observed_pattern_strings)
  n_k <- numeric(K_patterns)
  names(n_k) <- all_pattern_strings
  n_k[names(observed_counts_table)] <- observed_counts_table

  # Use smoothing for the phi estimate
  phi_hat_smooth <- (n_k + alpha) / (n + K_patterns * alpha)

  # Default values
  stderr <- NA
  ci <- c(NA, NA)

  # Initial kappa calculation (point estimate)
  p_estimate <- (1 / R) * (K_mat %*% phi_hat_smooth)
  obs_dis <- sum(constants$d_vec * phi_hat_smooth)
  chance_dis <- t(p_estimate) %*% loss_matrix %*% p_estimate
  kappa <- 1 - obs_dis / as.numeric(chance_dis)

  if (variance) {
    # Effective sample size for variance calculation
    n_effective <- n + K_patterns * alpha

    # Call the new inference function
    inference_results <- compute_inference(
      phi_estimate = phi_hat_smooth,
      n = n_effective,
      constants = constants,
      loss_matrix = loss_matrix
    )

    stderr <- sqrt(inference_results$kappa_var)
    z_crit <- qnorm(1 - (1 - conf.level) / 2)
    ci <- c(kappa - z_crit * stderr, kappa + z_crit * stderr)
  }

  results <- list(
    kappa = kappa,
    stderr = stderr,
    conf.int = ci,
    phi_estimate = phi_hat_smooth,
    observed_agreement = 1 - obs_dis,
    chance_agreement = 1 - as.numeric(chance_dis),
    params = list(n=n, R=R, C=C, alpha=alpha, conf.level=conf.level),
    # Optionally return the detailed inference components
    inference_details = if(variance) inference_results else NULL
  )
  return(results)
}








# --- Verification Script ---

# The dataset 'dat.fleiss1971' is assumed to be present in the environment.
# It is a 30x5 matrix of counts.
x_data <- as.matrix(dat.fleiss1971)

# 1. Run our estimator
our_model_results <- kappa_model(x = x_data, alpha = 0)

# 2. Run the reference package
irr_results <- irrCAC::fleiss.kappa.dist(x_data)

# 3. Compare the results
cat("--- Kappa Point Estimate Comparison ---\n")
cat(sprintf("Our Model-Based Kappa:      %.7f\n", our_model_results$kappa))
cat(sprintf("irrCAC Fleiss' Kappa:       %.7f\n\n", irr_results$coeff))

# (The rest of the comparison printouts are the same...)
