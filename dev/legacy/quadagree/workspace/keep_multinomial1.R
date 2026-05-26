#' MINIMALLY MODIFIED version of compute_kappa
#'
#' This function is identical to the original `compute_kappa`, but with one single
#' change: it replaces the variance-covariance matrix derived from the EM
#' algorithm's information matrix with the known, correct variance-covariance
#' matrix for a complete multinomial data MLE.
#'
#' This allows for a direct comparison to isolate errors in the variance calculation.
#'
#' @param x A matrix of count data.
#' @param raters Integer. The number of raters.
#' @param loss_matrix A C x C disagreement/loss matrix.
#'
compute_kappa_comp_var <- function(x, raters, loss_matrix = NULL) {

  # --- NO CHANGES in this section ---
  x <- as.matrix(x)
  n <- nrow(x)
  n_cat <- ncol(x)
  fit <- em_counts(x, raters) # We still run EM to get the correct MLE (phi)

  if (is.null(loss_matrix)) {
    loss_matrix <- matrix(1, nrow = n_cat, ncol = n_cat) - diag(n_cat)
  }

  constants <- kappa_constants(raters, n_cat, loss_matrix)

  phi_estimate <- fit$theta
  positions <- fit$positions

  # ********************************************************************
  # *                THE ONLY MODIFICATION IS HERE                     *
  # ********************************************************************
  #
  # We DISCARD the information matrix from the EM fit:
  #   information <- fit$information
  #   var_phi <- MASS::ginv(as.matrix(information))
  #
  # And REPLACE it with the known formula for complete data variance:
  var_phi <- (1/n) * (diag(phi_estimate) - tcrossprod(phi_estimate))
  #
  # ********************************************************************

  # --- NO CHANGES in this section ---
  K_mat <- constants$k_mat[, positions]
  d_vec <- constants$d_vec[positions]
  p_estimate <- (1 / raters) * (K_mat %*% phi_estimate)

  observed_disagreement <- sum(d_vec * phi_estimate)
  chance_disagreement <- t(p_estimate) %*% loss_matrix %*% p_estimate

  kappa <- 1 - observed_disagreement / as.numeric(chance_disagreement)

  # --- Gradient Calculation (with the factor of 2 fix) ---
  grad_num <- d_vec
  A_map <- (1/raters) * K_mat

  # Including the factor of 2 fix here as it's a known calculus error
  grad_den <- 2 * t(A_map) %*% loss_matrix %*% p_estimate

  num_val <- observed_disagreement
  den_val <- as.numeric(chance_disagreement)
  grad_kappa <- - (grad_num * den_val - as.numeric(num_val) * grad_den) / (den_val^2)

  # --- Final variance calculation using the SUBSTITUTED var_phi ---
  kappa_var <- t(grad_kappa) %*% var_phi %*% grad_kappa
  kappa_var <- as.numeric(kappa_var)

  list(
    kappa = kappa,
    grad_kappa = grad_kappa,
    var_phi = var_phi,
    kappa_var = kappa_var
  )
}

# Load data
data(dat.fleiss1971, package = "irrCAC")
x <- as.matrix(dat.fleiss1971)
n <- nrow(x)
raters <- 6

# Your original function's output
original_results <- compute_kappa(x, raters)
print(paste("Original Stderr:", sqrt(original_results$kappa_var)))

# The minimally modified function's output
comp_var_results <- compute_kappa_comp_var(x, raters)
print(paste("Comp Var Stderr:", sqrt(comp_var_results$kappa_var)))

# irrCAC for comparison (n-based)
irrCAC_stderr_n <- irrCAC::fleiss.kappa.dist(x)$stderr * sqrt(n - 1) / sqrt(n)
print(paste("irrCAC Stderr:  ", irrCAC_stderr_n))
