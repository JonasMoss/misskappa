#' Pre-compute Constant Matrices for Kappa Calculation
#'
#' This function generates the constant matrices used in the model-based kappa
#' calculation, which depend only on the dimensions of the problem.
#'
#' @param raters Integer. The number of raters.
#' @param n_cat Integer. The number of categories.
#' @param loss_matrix A n_cat x n_cat matrix where L_jl is the disagreement/loss
#'        for a pair of ratings (j, l). Defaults to unweighted (0/1) loss.
#' @return A list containing the constant matrices: K_mat (the count matrix),
#'         d_vec (the linear disagreement vector), and Q_mat (the quadratic
#'         chance disagreement matrix).
kappa_constants <- function(raters, n_cat, loss_matrix = NULL) {
  if (is.null(loss_matrix)) {
    loss_matrix <- matrix(1, nrow = n_cat, ncol = n_cat) - diag(n_cat)
  }
  all_k_vectors <- t(partitions::compositions(n = raters, m = n_cat, include.zero = TRUE))
  k_mat <- t(all_k_vectors)
  d_mat <- (1 / raters^2) * t(k_mat) %*% loss_matrix %*% k_mat
  d_vec <- apply(k_mat, 2, function(k) {
    quad_term <- t(k) %*% loss_matrix %*% k
    lin_term <- sum(diag(loss_matrix) * k)
    ((quad_term - lin_term) / (raters * (raters - 1)))
  })

  list(k_mat = k_mat,
       d_vec = d_vec,
       d_mat = d_mat)
}

#' Em algorithm on count data
#' @param x Count data.
#' @param raters Number of raters.
#' @param eps Tolereance for EM algorithm.
#' @param maxit Maximum number of iteratios.
#'
em_counts <- function(x, raters, eps = 1e-8, maxit = 1000) {

  state_space <- function(raters, C) {
    grd <- expand.grid(rep(list(0:raters), C))
    as.matrix(grd[rowSums(grd) == raters, , drop = FALSE])
  }

  if (any(is.na(x))) x[is.na(x)] <- 0

  n <- nrow(x);
  n_cat <- ncol(x)
  S <- state_space(raters, n_cat);
  Sdim <- nrow(S)

  theta <- rep(1 / Sdim, Sdim)

  compat <- lapply(
    1:n,
    function(i) which(apply(S >= matrix(x[i, ], nrow = Sdim, ncol = n_cat,
                                        byrow = TRUE), 1, all))
  )

  for (it in 1:maxit) {
    tau <- matrix(0, n, Sdim)
    for (i in 1:n) {
      idx        <- compat[[i]]
      w          <- theta[idx]
      tau[i,idx] <- w / sum(w)                 # E-step
    }
    new_theta <- colMeans(tau)                 # M-step
    if (max(abs(new_theta - theta)) < eps) break
    theta <- new_theta
  }

  eps  <- .Machine$double.eps
  positions <- which(theta > eps)
  theta_pos <- theta[positions]
  tau_pos   <- tau[, positions]

  # N_pos <- colSums(tau_pos)
  # information  <- Matrix::crossprod(Matrix::Matrix(tau_pos, sparse = TRUE)) /
  #   tcrossprod(theta_pos)

  S_matrix <- tau_pos / matrix(theta_pos, nrow = nrow(tau_pos), ncol = ncol(tau_pos), byrow = TRUE)
  information <- crossprod(S_matrix)


  list(theta = theta_pos,
       positions = positions,
       information  = information,
       iter   = it)
}

#' EM algorithm on count data with Constrained Variance Calculation
#'
#' This function performs the EM algorithm to estimate multinomial probabilities
#' from potentially incomplete count data.
#'
#' It calculates the variance-covariance matrix of the estimates using Louis's
#' method, correctly accounting for the sum-to-one constraint on the
#' probabilities. This is achieved by re-parameterizing to an unconstrained
#' space of m-1 parameters, calculating the variance there, and then
#' transforming back to the full m-dimensional space via the delta method.
#'
#' @param x A matrix of count data where rows are subjects and columns are categories.
#' @param raters Integer. The total number of raters for each subject.
#' @param eps Tolerance for EM algorithm convergence.
#' @param maxit Maximum number of iterations for the EM algorithm.
#'
#' @return A list containing:
#'         - theta: The vector of estimated probabilities for active patterns.
#'         - positions: The indices of the active patterns in the full state space.
#'         - var_phi: The correct, non-diagonal, m x m variance-covariance matrix.
#'         - iter: The number of iterations the EM algorithm took to converge.
#'
#' EM algorithm on count data with Constrained Variance Calculation
#'
#' (Full documentation from before...)
#'
em_counts_constrained <- function(x, raters, eps = 1e-8, maxit = 1000) {

  # --- This part is identical to your original function ---
  state_space <- function(raters, C) {
    grd <- expand.grid(rep(list(0:raters), C))
    as.matrix(grd[rowSums(grd) == raters, , drop = FALSE])
  }

  if (any(is.na(x))) x[is.na(x)] <- 0

  n <- nrow(x)
  n_cat <- ncol(x)
  S <- state_space(raters, n_cat)
  Sdim <- nrow(S)

  theta <- rep(1 / Sdim, Sdim)

  compat <- lapply(
    1:n,
    function(i) {
      which(apply(S >= matrix(x[i, ], nrow = Sdim, ncol = n_cat, byrow = TRUE), 1, all))
    }
  )

  for (it in 1:maxit) {
    tau <- matrix(0, n, Sdim)
    for (i in 1:n) {
      idx        <- compat[[i]]
      w          <- theta[idx]
      sum_w      <- sum(w)
      if (sum_w > 0) {
        tau[i, idx] <- w / sum_w
      }
    }
    new_theta <- colMeans(tau)
    if (max(abs(new_theta - theta)) < eps) break
    theta <- new_theta
  }

  # --- Get the active set ---
  eps_pos <- .Machine$double.eps
  positions <- which(theta > eps_pos)
  theta_pos <- theta[positions]
  tau_pos   <- tau[, positions, drop = FALSE]
  m_pos     <- length(theta_pos)

  if (m_pos <= 1) {
    return(list(
      theta = theta_pos,
      positions = positions,
      var_phi = matrix(0, nrow = m_pos, ncol = m_pos),
      iter = it
    ))
  }

  ref_idx <- m_pos
  non_ref_indices <- 1:(m_pos - 1)

  main_scores <- sweep(tau_pos[, non_ref_indices, drop = FALSE], 2, theta_pos[non_ref_indices], "/")

  ref_scores_vec <- tau_pos[, ref_idx] / theta_pos[ref_idx]
  s_star_matrix <- main_scores - ref_scores_vec


  info_star <- crossprod(s_star_matrix)

  var_star <- tryCatch(solve(info_star), error = function(e) {
    warning("Information matrix for reduced parameters was singular; using ginv().")
    MASS::ginv(info_star)
  })

  J <- matrix(0, nrow = m_pos, ncol = m_pos - 1)
  J[non_ref_indices, ] <- diag(m_pos - 1)
  J[ref_idx, ] <- -1

  var_phi_correct <- J %*% var_star %*% t(J)

  list(
    theta = theta_pos,
    positions = positions,
    var_phi = var_phi_correct,
    iter = it
  )
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
compute_kappa <- function(x, raters, loss_matrix = NULL) {

  x <- as.matrix(x)
  n <- nrow(x)
  n_cat <- ncol(x)
  fit <- em_counts(x, raters)

  if (is.null(loss_matrix)) {
    loss_matrix <- matrix(1, nrow = n_cat, ncol = n_cat) - diag(n_cat)
  }

  constants <- kappa_constants(raters, n_cat, loss_matrix)

  phi_estimate <- fit$theta
  positions <- fit$positions
  information <- fit$information
  var_phi <- MASS::ginv(as.matrix(information))

  K_mat <- constants$k_mat[, positions]
  d_vec <- constants$d_vec[positions]
  p_estimate <- (1 / raters) * (K_mat %*% phi_estimate)

  observed_disagreement <- sum(d_vec * phi_estimate)
  chance_disagreement <- t(p_estimate) %*% loss_matrix %*% p_estimate

  kappa <- 1 - observed_disagreement / as.numeric(chance_disagreement)
  # This kappa is correct with complete data.




  ## Old grads.
  grad_num <- d_vec
  A_map <- (1/raters) * K_mat


  grad_den <- t(A_map) %*% (loss_matrix) %*% p_estimate
  num_val <- observed_disagreement
  den_val <- as.numeric(chance_disagreement)
  grad_kappa <- - (grad_num * den_val - 2 * num_val * grad_den) / (den_val^2)

  kappa_var <- t(grad_kappa) %*% var_phi %*% grad_kappa
  kappa_var <- as.numeric(kappa_var) # Ensure it's a scalar

  # --- 6. Return All Components ---
  results <- list(
    kappa = kappa,
    grad_kappa = grad_kappa,
    var_phi = var_phi,
    kappa_var = kappa_var
  )
  return(results)
}

# Rename your original compute_kappa to use the new EM function
compute_kappa_final <- function(x, raters, loss_matrix = NULL) {

  x <- as.matrix(x)
  n <- nrow(x)
  n_cat <- ncol(x)

  # **** THE ONLY CHANGE IS HERE ****
  fit <- em_counts_constrained(x, raters)
  # **********************************

  if (is.null(loss_matrix)) {
    loss_matrix <- matrix(1, nrow = n_cat, ncol = n_cat) - diag(n_cat)
  }

  constants <- kappa_constants(raters, n_cat, loss_matrix)

  phi_estimate <- fit$theta
  positions <- fit$positions
  var_phi <- fit$var_phi # Use the corrected variance matrix

  K_mat <- constants$k_mat[, positions]
  d_vec <- constants$d_vec[positions]
  p_estimate <- (1 / raters) * (K_mat %*% phi_estimate)

  observed_disagreement <- sum(d_vec * phi_estimate)
  chance_disagreement <- t(p_estimate) %*% loss_matrix %*% p_estimate

  kappa <- 1 - observed_disagreement / as.numeric(chance_disagreement)

  # Gradient with the factor of 2 fix
  grad_num <- d_vec
  A_map <- (1/raters) * K_mat
  grad_den <- 2 * t(A_map) %*% loss_matrix %*% p_estimate
  num_val <- observed_disagreement
  den_val <- as.numeric(chance_disagreement)
  grad_kappa <- - (grad_num * den_val - as.numeric(num_val) * grad_den) / (den_val^2)

  kappa_var <- t(grad_kappa) %*% var_phi %*% grad_kappa
  kappa_var <- as.numeric(kappa_var)

  list(
    kappa = kappa,
    kappa_var = kappa_var
  )
}

#' Calculate Fleiss' Kappa with Support for Missing Data
#'
#' @description
#' Provides a robust, likelihood-based estimate of Fleiss' Kappa and its
#' confidence intervals. This implementation uses an Expectation-Maximization (EM)
#' algorithm to find the maximum likelihood estimate of the distribution of
#' rating patterns, which naturally handles datasets where some subjects have
#' been rated by fewer than the total number of raters.
#'
#' @param x An N x c matrix or data.frame of counts, where N is the number
#'   of subjects and c is the number of categories.
#' @param weights A character string specifying the weighting scheme (one of
#'   "unweighted", "linear", "quadratic") or a single custom c x c disagreement
#'   (loss) matrix. Defaults to "unweighted".
#' @param conf_level The confidence level for the returned intervals.
#' @param control A list of advanced control options. See Details.
#'
#' @details
#' \strong{Methodology}:
#' The core of this function is a likelihood-based approach to agreement. It assumes
#' the observed rating patterns are incomplete observations from a single multinomial
#' distribution over all possible complete rating patterns. The function estimates
#' the parameters of this distribution using an EM algorithm, as described by
#' Dempster, Laird, and Rubin (1977). This approach is valid under the Missing
#' at Random (MAR) assumption.
#'
#' The variance of the kappa estimate is calculated using the delta method, with the
#' variance-covariance matrix of the multinomial parameters estimated via Louis's
#' (1982) method for observed information from the EM algorithm.
#'
#' \strong{Control Options}:
#' The `control` list provides advanced control over the analysis:
#' \itemize{
#'   \item \code{r}: The total number of raters in the study. If NULL (the
#'     default), it is automatically set to the maximum number of ratings
#'     observed for any single subject.
#'   \item \code{category_values}: A numeric vector of values for each
#'     category, required for "linear" and "quadratic" weights. Defaults to `1:c`.
#'   \item \code{transform}: The variance-stabilizing transformation for the
#'     asymptotic CI. Defaults to "fisher", which is recommended for kappa.
#'     Can be one of "fisher", "none", "log", "arcsin".
#'   \item \code{bootstrap}: The type of bootstrap to perform. Defaults to "none".
#'     Can be one of "none", "nonparametric", "parametric", "nonparametric-t" (studentized),
#'     or "parametric-t" (studentized). See Efron (1994) for a discussion of
#'     bootstrapping with missing data.
#'   \item \code{bootstrap_reps}: The number of bootstrap replications. Defaults
#'     to 1000.
#'   \item \code{seed}: An integer to seed the random number generator for
#'     reproducible bootstrap results.
#' }
#'
#' @return
#' An object of class `fleiss_kappa` containing two main components:
#' \item{estimates}{A data frame with the kappa estimate, its standard error,
#'   and asymptotic and (if requested) bootstrap confidence intervals.}
#' \item{diagnostics}{A list containing details from the estimation process,
#'   including EM convergence status, iterations, the final parameter vector
#'   (`theta`), its variance-covariance matrix, and the full set of bootstrap
#'   replicates if a bootstrap was performed.}
#'
#' @references
#' Dempster, A. P., Laird, N. M., & Rubin, D. B. (1977). Maximum Likelihood from
#'   Incomplete Data via the EM Algorithm. *Journal of the Royal Statistical
#'   Society: Series B (Methodological), 39*(1), 1–22.
#'
#' Efron, B. (1994). Missing Data, Imputation, and the Bootstrap.
#'   *Journal of the American Statistical Association, 89*(426), 463–475.
#'
#' Fleiss, J. L. (1971). Measuring nominal scale agreement among many raters.
#'   *Psychological Bulletin, 76*(5), 378–382.
#'
#' Schafer, J. L. (1997). *Analysis of Incomplete Multivariate Data*.
#'   Chapman and Hall/CRC. (See Chapter 4 for bootstrap methods with missing data).
#'
#' @export
#' @examples
#' # --- Basic Usage ---
#' # Create a sample dataset with missing data (6 total raters)
#' set.seed(123)
#' x_complete <- t(rmultinom(30, size = 6, prob = c(0.4, 0.3, 0.1, 0.1, 0.1)))
#' x_missing <- x_complete
#' x_missing[1:5, ] <- t(rmultinom(5, size = 5, prob = c(0.4, 0.3, 0.1, 0.1, 0.1)))
#' x_missing[6:10, ] <- t(rmultinom(5, size = 4, prob = c(0.4, 0.3, 0.1, 0.1, 0.1)))
#'
#' # Default call with unweighted kappa and Fisher-transformed CI
#' fleiss_kappa(x_missing)
#'
#' # --- Advanced Usage ---
#' # Quadratic weights for ordered categories
#' fleiss_kappa(x_missing, weights = "quadratic")
#'
#' # Non-parametric bootstrap with 500 replicates for more robust CIs
#' \dontrun{
#'   fleiss_kappa(x_missing, bootstrap = "nonparametric", bootstrap_reps = 500, seed = 42)
#' }
#'
#' # Specifying `r` manually for the theoretical shrinkage case
#' # (Note: This is for demonstration; `r` should typically be the true total raters)
#' fleiss_kappa(x_complete, control = list(r = 8))
#'
fleiss_kappa <- function(
    x,
    weights = "unweighted",
    conf_level = 0.95,
    control = list()
) {

  # --- 1. Argument Validation and Control List Handling ---
  if (!is.matrix(x) && !is.data.frame(x)) stop("`x` must be a matrix or data.frame.")
  x_mat <- as.matrix(x)
  if (!is.numeric(x_mat)) stop("`x` must be numeric.")
  storage.mode(x_mat) <- "integer"

  # Define defaults and merge with user-provided control list
  ctrl <- list(
    r = NULL,
    category_values = NULL,
    transform = "fisher",
    bootstrap = "none",
    bootstrap_reps = 1000,
    seed = NULL
  )
  # Overwrite defaults with any user-supplied values
  ctrl[names(control)] <- control

  # Use match.arg for robust argument checking
  transform_type <- match.arg(tolower(ctrl$transform), c("fisher", "none", "log", "arcsin"))
  bootstrap_method <- match.arg(tolower(ctrl$bootstrap), c("none", "nonparametric", "parametric", "nonparametric-t", "parametric-t"))

  # --- 2. Prepare Model Parameters ---
  if (is.null(ctrl$r)) {
    r_model <- max(rowSums(x_mat, na.rm = TRUE))
    message(paste("`r` not specified in control list. Setting to max observed raters:", r_model))
  } else {
    r_model <- as.integer(ctrl$r)
  }

  c <- ncol(x_mat)
  category_values <- if (is.null(ctrl$category_values)) 1:c else ctrl$category_values

  # --- 3. Build Loss Matrix (singular, for simplicity) ---
  weight_keys <- character()
  loss_matrices <- list()

  if (is.matrix(weights)) {
    if(nrow(weights) != c || ncol(weights) !=c) stop("Custom weight matrix must be c x c.")
    weight_keys <- "custom"
    loss_matrices <- list(weights)
  } else {
    weight_type <- match.arg(weights, c("unweighted", "linear", "quadratic"))
    weight_keys <- weight_type
    loss_matrices <- list(generate_loss_matrix_rcpp(weight_type, c, category_values))
  }

  # --- 4. Build Configuration List for C++ ---
  config <- list(
    x_data = x_mat,
    r_model = r_model,
    c = c,
    weight_keys = weight_keys,
    loss_matrices = loss_matrices,
    conf_level = conf_level,
    transform_type = transform_type,
    bootstrap_method = bootstrap_method,
    bootstrap_reps = as.integer(ctrl$bootstrap_reps),
    seed = if (!is.null(ctrl$seed)) as.integer(ctrl$seed) else NULL
  )

  # --- 5. Call the Rcpp function ---
  result_list <- run_analysis_rcpp(config)

  # --- 6. Format and Return Output ---
  names(result_list) <- c("estimates", "diagnostics")
  class(result_list) <- "fleiss_kappa"
  return(result_list)
}

# Define a new print method for the new class name
print.fleiss_kappa <- function(x, ...) {
  cat("Fleiss' Kappa with Missing Data\n\n")
  cat("Agreement Estimates:\n")
  print(x$estimates, row.names = FALSE, digits=3)
  cat("\nEM Diagnostics:\n")
  if(!is.null(x$diagnostics$iterations)) {
    cat(paste("  Converged:", x$diagnostics$converged, "in", x$diagnostics$iterations, "iterations\n"))
  }
}

