set.seed(313)
x <- as.matrix(dat.zapf2016)
#x[sample(length(x), 40)] <- NA
storage.mode(x) <- "integer"
#x <- x[20:30, ]
x

# -------------------------------------------------------------------
# Step 1: Preprocess the data
# -------------------------------------------------------------------
prep_results <- preprocess_conger_cpp(x)

cat("\n--- Preprocessing Results ---\n")
str(prep_results)


# -------------------------------------------------------------------
# Step 2: Define the full configuration list
# -------------------------------------------------------------------
# This now matches the structure your C++ code expects, even if not all
# parameters are used by the Conger's kappa EM function yet.
config <- list(
  # Data-related (not used by run_em_conger_cpp, but good for completeness)
  x_data = x,
  r_model = NA_integer_, # Not applicable for Conger's
  c = ncol(x),           # Number of categories/raters (depends on context)

  # Weighting-related (will be used in the final kappa calculation step)
  weight_keys = c("unweighted", "linear", "quadratic"),
  loss_matrices = list(
    generate_loss_matrix_rcpp("unweighted", prep_results$K_categories, values = prep_results$original_categories),
    generate_loss_matrix_rcpp("linear", prep_results$K_categories, values = prep_results$original_categories),
    generate_loss_matrix_rcpp("quadratic", prep_results$K_categories, values = prep_results$original_categories)
  ),

  # EM Algorithm parameters
  em_tol = 1e-8,
  max_iter = 10000,
  em_prune_tol = 1e-9,
  start_alpha = 0.001,

  # Bootstrap and CI parameters (will be used by the final analysis function)
  seed = 123,
  conf_level = 0.95,
  transform_type = "none", # e.g., "none", "fisher", "log", "arcsin"
  alternative_hypothesis = "two.sided", # "two.sided", "greater", "less"
  bootstrap_method = "none", # "none", "nonparametric", etc.
  bootstrap_reps = 1000
)

cat("\n--- Full Configuration List ---\n")
str(config)


# -------------------------------------------------------------------
# Step 3: Run the EM algorithm
# -------------------------------------------------------------------
cat("\n--- Running EM Algorithm ---\n")
em_results <- tryCatch({
  run_em_conger_cpp(prep_results, config)
}, error = function(e) {
  cat("Error during EM algorithm execution:\n")
  print(e)
  NULL
})


# -------------------------------------------------------------------
# Step 4: Perform Sanity Checks on the Output
# -------------------------------------------------------------------
if (!is.null(em_results)) {
  cat("\n--- EM Algorithm Results ---\n")
  print(em_results)

  cat("\n--- Sanity Checks on Variance Matrix ---\n")
  var_mat <- em_results$var_theta_final
  theta_vec <- em_results$theta_final

  cat("Dimensions of theta_final:", length(theta_vec), "\n")
  cat("Dimensions of var_theta_final:", dim(var_mat), "\n")

  if (length(theta_vec) > 1 && all(dim(var_mat) == length(theta_vec))) {
    # The variance matrix should be singular. Its rank should be n-1.
    # We can check this by looking at the eigenvalues. One should be near zero.
    eigenvalues <- eigen(var_mat, only.values = TRUE)$values
    cat("\nEigenvalues of the variance matrix (one should be ~0):\n")
    print(eigenvalues, digits = 4)

    # The sum of each row/column should be close to zero, due to the sum-to-one constraint
    cat("\nRow sums of the variance matrix (should be ~0):\n")
    print(rowSums(var_mat), digits = 4)

    # The variance of the sum of thetas should be zero.
    # Var(sum(theta)) = 1' * Var(theta) * 1
    var_of_sum <- t(rep(1, length(theta_vec))) %*% var_mat %*% rep(1, length(theta_vec))
    cat("\nVariance of the sum of probabilities (should be ~0):", var_of_sum, "\n")
  } else {
    cat("\nCould not perform variance matrix sanity checks (theta vector length <= 1 or matrix dimensions mismatch).\n")
  }
}

# ... (Previous R script code is the same) ...

# -------------------------------------------------------------------
# Step 5: Calculate Kappa from EM Results
# -------------------------------------------------------------------
if (!is.null(em_results) && length(em_results$theta_final) > 1) {
  cat("\n--- Calculating Unweighted Kappa ---\n")

  # Use the unweighted loss matrix from our config list
  unweighted_loss_matrix <- config$loss_matrices[[1]]

  kappa_results <- calculate_kappa_conger_cpp(
    prep_list = prep_results,
    em_results_list = em_results,
    loss_matrix = unweighted_loss_matrix
  )

  print(kappa_results)

  # As a basic verification, let's calculate kappa on the complete data subset
  # and see if it's in the same ballpark.
  cat("\n--- Comparison with Complete Data Subset ---\n")
}





irr <- conger(x)

kappa_results$kappa
irr[1]

sqrt(kappa_results$kappa_var)
irr[2]
