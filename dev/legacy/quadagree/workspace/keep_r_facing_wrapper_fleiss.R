# Load Rcpp to compile and load the C++ functions
library(Rcpp)

# =============================================================================
# R USER-FACING WRAPPER
# =============================================================================

#' Calculate Fleiss' Kappa with Missing Data
#'
#' A flexible and powerful function to calculate Fleiss' Kappa and its
#' confidence intervals from data with missing ratings. It supports various
#' weighting schemes and bootstrap methods.
#'
#' @param data An N x c matrix or data.frame of counts, where N is the number
#'   of subjects and c is the number of categories.
#' @param r_model The total number of raters in the study. If NULL, it defaults
#'   to the maximum number of ratings observed for any single subject.
#' @param weights A character vector specifying weighting schemes (e.g.,
#'   "unweighted", "linear", "quadratic") or a single custom weight matrix.
#' @param category_values A numeric vector of values for each category,
#'   required for linear and quadratic weights. Defaults to `1:c`.
#' @param conf_level The confidence level for intervals.
#' @param transform A variance-stabilizing transformation for the asymptotic
#'   confidence interval. One of "none", "fisher", "log", "arcsin".
#' @param bootstrap The type of bootstrap to perform. One of "none",
#'   "nonparametric", or "parametric".
#' @param bootstrap_reps The number of bootstrap replications.
#' @param seed An optional integer to seed the random number generator for
#'   reproducible bootstrap results.
#' @return An object of class `missfleiss` containing the results and diagnostics.
fleiss_kappa_missing <- function(
    data,
    r_model = NULL,
    weights = "unweighted",
    category_values = NULL,
    conf_level = 0.95,
    transform = "none",
    bootstrap = "none",
    bootstrap_reps = 1000,
    seed = NULL
) {

  # --- 1. Input Validation and Preparation ---
  if (!is.matrix(data) && !is.data.frame(data)) stop("`data` must be a matrix or data.frame.")
  X <- as.matrix(data)
  if (!is.numeric(X)) stop("`data` must be numeric.")
  storage.mode(X) <- "integer"

  if (is.null(r_model)) {
    r_model <- max(rowSums(X, na.rm = TRUE))
    message(paste("`r_model` not specified. Setting to max observed raters:", r_model))
  }

  c <- ncol(X)
  if (is.null(category_values)) {
    category_values <- 1:c
  }

  # --- 2. Build Loss Matrices ---
  weight_keys <- character()
  loss_matrices <- list()

  if (is.matrix(weights)) {
    if(nrow(weights) != c || ncol(weights) !=c) stop("Custom weight matrix must be c x c.")
    weight_keys <- c("custom")
    loss_matrices <- list(weights)
  } else {
    valid_weights <- c("unweighted", "linear", "quadratic")
    if (!all(weights %in% valid_weights)) stop("Invalid standard weight specified.")
    for (w in weights) {
      weight_keys <- c(weight_keys, w)
      loss_matrices <- c(loss_matrices, list(generate_loss_matrix_rcpp(w, c, category_values)))
    }
  }

  # --- 3. Build Configuration List for C++ ---
  config <- list(
    X = X,
    r_model = as.integer(r_model),
    c = as.integer(c),
    weight_keys = weight_keys,
    loss_matrices = loss_matrices,
    conf_level = conf_level,
    transform_type = tolower(transform),
    bootstrap_method = tolower(bootstrap),
    bootstrap_reps = as.integer(bootstrap_reps),
    seed = if (!is.null(seed)) as.integer(seed) else NULL
  )

  # --- 4. Call the Rcpp function ---
  result_list <- run_analysis_rcpp(config)

  # --- 5. Format and Return Output ---
  class(result_list) <- "missfleiss"
  return(result_list)
}

# Custom print method for our result object
print.missfleiss <- function(x, ...) {
  cat("Fleiss' Kappa with Missing Data\n\n")
  cat("Agreement Estimates:\n")
  print(x$results, row.names = FALSE, digits=3)
  cat("\nEM Diagnostics:\n")
  if(!is.null(x$diagnostics$iterations)) {
    cat(paste("  Converged:", x$diagnostics$converged, "in", x$diagnostics$iterations, "iterations\n"))
  } else {
    cat("  (Diagnostics not available for this run)\n")
  }
}


# =============================================================================
# TEST SUITE
# =============================================================================
cat("===========================================\n")
cat("          MISSFLEISS TEST SUITE          \n")
cat("===========================================\n\n")

# Create a sample dataset `x_test` with missing data
# 30 subjects, 5 categories. Total raters in study is 6.
set.seed(123)
# Base probabilities for categories
base_probs <- c(0.4, 0.3, 0.15, 0.1, 0.05)
# Generate complete data first
x_complete <- t(rmultinom(30, size = 6, prob = base_probs))
# Create the test dataset with missingness
x_test <- x_complete
# Subjects 1-5 have 5 raters
x_test[1:5, ] <- t(rmultinom(5, size = 5, prob = base_probs))
# Subjects 6-10 have 4 raters
x_test[6:10, ] <- t(rmultinom(5, size = 4, prob = base_probs))


cat("--- TEST 1: Default run (unweighted kappa) ---\n")
test1 <- fleiss_kappa_missing(x_test, r_model = 6)
print(test1)
cat("\n\n")


cat("--- TEST 2: Multiple weights (unweighted and quadratic) ---\n")
test2 <- fleiss_kappa_missing(x_test, r_model = 6, weights = c("unweighted", "quadratic"))
print(test2)
cat("\n\n")


cat("--- TEST 3: Fisher transform for CIs ---\n")
test3 <- fleiss_kappa_missing(x_test, r_model = 6, transform = "fisher")
print(test3)
cat("\n\n")


cat("--- TEST 4: Non-parametric bootstrap ---\n")
# NOTE: This is slow. Using fewer reps for the test.
test4 <- fleiss_kappa_missing(x_test, r_model = 6, bootstrap = "nonparametric", bootstrap_reps = 200, seed = 42)
print(test4)
cat("Structure of bootstrap replicates:\n")
str(test4$diagnostics$bootstrap_replicates)
cat("\n\n")

cat("--- TEST 5: Parametric bootstrap ---\n")
test5 <- fleiss_kappa_missing(x_test, r_model = 6, bootstrap = "parametric", bootstrap_reps = 1000, seed = 42)
print(test5)
cat("\n\n")

cat("--- TEST 6: Theoretical shrinkage case (r_model > r_obs) ---\n")
# Using the complete data for a clean demonstration
kappa_at_6 <- fleiss_kappa_missing(x_complete, r_model = 6)$results$kappa
test6 <- fleiss_kappa_missing(x_complete, r_model = 8)
print(test6)
cat("Theoretical Shrinkage Factor:", (6*5)/(8*7), "\n")
cat("Theoretical Kappa at r=8:", kappa_at_6 * (6*5)/(8*7), "\n")
cat("\n\n")

cat("--- TEST 7: Custom weight matrix ---\n")
# A weight matrix that heavily penalizes disagreement with category 1
custom_w <- matrix(1, nrow=5, ncol=5)
diag(custom_w) <- 0
custom_w[1, ] <- 5
custom_w[, 1] <- 5
custom_w[1,1] <- 0
cat("Custom Weight Matrix:\n")
print(custom_w)
test7 <- fleiss_kappa_missing(x_test, r_model = 6, weights = custom_w)
print(test7)
cat("\n\n")

cat("--- TEST 8: Custom weight matrix WITH bootstrap ---\n")
# This is the true test to ensure custom weights work with the bootstrap engine.
# Using the same custom matrix as before.
test8 <- fleiss_kappa_missing(
  x_test,
  r_model = 6,
  weights = custom_w,
  bootstrap = "nonparametric",
  bootstrap_reps = 200,
  seed = 42
)
print(test8)
cat("\n\n")
