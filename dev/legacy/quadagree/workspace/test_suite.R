
# =============================================================================
# R USER-FACING WRAPPER
# =============================================================================

#' Calculate Fleiss' Kappa with Missing Data
#'
#' @param bootstrap The type of bootstrap. One of "none", "nonparametric",
#'   "parametric", "nonparametric-t", or "parametric-t".
#' @inheritParams other_params
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
      # *** THE FIX IS HERE ***
      loss_matrices <- c(loss_matrices, list(generate_loss_matrix_rcpp(w, c, category_values)))
    }
  }

  # --- 3. Build Configuration List for C++ ---
  valid_boot <- c("none", "nonparametric", "parametric", "nonparametric-t", "parametric-t")
  if (!tolower(bootstrap) %in% valid_boot) stop("Invalid bootstrap method specified.")

  config <- list(
    X = X, # Key is "X"
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

# ... (print method and test data setup are unchanged) ...
# Custom print method
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

# Create a sample dataset `x_test`
set.seed(123)
base_probs <- c(0.4, 0.3, 0.15, 0.1, 0.05)
x_complete <- t(rmultinom(30, size = 6, prob = base_probs))
x_test <- x_complete
x_test[1:5, ] <- t(rmultinom(5, size = 5, prob = base_probs))
x_test[6:10, ] <- t(rmultinom(5, size = 4, prob = base_probs))

cat("--- TEST 1: Default run ---\n")
test1 <- fleiss_kappa_missing(x_test, r_model = 6)
print(test1)
cat("\n\n")

cat("--- TEST 2: Multiple weights ---\n")
test2 <- fleiss_kappa_missing(x_test, r_model = 6, weights = c("unweighted", "quadratic"))
print(test2)
cat("\n\n")

cat("--- TEST 3: Fisher transform ---\n")
test3 <- fleiss_kappa_missing(x_test, r_model = 6, transform = "fisher")
print(test3)
cat("\n\n")

cat("--- TEST 4: Non-parametric percentile bootstrap ---\n")
test4 <- fleiss_kappa_missing(x_test, r_model = 6, bootstrap = "nonparametric", bootstrap_reps = 200, seed = 42)
print(test4)
cat("\n\n")

cat("--- TEST 5: Parametric percentile bootstrap ---\n")
test5 <- fleiss_kappa_missing(x_test, r_model = 6, bootstrap = "parametric", bootstrap_reps = 200, seed = 42)
print(test5)
cat("\n\n")

cat("--- TEST 6: Non-parametric STUDENTIZED bootstrap ---\n")
# NOTE: This is very slow. Using few reps.
test6 <- fleiss_kappa_missing(x_test, r_model = 6, bootstrap = "nonparametric-t", bootstrap_reps = 50, seed = 42)
print(test6)
cat("\n\n")

cat("--- TEST 7: Parametric STUDENTIZED bootstrap ---\n")
test7 <- fleiss_kappa_missing(x_test, r_model = 6, bootstrap = "parametric-t", bootstrap_reps = 200, seed = 42)
print(test7)
cat("\n\n")

cat("--- TEST 8: Custom weight matrix with bootstrap ---\n")
custom_w <- matrix(1, nrow=5, ncol=5)
diag(custom_w) <- 0
test8 <- fleiss_kappa_missing(x_test, r_model = 6, weights = custom_w, bootstrap = "nonparametric", bootstrap_reps = 200, seed = 42)
print(test8)
cat("\n\n")
