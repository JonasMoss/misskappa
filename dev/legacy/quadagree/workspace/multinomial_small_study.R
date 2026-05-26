# =============================================================================
#              Simulation Study for Fleiss' Kappa with Missing Data
# =============================================================================

library(Rcpp)
library(knitr)
library(dplyr)
library(tidyr)

# The user-facing R wrapper function (from our previous work)
fleiss_kappa_missing <- function(
    data, r_model = NULL, weights = "unweighted", category_values = NULL,
    conf_level = 0.95, transform = "none", bootstrap = "none",
    bootstrap_reps = 1000, seed = NULL) {
  # ... (full implementation from previous response) ...
  if (!is.matrix(data) && !is.data.frame(data)) stop("`data` must be a matrix or data.frame.")
  X <- as.matrix(data)
  if (!is.numeric(X)) stop("`data` must be numeric.")
  storage.mode(X) <- "integer"
  if (is.null(r_model)) {
    r_model <- max(rowSums(X, na.rm = TRUE))
    message(paste("`r_model` not specified. Setting to max observed raters:", r_model))
  }
  c <- ncol(X)
  if (is.null(category_values)) { category_values <- 1:c }
  weight_keys <- character(); loss_matrices <- list()
  if (is.matrix(weights)) {
    if(nrow(weights) != c || ncol(weights) !=c) stop("Custom weight matrix must be c x c.")
    weight_keys <- c("custom"); loss_matrices <- list(weights)
  } else {
    valid_weights <- c("unweighted", "linear", "quadratic")
    if (!all(weights %in% valid_weights)) stop("Invalid standard weight specified.")
    for (w in weights) {
      weight_keys <- c(weight_keys, w)
      loss_matrices <- c(loss_matrices, list(generate_loss_matrix_rcpp(w, c, category_values)))
    }
  }
  valid_boot <- c("none", "nonparametric", "parametric", "nonparametric-t", "parametric-t")
  if (!tolower(bootstrap) %in% valid_boot) stop("Invalid bootstrap method specified.")
  config <- list(
    X = X, r_model = as.integer(r_model), c = as.integer(c), weight_keys = weight_keys, loss_matrices = loss_matrices,
    conf_level = conf_level, transform_type = tolower(transform), bootstrap_method = tolower(bootstrap),
    bootstrap_reps = as.integer(bootstrap_reps), seed = if (!is.null(seed)) as.integer(seed) else NULL
  )
  result_list <- run_analysis_rcpp(config)
  class(result_list) <- "missfleiss"
  return(result_list)
}


# --- 2. Simulation Parameters ---

R_MODEL <- 6       # Total number of raters in the study
C_CATS  <- 5       # Number of categories
N_SUBJ  <- 100     # Number of subjects per simulated dataset
N_REPS  <- 1000    # Number of replications per scenario (use a smaller number like 100 for quick tests)
N_BOOT  <- 500     # Number of bootstrap samples within each replication

# --- 3. Define the "True Worlds" (Theta Vectors) ---

# We need the full list of possible patterns for R=6, C=5
# This requires a function to generate all compositions
generate_all_compositions <- function(n, k) {
  if (k == 1) return(matrix(n, 1, 1))
  do.call(rbind, lapply(0:n, function(i) {
    cbind(i, generate_all_compositions(n - i, k - 1))
  }))
}

all_patterns <- generate_all_compositions(R_MODEL, C_CATS)
n_patterns <- nrow(all_patterns)

# Scenario A: High Agreement World
theta_high_agreement <- rep(0, n_patterns)
# Put most probability mass on patterns of perfect agreement
perfect_agreement_indices <- which(apply(all_patterns, 1, function(p) any(p == R_MODEL)))
theta_high_agreement[perfect_agreement_indices] <- 1 / length(perfect_agreement_indices)

# Scenario B: Low Agreement World
theta_low_agreement <- rep(0, n_patterns)
# Put most mass on patterns of high disagreement (e.g., ratings spread out)
high_disagreement_indices <- which(apply(all_patterns, 1, function(p) max(p) <= 2))
theta_low_agreement[high_disagreement_indices] <- 1 / length(high_disagreement_indices)
theta_low_agreement <- theta_low_agreement / sum(theta_low_agreement) # Re-normalize


# --- 4. Calculate the True Kappa for Each World ---

# To do this, we generate one massive, perfect dataset for each theta
# and calculate its kappa. This gives us our ground truth.
message("Calculating true kappa values for simulation worlds...")
truth_data_high <- all_patterns[sample(1:n_patterns, 10000, replace = TRUE, prob = theta_high_agreement), ]
truth_data_low <- all_patterns[sample(1:n_patterns, 10000, replace = TRUE, prob = theta_low_agreement), ]

kappa_true_high <- fleiss_kappa_missing(truth_data_high, r_model = R_MODEL)$results$kappa
kappa_true_low <- fleiss_kappa_missing(truth_data_low, r_model = R_MODEL)$results$kappa

scenarios <- list(
  high_agreement = list(theta = theta_high_agreement, kappa_true = kappa_true_high),
  low_agreement  = list(theta = theta_low_agreement,  kappa_true = kappa_true_low)
)

# --- 5. Helper Function to Introduce Missingness ---

degrade_data <- function(complete_data, n_missing_subjects, n_raters_min, n_raters_max) {
  degraded_data <- complete_data
  if (n_missing_subjects == 0) return(degraded_data)

  # Select subjects to have missing data
  missing_indices <- sample(1:nrow(complete_data), n_missing_subjects)

  for (i in missing_indices) {
    # Get the vector of individual ratings, e.g., (3,2,1,0,0) -> c(1,1,1,2,2,3)
    original_pattern <- complete_data[i, ]
    ratings_long <- rep(1:length(original_pattern), original_pattern)

    # Decide how many raters this subject will have
    n_raters_new <- sample(n_raters_min:n_raters_max, 1)

    # Sample a subset of ratings and re-tabulate
    ratings_new_long <- sample(ratings_long, n_raters_new)
    degraded_data[i, ] <- as.integer(tabulate(ratings_new_long, nbins = length(original_pattern)))
  }
  return(degraded_data)
}


# --- 6. The Main Simulation Loop ---

# Conditions to test
missingness_conditions <- list(
  complete = list(n_miss = 0),
  moderate = list(n_miss = round(N_SUBJ * 0.25), min_r = 4, max_r = 5),
  high     = list(n_miss = round(N_SUBJ * 0.50), min_r = 3, max_r = 5)
)

# Store results
all_results <- list()

for (scen_name in names(scenarios)) {
  for (miss_name in names(missingness_conditions)) {

    message(paste("\nRunning Scenario:", scen_name, "| Missingness:", miss_name))

    kappa_true <- scenarios[[scen_name]]$kappa_true
    theta_true <- scenarios[[scen_name]]$theta
    miss_cond <- missingness_conditions[[miss_name]]

    # Pre-allocate storage for this scenario's results
    results_collector <- data.frame(
      kappa_est = rep(NA, N_REPS),
      se_asym = rep(NA, N_REPS),
      ci_low_asym = rep(NA, N_REPS),
      ci_high_asym = rep(NA, N_REPS),
      ci_low_boot = rep(NA, N_REPS),
      ci_high_boot = rep(NA, N_REPS)
    )

    pb <- txtProgressBar(min = 0, max = N_REPS, style = 3)

    for (i in 1:N_REPS) {

      # A. Generate complete data from the true theta
      complete_data <- all_patterns[sample(1:n_patterns, N_SUBJ, replace = TRUE, prob = theta_true), ]

      # B. Introduce missingness
      degraded_data <- degrade_data(complete_data, miss_cond$n_miss, miss_cond$min_r, miss_cond$max_r)

      # C. Analyze the data
      analysis_result <- tryCatch({
        fleiss_kappa_missing(
          data = degraded_data,
          r_model = R_MODEL,
          transform = "fisher", # Use a good transform
          bootstrap = "nonparametric",
          bootstrap_reps = N_BOOT,
          seed = i # Use loop index for reproducible seeds
        )
      }, error = function(e) {
        warning(paste("Replication", i, "failed:", e$message))
        return(NULL) # Return NULL on failure
      })

      # D. Store results if successful
      if (!is.null(analysis_result)) {
        res <- analysis_result$results
        results_collector$kappa_est[i] <- res$kappa
        results_collector$se_asym[i] <- res$se
        results_collector$ci_low_asym[i] <- res$ci.lower
        results_collector$ci_high_asym[i] <- res$ci.upper
        results_collector$ci_low_boot[i] <- res$ci.lower.boot
        results_collector$ci_high_boot[i] <- res$ci.upper.boot
      }

      setTxtProgressBar(pb, i)
    }
    close(pb)

    # Store the results for this condition
    all_results[[paste(scen_name, miss_name, sep = "_")]] <- list(
      kappa_true = kappa_true,
      results = results_collector
    )
  }
}

# --- 7. Summarize and Display Results ---

summary_list <- lapply(names(all_results), function(name) {

  data <- all_results[[name]]
  res <- na.omit(data$results) # Remove failed replications
  kappa_true <- data$kappa_true

  # Calculate summary stats
  n_success <- nrow(res)
  bias <- mean(res$kappa_est) - kappa_true
  empirical_se <- sd(res$kappa_est)
  avg_asym_se <- mean(res$se_asym)
  coverage_asym <- mean(res$ci_low_asym <= kappa_true & res$ci_high_asym >= kappa_true)
  coverage_boot <- mean(res$ci_low_boot <= kappa_true & res$ci_high_boot >= kappa_true)

  data.frame(
    Scenario = name,
    N_Success = n_success,
    True_Kappa = round(kappa_true, 3),
    Bias = round(bias, 4),
    Empirical_SE = round(empirical_se, 4),
    Avg_Asym_SE = round(avg_asym_se, 4),
    Coverage_Asym = paste0(round(coverage_asym * 100, 1), "%"),
    Coverage_Boot = paste0(round(coverage_boot * 100, 1), "%")
  )
})

summary_table <- do.call(rbind, summary_list)

cat("\n\n===========================================\n")
cat("          SIMULATION RESULTS             \n")
cat("===========================================\n\n")

print(kable(summary_table, format = "pipe", align = 'l'))
