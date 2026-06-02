# This script runs a master simulation study to test the consistency of the
# available-case (AC) Conger's kappa estimator under four critical missingness
# mechanisms. The goal is to provide empirical validation for the theoretical
# conditions under which the AC estimator is consistent.

# We only need the irrCAC for the analytical variance calculation
library(irrCAC)
library(tidyr)
library(knitr)

conger_ac_use <- \(x) {
  res <- calculate_kappas_binned_cpp(as.matrix(x), use_ipw = TRUE)
  c(kappa = as.numeric(res$kappa_estimates[2]), se = as.numeric((sqrt(diag(res$kappa_cov_matrix)[2]))))
}

#conger_ac_use <- conger_ac_kappa

construct_fisher_ci <- function(kappa, se, conf.level = 0.95) {
  if (is.na(kappa) || is.na(se) || se <= 1e-9 || abs(kappa) >= 1) return(c(NA, NA))
  z_crit <- qnorm(1 - (1 - conf.level) / 2)
  z_val <- 0.5 * log((1 + kappa) / (1 - kappa))
  se_z <- se / (1 - kappa^2)
  ci_z <- z_val + c(-1, 1) * z_crit * se_z
  ci_kappa <- (exp(2 * ci_z) - 1) / (exp(2 * ci_z) + 1)
  return(ci_kappa)
}

# --- MISSINGNESS FUNCTIONS ---
introduce_missingness <- function(data, mechanism_def) {
  n_subjects <- nrow(data)
  n_raters <- ncol(data)
  if (is.array(mechanism_def) && length(dim(mechanism_def)) == n_raters) {
    if (any(dim(mechanism_def) != 2)) stop("Tensor dimensions must be 2.")
    if (abs(sum(mechanism_def) - 1) > 1e-9) stop("Tensor probabilities must sum to 1.")
    all_patterns <- expand.grid(replicate(n_raters, c(1, 0), simplify = FALSE))
    sampled_indices <- sample.int(2^n_raters, size = n_subjects, replace = TRUE, prob = c(mechanism_def))
    chosen_patterns <- all_patterns[sampled_indices, , drop = FALSE]
    na_mask <- (chosen_patterns == 0)
    data[na_mask] <- NA
  } else if (is.vector(mechanism_def) && !is.list(mechanism_def)) {
    for (j in 1:n_raters) {
      na_mask <- sample(c(TRUE, FALSE), size = n_subjects, replace = TRUE, prob = c(mechanism_def[j], 1 - mechanism_def[j]))
      data[na_mask, j] <- NA
    }
  } else {
    stop("Invalid 'mechanism_def' format.")
  }
  return(data)
}
# --- SIMULATION TRIAL (NOW RETURNS MORE DATA) ---
run_simulation_trial <- function(n_subjects, skill_vec, mechanism_def, kappa_population, seed) {
  set.seed(seed)
  data_complete <- simulate_jsm(n = n_subjects, s = skill_vec, model = "fleiss", true_dist = c(0.5, 0.3, 0.2))

  # Introduce missingness
  if (is.function(mechanism_def)) {
    data_missing <- mechanism_def(data_complete)
  } else {
    data_missing <- introduce_missingness(data_complete, mechanism_def)
  }

  # Get the AC estimate and its analytical SE
  fit_ac <- tryCatch(conger_ac_use(data_missing), error = function(e) c(kappa=NA, se=NA))
  kappa_est <- fit_ac["kappa"]
  se_est <- fit_ac["se"]

  # Calculate CI coverage
  coverage <- NA
  if (all(!is.na(fit_ac))) {
    ci_ac <- construct_fisher_ci(kappa_est, se_est)
    if (all(!is.na(ci_ac))) {
      coverage <- kappa_population >= ci_ac[1] && kappa_population <= ci_ac[2]
    }
  }

  # Return everything we need
  return(data.frame(
    estimate = kappa_est,
    se_analytic = se_est,
    coverage = coverage
  ))
}


# --- MASTER SIMULATION ENGINE (NOW WITH FULL ANALYSIS) ---
run_master_performance_study <- function(reps = 1000) {
  n_raters <- 3

  # --- Mechanism Definitions (same as before) ---
  # 1. Uniform PMCAR Tensor
  c1 <- 0.8; c12 <- 0.7; c123 <- 0.6
  p_111 <- c123; p_110 <- c12 - c123; p_101 <- c12 - c123; p_011 <- c12 - c123
  p_100 <- c1 - p_110 - p_101 - p_111; p_010 <- c1 - p_110 - p_011 - p_111; p_001 <- c1 - p_101 - p_011 - p_111
  p_000 <- 1 - (p_111+p_110+p_101+p_011+p_100+p_010+p_001)
  uniform_pmcar_tensor <- array(0, dim = rep(2, n_raters))
  uniform_pmcar_tensor[1,1,1] <- p_111; uniform_pmcar_tensor[1,1,2] <- p_110; uniform_pmcar_tensor[1,2,1] <- p_101; uniform_pmcar_tensor[2,1,1] <- p_011
  uniform_pmcar_tensor[1,2,2] <- p_100; uniform_pmcar_tensor[2,1,2] <- p_010; uniform_pmcar_tensor[2,2,1] <- p_001; uniform_pmcar_tensor[2,2,2] <- p_000

  # 2. Non-Uniform PMCAR Tensor
  non_uniform_pmcar_tensor <- array(0, dim = rep(2, n_raters))
  non_uniform_pmcar_tensor[1,1,1] <- 0.40; non_uniform_pmcar_tensor[1,1,2] <- 0.10; non_uniform_pmcar_tensor[1,2,1] <- 0.15; non_uniform_pmcar_tensor[2,1,1] <- 0.05
  non_uniform_pmcar_tensor[1,2,2] <- 0.10; non_uniform_pmcar_tensor[2,1,2] <- 0.10; non_uniform_pmcar_tensor[2,2,1] <- 0.05; non_uniform_pmcar_tensor[2,2,2] <- 0.05

  mar_fn <- function(data) {
    n_subjects <- nrow(data); data_missing <- data
    data_missing[sample(n_subjects, floor(0.2*n_subjects)), 2] <- NA
    for (i in 1:n_subjects) {
      if (data[i, 1] == 1 && runif(1) < 0.5) data_missing[i, 3] <- NA
      else if (data[i, 1] != 1 && runif(1) < 0.05) data_missing[i, 3] <- NA
    }
    return(data_missing)
  }

  # --- Case Definitions (same as before) ---
  case_definitions <- list(
    #"1_UnifPMCAR_Het"    = list(skills = c(1, 0.5, 0.6), missing = uniform_pmcar_tensor),
    #"2_NonUnifPMCAR_Exch" = list(skills = rep(0.8, 3),      missing = non_uniform_pmcar_tensor),
    "3_NonUnifPMCAR_Het"  = list(skills = c(1, 0.5, 0.6), missing = non_uniform_pmcar_tensor)
    #"4_MAR_Exch"         = list(skills = rep(0.8, 3),      missing = mar_fn),
    #"5_MAR_Unexch"         = list(skills = c(1, 0.5, 0.6),      missing = mar_fn)
  )

  # --- Simulation Grid and Loop ---
  n_subjects <- c(50, 100, 500, 1000, 2000)
  n_subjects <- c(2000, 5000)
  param_grid <- tidyr::crossing(case_name = names(case_definitions), n_subjects = n_subjects)
  cat("Running Full Performance Study...\n")
  full_results <- list()

  for (i in 1:nrow(param_grid)) {
    params <- param_grid[i, ]; case_params <- case_definitions[[params$case_name]]
    kappa_pop <- attr(simulate_jsm(n=1, s=case_params$skills, model="fleiss", true_dist=c(0.5,0.3,0.2)), "kappa")

    cat(sprintf("  - Running: %s, n = %d\n", params$case_name, params$n_subjects))
    reps_list <- lapply(1:reps, function(rep_seed) {
      run_simulation_trial(n_subjects = params$n_subjects, skill_vec = case_params$skills,
                           mechanism_def = case_params$missing, kappa_population = kappa_pop, seed = (i * reps) + rep_seed)
    })

    # Aggregate all results for this condition
    results_df <- do.call(rbind, reps_list)
    results_df <- na.omit(results_df) # Remove failed runs

    # Calculate all performance metrics
    if (nrow(results_df) > 1) {
      bias <- mean(results_df$estimate) - kappa_pop
      variance <- var(results_df$estimate)
      mse <- bias^2 + variance
      se_empirical <- sd(results_df$estimate)
      se_analytic_mean <- mean(results_df$se_analytic)
      coverage <- mean(results_df$coverage)
      n_valid_reps <- nrow(results_df)
    } else {
      # Handle cases where all reps fail
      bias <- variance <- mse <- se_empirical <- se_analytic_mean <- coverage <- NA
      n_valid_reps <- 0
    }

    # Store the summarized results
    summary_row <- data.frame(
      case_name = params$case_name,
      n_subjects = params$n_subjects,
      true_kappa = kappa_pop,
      bias = bias,
      variance = variance,
      mse = mse,
      se_empirical = se_empirical,
      se_analytic = se_analytic_mean,
      coverage = coverage,
      n_valid = n_valid_reps
    )
    full_results[[i]] <- summary_row
  }

  final_summary_df <- do.call(rbind, full_results)
  return(final_summary_df)
}

# --- Format and Print Results ---
format_full_results_table <- function(summary_df) {
  # Select and format columns for a clear table
  table_data <- summary_df[, c("case_name", "n_subjects", "bias", "se_empirical", "se_analytic", "coverage")]
  table_data$bias <- sprintf("%.4f", table_data$bias)
  table_data$se_empirical <- sprintf("%.4f", table_data$se_empirical)
  table_data$se_analytic <- sprintf("%.4f", table_data$se_analytic)
  table_data$coverage <- sprintf("%.3f", table_data$coverage)

  # You can create multiple tables, e.g., one for bias/SE and one for coverage
  kable(table_data, format = "pipe", align = "lccccc",
        caption = "Full Performance Metrics for the Available-Case Estimator")
}

# --- Run the study and print the results table ---
set.seed(456)
full_performance_results <- run_master_performance_study(reps = 1000)
print(format_full_results_table(full_performance_results))
