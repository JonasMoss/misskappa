# This script conducts a confidence interval (CI) coverage study to demonstrate
# the practical consequences of using a biased estimator. It compares the CIs
# from misskappa's MLE to those from a naive available-case (AC) estimator.
#
# The core hypothesis is that because the AC estimator is biased, its CIs will
# fail to achieve the nominal coverage rate (e.g., 95%), even if its standard
# error is estimated well. The MLE's CIs, being centered on a consistent
# estimate, should achieve the correct coverage.
#
# The simulation uses the same challenging scenario as the recovery study:
# non-exchangeable raters and non-exchangeable MAR.

library(misskappa)
library(ggplot2)

#' Naive Conger's Kappa for Available-Case Data
#' @keywords internal
conger_kappa_naive <- function(x) {
  r <- ncol(x); categories <- sort(unique(na.omit(as.vector(x)))); q <- length(categories)
  total_pairs <- 0; total_agreements <- 0
  for (j1 in 1:(r - 1)) for (j2 in (j1 + 1):r) {
    complete_cases <- !is.na(x[, j1]) & !is.na(x[, j2]); n_complete <- sum(complete_cases)
    if (n_complete > 0) {
      agreements <- sum(x[complete_cases, j1] == x[complete_cases, j2])
      total_pairs <- total_pairs + n_complete; total_agreements <- total_agreements + agreements
    }
  }
  pa <- if (total_pairs > 0) total_agreements / total_pairs else 0
  p_j <- matrix(0, nrow = r, ncol = q)
  for (j in 1:r) {
    rater_obs <- na.omit(x[, j])
    if (length(rater_obs) > 0) p_j[j, ] <- table(factor(rater_obs, levels = categories)) / length(rater_obs)
  }
  pe_sum <- sum(sapply(1:(r - 1), function(j1) sum(sapply((j1 + 1):r, function(j2) sum(p_j[j1, ] * p_j[j2, ])))))
  pe <- pe_sum / (r * (r - 1) / 2)
  return(if (abs(1 - pe) < 1e-9) 0 else (pa - pe) / (1 - pe))
}

#' Introduce MAR by Rater
#' @keywords internal
introduce_mar_by_rater <- function(data, props_by_rater) {
  data_missing <- data
  for (j in seq_len(ncol(data))) {
    n_to_na <- floor(nrow(data) * props_by_rater[j])
    if (n_to_na > 0) {
      na_indices <- sample(nrow(data), size = n_to_na, replace = FALSE)
      data_missing[na_indices, j] <- NA
    }
  }
  return(data_missing)
}

#' Run a Single Simulation Trial for CI Coverage
#' @return A data frame with coverage indicators (TRUE/FALSE) for both estimators.
#' @keywords internal
run_ci_trial <- function(n_subjects, skill_vec, missing_props, kappa_population, seed) {
  set.seed(seed)
  z_crit <- qnorm(0.975) # For 95% CI

  # 1. Generate data
  data_complete <- simulate_jsm(
    n = n_subjects, s = skill_vec, model = "fleiss", true_dist = c(0.5, 0.3, 0.2)
  )
  data_missing <- introduce_mar_by_rater(data_complete, missing_props)

  # 2. MLE Confidence Interval
  fit_mle <- tryCatch(conger_kappa(data_missing), error = function(e) NULL)
  coverage_mle <- NA
  if (!is.null(fit_mle)) {
    ci_mle <- fit_mle$estimate + c(-1, 1) * z_crit * fit_mle$se
    coverage_mle <- kappa_population >= ci_mle[1] && kappa_population <= ci_mle[2]
  }

  # 3. Available-Case Confidence Interval
  # We must bootstrap to get a standard error for the naive estimator.
  # This gives it the "best possible" chance to have correct coverage.
  est_ac <- tryCatch(conger_kappa_naive(data_missing), error = function(e) NA)
  coverage_ac <- NA
  if (!is.na(est_ac)) {
    boot_reps_ac <- replicate(200, {
      idx <- sample(1:n_subjects, size = n_subjects, replace = TRUE)
      conger_kappa_naive(data_missing[idx, ])
    })
    se_ac <- sd(boot_reps_ac, na.rm = TRUE)
    ci_ac <- est_ac + c(-1, 1) * z_crit * se_ac
    coverage_ac <- kappa_population >= ci_ac[1] && kappa_population <= ci_ac[2]
  }

  data.frame(
    estimator = c("MLE", "Available-Case"),
    coverage = c(coverage_mle, coverage_ac)
  )
}

#' Run the Full CI Coverage Study
#' @param reps Number of Monte Carlo repetitions per setting.
#' @return A data frame summarizing CI coverage for each estimator and sample size.
#' @keywords internal
run_ci_coverage_study <- function(reps = 1000) {
  skill_vec <- c(0.9, 0.9, 0.5, 0.5)
  missing_props <- c(0.1, 0.1, 0.4, 0.4)
  kappa_population <- attr(simulate_jsm(
    n = 1, s = skill_vec, model = "fleiss", true_dist = c(0.5, 0.3, 0.2)
  ), "kappa")

  param_grid <- expand.grid(n_subjects = c(50, 100, 500))

  cat("Running CI coverage study...\n")

  full_results <- list()
  for (i in 1:nrow(param_grid)) {
    n_sub <- param_grid$n_subjects[i]
    cat(sprintf("  - n = %d\n", n_sub))

    reps_list <- lapply(1:reps, function(rep_seed) {
      run_ci_trial(n_subjects = n_sub, skill_vec = skill_vec,
                   missing_props = missing_props,
                   kappa_population = kappa_population, seed = rep_seed * i)
    })

    results_df <- do.call(rbind, reps_list)
    results_df$n_subjects <- n_sub
    full_results[[i]] <- results_df
  }

  final_df <- do.call(rbind, full_results)

  # Calculate coverage probability
  aggregate(coverage ~ estimator + n_subjects, data = final_df, FUN = mean, na.rm = TRUE)
}

#' Plot the Results of the CI Coverage Study
#' @return A ggplot object.
#' @keywords internal
plot_ci_coverage_results <- function(summary_df) {

  summary_df$estimator <- factor(summary_df$estimator, levels = c("MLE", "Available-Case"))

  ggplot(summary_df, aes(x = n_subjects, y = coverage, color = estimator, shape = estimator)) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "gray20") +
    geom_line(linewidth = 1) +
    geom_point(size = 4, fill = "white") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
    scale_x_continuous(breaks = summary_df$n_subjects) +
    scale_color_manual(values = c("MLE" = "black", "Available-Case" = "gray50")) +
    scale_shape_manual(values = c("MLE" = 16, "Available-Case" = 17)) +
    labs(
      title = "Confidence Interval Coverage",
      subtitle = "Scenario: Non-exchangeable raters and MAR",
      x = "Sample Size (n)",
      y = "Empirical Coverage Probability",
      color = "Estimator",
      shape = "Estimator"
    ) +
    theme_bw(base_size = 14) +
    theme(legend.position = "bottom")
}

# --- Run the study and plot the results ---
coverage_results <- run_ci_coverage_study(reps = 1000)
print(coverage_results)
plot_ci_coverage_results(coverage_results)
