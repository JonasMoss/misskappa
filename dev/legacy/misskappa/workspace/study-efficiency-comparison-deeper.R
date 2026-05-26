# This script is a "deep dive" into the efficiency comparison between the
# misskappa MLE and the naive Available-Case (AC) estimator. It is designed to
# create a scenario that magnifies the efficiency gain of the MLE.
#
# The core hypothesis is that the MLE's efficiency advantage, while minimal for
# perfectly exchangeable raters, becomes significant and practically important
# as soon as rater heterogeneity is introduced.
#
# The simulation compares two conditions, both under Exchangeable MCAR to ensure
# both estimators are consistent:
# 1. Moderate Rater Heterogeneity: A baseline with a linear spread of skills.
# 2. High Rater Heterogeneity: A challenging scenario with "expert" and
#    "novice" raters, designed to highlight the MLE's superior use of information.

library(misskappa)
library(ggplot2)
library(tidyr)
library(patchwork)

#' Naive Conger's Kappa for Available-Case Data (from bias study)
#' @keywords internal
conger_kappa_naive <- function(x) {
  r <- ncol(x)
  categories <- sort(unique(na.omit(as.vector(x))))
  q <- length(categories)
  total_pairs <- 0; total_agreements <- 0
  for (j1 in 1:(r - 1)) {
    for (j2 in (j1 + 1):r) {
      complete_cases <- !is.na(x[, j1]) & !is.na(x[, j2])
      n_complete <- sum(complete_cases)
      if (n_complete > 0) {
        agreements <- sum(x[complete_cases, j1] == x[complete_cases, j2])
        total_pairs <- total_pairs + n_complete
        total_agreements <- total_agreements + agreements
      }
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

#' Run a Single Simulation Trial for the Deep Dive Study
#' @keywords internal
run_deep_dive_trial <- function(n_subjects, missing_prop, skill_vec, seed) {
  set.seed(seed)
  sim_data_complete <- simulate_jsm(
    n = n_subjects, s = skill_vec, model = "fleiss", true_dist = c(0.6, 0.3, 0.1)
  )
  true_kappa <- attr(sim_data_complete, "kappa")
  sim_data_missing <- sim_data_complete
  n_to_na <- floor(length(sim_data_missing) * missing_prop)
  sim_data_missing[sample(length(sim_data_missing), n_to_na)] <- NA
  est_misskappa <- tryCatch(conger_kappa(sim_data_missing)$estimate, error = function(e) NA)
  est_naive <- tryCatch(conger_kappa_naive(sim_data_missing), error = function(e) NA)
  data.frame(
    true_kappa = true_kappa, est_misskappa = est_misskappa, est_naive = est_naive
  )
}

#' Run the Full Efficiency Deep Dive Study
#' @keywords internal
run_efficiency_deep_dive <- function(reps = 1500) {
  # Define the two scenarios for rater skills
  skill_scenarios <- list(
    "Moderate Heterogeneity" = c(0.9, 0.8, 0.7, 0.6),
    "High Heterogeneity (Experts/Novices)" = c(0.9, 0.7, 0.5, 0.3)
  )

  param_grid <- expand.grid(
    n_subjects = c(50, 100, 200, 400, 800),
    missing_prop = 0.35,
    skill_scenario = names(skill_scenarios),
    stringsAsFactors = FALSE
  )

  full_results <- list()
  cat("Running efficiency deep dive simulation...\n")
  for (i in 1:nrow(param_grid)) {
    params <- param_grid[i, ]
    skill_vec <- skill_scenarios[[params$skill_scenario]]
    cat(sprintf("  - Scenario: %s, n=%-4d\n", params$skill_scenario, params$n_subjects))

    reps_list <- lapply(1:reps, function(rep_seed) {
      trial_res <- run_deep_dive_trial(
        n_subjects = params$n_subjects, missing_prop = params$missing_prop,
        skill_vec = skill_vec, seed = rep_seed * i
      )
      # Add parameters to result
      trial_res$n_subjects <- params$n_subjects
      trial_res$skill_scenario <- params$skill_scenario
      trial_res
    })
    full_results[[i]] <- do.call(rbind, reps_list)
  }

  results_df <- do.call(rbind, full_results)

  summary_list <- by(results_df, results_df[, c("n_subjects", "skill_scenario")], function(group) {
    n <- group$n_subjects[1]
    true_k <- mean(group$true_kappa, na.rm=TRUE)
    calc_metrics <- function(est) {
      error <- est - true_k
      c(bias = mean(error, na.rm = TRUE),
        var = var(est, na.rm = TRUE),
        mse = mean(error^2, na.rm = TRUE))
    }
    metrics_mle <- calc_metrics(group$est_misskappa)
    metrics_ac <- calc_metrics(group$est_naive)
    data.frame(
      n_subjects = n, skill_scenario = group$skill_scenario[1],
      estimator = c("misskappa (MLE)", "Naive (AC)"),
      bias = c(metrics_mle["bias"], metrics_ac["bias"]),
      rescaled_var = c(n * metrics_mle["var"], n * metrics_ac["var"]),
      rescaled_mse = c(n * metrics_mle["mse"], n * metrics_ac["mse"])
    )
  })

  do.call(rbind, summary_list)
}

#' Plot the Results of the Deep Dive Study
#' @keywords internal
plot_deep_dive_results <- function(summary_df) {

  summary_df$skill_scenario <- factor(
    summary_df$skill_scenario,
    levels = c("Moderate Heterogeneity", "High Heterogeneity (Experts/Novices)")
  )

  # Plot 1: Bias
  p_bias <- ggplot(summary_df, aes(x = n_subjects, y = bias, color = estimator, shape = estimator)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    geom_line() + geom_point(size = 2.5) +
    facet_wrap(~ skill_scenario) +
    labs(subtitle = "A) Estimator Bias", y = "Bias") +
    theme(legend.position = "none")

  # Plot 2: Asymptotic Variance
  p_var <- ggplot(summary_df, aes(x = n_subjects, y = rescaled_var, color = estimator, shape = estimator)) +
    geom_line() + geom_point(size = 2.5) +
    facet_wrap(~ skill_scenario, scales = "free_y") +
    labs(subtitle = "B) Asymptotic Variance", y = expression(n %*% "Variance")) +
    theme(legend.position = "none")

  # Plot 3: Asymptotic MSE
  p_mse <- ggplot(summary_df, aes(x = n_subjects, y = rescaled_mse, color = estimator, shape = estimator)) +
    geom_line() + geom_point(size = 2.5) +
    facet_wrap(~ skill_scenario, scales = "free_y") +
    labs(subtitle = "C) Asymptotic MSE", x = "Sample Size (n)") +
    theme(legend.position = "bottom")

  (p_bias / p_var / p_mse) +
    plot_annotation(
      title = "Efficiency Deep Dive: Comparing Consistent Estimators",
      subtitle = "Conditions: Exchangeable MCAR with 35% missingness",
      theme = theme(plot.title = element_text(hjust = 0.5, size = 16),
                    plot.subtitle = element_text(hjust = 0.5, size = 12))
    ) &
    theme_bw(base_size = 11) &
    labs(color = "Estimator", shape = "Estimator") &
    scale_x_continuous(breaks = unique(summary_df$n_subjects)) &
    expand_limits(y = 0)
}


# --- Run the study and plot the results ---
# Note: Using a higher 'reps' value (e.g., 2000+) is recommended for final paper figures.
efficiency_results <- run_efficiency_deep_dive(reps = 2000)
print(efficiency_results)
plot_deep_dive_results(efficiency_results)
