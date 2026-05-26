# This script provides an empirical demonstration of the superior efficiency of
# misskappa's MLE estimator compared to a naive available-case (AC) estimator.
#
# The simulation operates under Exchangeable MCAR conditions, a scenario where
# BOTH the MLE and the AC estimators are statistically consistent. This allows
# for a direct comparison of their bias, variance, and Mean Squared Error (MSE).
#
# Key metrics are rescaled by sample size `n` to directly visualize the
# asymptotic properties of the estimators.

library(misskappa)
library(ggplot2)
library(tidyr)

#' Naive Conger's Kappa for Available-Case Data (from bias study)
#' @keywords internal
conger_kappa_naive <- function(x) {
  r <- ncol(x)
  categories <- sort(unique(na.omit(as.vector(x))))
  q <- length(categories)
  total_pairs <- 0
  total_agreements <- 0
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
    if (length(rater_obs) > 0) {
      p_j[j, ] <- table(factor(rater_obs, levels = categories)) / length(rater_obs)
    }
  }
  pe_sum <- sum(sapply(1:(r - 1), function(j1) {
    sum(sapply((j1 + 1):r, function(j2) sum(p_j[j1, ] * p_j[j2, ])))
  }))
  pe <- pe_sum / (r * (r - 1) / 2)
  return(if (abs(1 - pe) < 1e-9) 0 else (pa - pe) / (1 - pe))
}


#' Run a Single Simulation Trial for Efficiency Comparison
#' @return A single-row data frame with estimates, true kappa, and parameters.
#' @keywords internal
run_efficiency_trial <- function(n_subjects, missing_prop, exchangeable_raters, seed) {
  set.seed(seed)
  n_raters <- 4
  skill_vec <- if (exchangeable_raters) rep(0.8, n_raters) else seq(0.6, 0.9, length.out = n_raters)

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
    n_subjects = n_subjects, missing_prop = missing_prop,
    exchangeable_raters = exchangeable_raters, true_kappa = true_kappa,
    est_misskappa = est_misskappa, est_naive = est_naive
  )
}

#' Run the Full Efficiency Simulation Study
#' @return A data frame summarizing performance metrics for each estimator.
#' @keywords internal
run_efficiency_study <- function(reps = 1000) {
  param_grid <- expand.grid(
    n_subjects = c(100, 250, 500, 1000),
    missing_prop = c(0.15, 0.30),
    exchangeable_raters = c(TRUE, FALSE)
  )

  full_results <- list()
  cat("Running efficiency simulation study...\n")
  for (i in 1:nrow(param_grid)) {
    params <- param_grid[i, ]
    cat(sprintf("  - n=%-4d, missing=%.2f, exchangeable=%-5s\n",
                params$n_subjects, params$missing_prop, params$exchangeable_raters))

    reps_list <- lapply(1:reps, function(rep_seed) {
      run_efficiency_trial(
        n_subjects = params$n_subjects, missing_prop = params$missing_prop,
        exchangeable_raters = params$exchangeable_raters, seed = rep_seed * i
      )
    })
    full_results[[i]] <- do.call(rbind, reps_list)
  }

  results_df <- do.call(rbind, full_results)

  summary_list <- by(results_df, results_df[, c("n_subjects", "missing_prop", "exchangeable_raters")], function(group) {
    n <- group$n_subjects[1]
    true_k <- mean(group$true_kappa)

    calc_metrics <- function(est) {
      error <- est - true_k
      bias <- mean(error, na.rm = TRUE)
      variance <- var(est, na.rm = TRUE)
      mse <- mean(error^2, na.rm = TRUE)
      c(bias = bias, var = variance, mse = mse,
        rescaled_var = n * variance, rescaled_mse = n * mse)
    }

    metrics_mle <- calc_metrics(group$est_misskappa)
    metrics_ac <- calc_metrics(group$est_naive)

    data.frame(
      n_subjects = n, missing_prop = group$missing_prop[1],
      exchangeable_raters = group$exchangeable_raters[1],
      estimator = c("misskappa (MLE)", "Naive (AC)"),
      bias = c(metrics_mle["bias"], metrics_ac["bias"]),
      variance = c(metrics_mle["var"], metrics_ac["var"]),
      mse = c(metrics_mle["mse"], metrics_ac["mse"]),
      rescaled_var = c(metrics_mle["rescaled_var"], metrics_ac["rescaled_var"]),
      rescaled_mse = c(metrics_mle["rescaled_mse"], metrics_ac["rescaled_mse"])
    )
  })

  do.call(rbind, summary_list)
}

#' Plot the Results of the Efficiency Study
#'
#' @param summary_df The summary data frame from `run_efficiency_study`.
#' @return A list of ggplot objects.
#' @keywords internal
plot_efficiency_results <- function(summary_df) {

  summary_df$exchangeable_label <- factor(summary_df$exchangeable_raters,
                                          levels = c(TRUE, FALSE),
                                          labels = c("Raters: Exchangeable", "Raters: Non-Exchangeable"))
  summary_df$missing_label <- paste0("Missingness: ", summary_df$missing_prop * 100, "%")

  # Plot 1: Bias
  p_bias <- ggplot(summary_df, aes(x = n_subjects, y = bias, color = estimator, shape = estimator)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray20") +
    geom_line(linewidth = 0.8) + geom_point(size = 3) +
    facet_grid(exchangeable_label ~ missing_label) +
    labs(subtitle = "A) Estimator Bias", y = "Bias") +
    theme_bw(base_size = 12) + theme(legend.position = "none")

  # Plot 2: Rescaled Variance
  p_var <- ggplot(summary_df, aes(x = n_subjects, y = rescaled_var, color = estimator, shape = estimator)) +
    geom_line(linewidth = 0.8) + geom_point(size = 3) +
    facet_grid(exchangeable_label ~ missing_label, scales = "free_y") +
    labs(subtitle = "B) Asymptotic Variance", y = expression(n %*% "Variance")) +
    theme_bw(base_size = 12) + theme(legend.position = "none")

  # Plot 3: Rescaled MSE
  p_mse <- ggplot(summary_df, aes(x = n_subjects, y = rescaled_mse, color = estimator, shape = estimator)) +
    geom_line(linewidth = 0.8) + geom_point(size = 3) +
    facet_grid(exchangeable_label ~ missing_label, scales = "free_y") +
    labs(subtitle = "C) Asymptotic MSE", y = expression(n %*% "MSE"), x = "Sample Size (n)") +
    theme_bw(base_size = 12) + theme(legend.position = "bottom")

  # Combine plots
  if (requireNamespace("patchwork", quietly = TRUE)) {
    (p_bias / p_var / p_mse) +
      patchwork::plot_layout(guides = 'collect') &
      theme(legend.position = "bottom") &
      labs(x = "Sample Size (n)") &
      patchwork::plot_annotation(
        title = "Performance of Consistent Estimators under Exchangeable MCAR",
        theme = theme(plot.title = element_text(hjust = 0.5, size = 16))
      )
  } else {
    warning("Install 'patchwork' for combined plotting.")
    list(bias = p_bias, variance = p_var, mse = p_mse)
  }
}


# --- Run the study and plot the results ---
efficiency_summary <- run_efficiency_study(reps = 1000)
print(efficiency_summary)
plot_efficiency_results(efficiency_summary)
