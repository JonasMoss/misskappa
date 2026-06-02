# This script conducts a definitive parameter recovery study using simulated data.
# It is designed to create a clear, verifiable example where misskappa's MLE
# succeeds and the naive available-case (AC) method fails demonstrably.
#
# The simulation establishes a "ground truth" by:
# 1. Defining a population model with high rater heterogeneity
#    (i.e., "expert" and "novice" raters). The true kappa for this population
#    is calculated theoretically.
# 2. Repeatedly sampling from this population and introducing a non-exchangeable
#    MAR pattern where the less-skilled "novice" raters are also more likely to
#    have missing data.
# 3. Comparing the sampling distribution of the MLE and AC estimators against
#    the single, known population kappa.

library(misskappa)
library(ggplot2)
library(tidyr)

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

#' Run the Full JSM-Based Parameter Recovery Study
#' @param reps Number of Monte Carlo repetitions.
#' @return A list containing the results data frame and the single ground truth kappa.
#' @keywords internal
run_jsm_recovery_study <- function(reps = 1000) {
  # --- FIX: Define the population parameters and ground truth ONCE ---
  n_subjects <- 500
  # Heterogeneous raters: two experts, two novices
  skill_vec <- c(1, 1, 0.5, 0.4)
  # Non-exchangeable missingness: novices are "lazy"
  missing_props <- c(0.1, 0.1, 0.4, 0.4)

  # Calculate the single, theoretical ground truth kappa from the model
  kappa_population <- attr(simulate_jsm(n=1, s=skill_vec, model = "fleiss", true_dist = c(0.5, 0.3, 0.2)), "kappa")
  # --- END FIX ---

  cat("Running JSM-based parameter recovery study...\n")

  results_list <- lapply(1:reps, function(i) {
    set.seed(i)
    if (i %% 100 == 0) cat(sprintf("  - Repetition %d of %d\n", i, reps))

    # 1. Simulate a new sample from the defined population
    data_complete <- simulate_jsm(
      n = n_subjects, s = skill_vec, model = "fleiss", true_dist = c(0.5, 0.3, 0.2)
    )

    # 2. Introduce the MAR pattern
    data_missing <- introduce_mar_by_rater(data_complete, missing_props)

    # 3. Estimate kappa with both methods
    est_mle <- tryCatch(conger_kappa(data_missing)$estimate, error = function(e) NA)
    est_ac <- tryCatch(conger_kappa_naive(data_missing), error = function(e) NA)

    data.frame(
      estimator = c("MLE", "Available-Case"),
      estimate = c(est_mle, est_ac)
    )
  })

  results_df <- do.call(rbind, results_list)

  list(results = results_df,
       ground_truth = kappa_population)
}

#' Plot the Results of the JSM Recovery Study
#' @param study_output The list returned by `run_jsm_recovery_study`.
#' @return A ggplot object.
#' @keywords internal
plot_jsm_recovery_results <- function(study_output) {

  results <- study_output$results
  kappa_true <- study_output$ground_truth

  results$estimator <- factor(results$estimator, levels = c("MLE", "Available-Case"))

  summary_df <- aggregate(estimate ~ estimator, data = results, FUN = mean, na.rm = TRUE)
  names(summary_df)[names(summary_df) == "estimate"] <- "mean_est"

  print_summary <- summary_df
  print_summary$bias <- print_summary$mean_est - kappa_true
  rownames(print_summary) <- print_summary$estimator
  print_summary <- print_summary[c("mean_est", "bias")]

  cat("\n--- JSM-Based Parameter Recovery Study Summary ---\n")
  cat(sprintf("Population Ground Truth Kappa: %.4f\n", kappa_true))
  print(print_summary)
  cat("--------------------------------------------------\n\n")

  ggplot(results, aes(x = estimate)) +
    geom_density(aes(fill = estimator, linetype = estimator), alpha = 0.5, bw = 0.01) +

    geom_vline(xintercept = kappa_true, linetype = "solid", color = "black", linewidth = 1) +
    geom_vline(data = summary_df, aes(xintercept = mean_est, linetype = estimator),
               color = "gray30", linewidth = 0.8) +

    scale_fill_manual(values = c("MLE" = "white", "Available-Case" = "white")) +
    scale_linetype_manual(values = c("MLE" = "solid", "Available-Case" = "dashed")) +

    annotate("text", x = kappa_true, y = Inf, label = "Ground Truth  ", vjust = 1.5, hjust = 1,
             fontface = "plain", size = 4) +

    labs(
      x = "Conger's Kappa Estimate",
      y = "Density",
      fill = "Estimator",
      linetype = "Estimator"
    ) +
    coord_cartesian(expand = FALSE, clip = "off") +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = c(0.18, 0.85),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.title = element_blank(),
      axis.line = element_line(colour = "black"),
      plot.margin = margin(10, 10, 10, 10)
    )
}

# --- Run the study and plot the results ---
recovery_results_jsm <- run_jsm_recovery_study(reps = 1000)
#plot_jsm_recovery_results(recovery_results_jsm)

mle <- subset(recovery_results_jsm$results, estimator == "MLE")$estimate
ac  <- subset(recovery_results_jsm$results, estimator == "Available-Case")$estimate

theta <- recovery_results_jsm$ground_truth
plot(density(ac, adjust = 2), lwd = 2, lty = 2, main = "", bty = "l", xlab = "Conger's kappa estimate",
     xlim = c(0.4, 0.65))
lines(density(mle, adjust = 2), lwd = 2)
#axis(side = 1, at = theta, labels = format(theta, digits = 3, nsmall = 3), lwd = 1, lty = 3)
abline(v = theta, lty = 3)
legend("topleft",
       legend = c("MLE", "Available case", "True kappa"),
       lty = c(1, 2, 3),
       lwd = c(2, 2, 1),
       bty = "n")

#abline(v = recovery_results_jsm$ground_truth)
#abline(v = mean(ac), lty = 2)




