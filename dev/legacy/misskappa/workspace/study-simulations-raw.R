#' Run a Full Consistency and Variance Study for Raw and Counts Data Estimators
#'
#' Evaluates consistency, standard error accuracy, and CI coverage of agreement
#' estimators under various conditions.
#'
#' @param reps Number of simulation replications for each parameter setting.
#' @param n_raters Number of raters to simulate. Defaults to 3.
#' @param n_categories Number of rating categories. Defaults to 3.
#' @return A data frame summarizing performance metrics for each estimator.
#' @keywords internal
run_consistency_study <- function(reps = 100, n_raters = 3, n_categories = 3) {

  param_grid <- expand.grid(
    n_subjects = c(100, 200, 800),
    missing_prop = c(0.2),
    missing_func_name = c("mcar_uniform", "mcar_by_rater"),
    exchangeable_raters = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )

  all_results <- list()

  for (i in seq_len(nrow(param_grid))) {
    params <- param_grid[i, ]
    cat(sprintf("Running setting %d of %d: n=%d, prop=%.1f, method=%s, exch=%s\n",
                i, nrow(param_grid), params$n_subjects, params$missing_prop,
                params$missing_func_name, params$exchangeable_raters))

    skill_vec <- if (params$exchangeable_raters) rep(0.8, n_raters) else seq(0.6, 0.9, length.out = n_raters)
    missing_func <- get(paste0("introduce_missing_", params$missing_func_name))
    missing_params <- if (params$missing_func_name == "mcar_by_rater") {
      list(props_by_rater = seq(0.1, params$missing_prop * 1.5, length.out = n_raters))
    } else {
      list(prop = params$missing_prop)
    }

    setting_results <- lapply(1:reps, function(rep_seed) {
      run_consistency_trial(
        n_subjects = params$n_subjects, skill_vec = skill_vec, n_categories = n_categories,
        missing_func = missing_func, missing_params = missing_params, seed = rep_seed * i
      )
    })

    combined_df <- do.call(rbind, setting_results)
    combined_df$n_subjects <- params$n_subjects
    combined_df$missing_prop <- params$missing_prop
    combined_df$missing_method <- params$missing_func_name
    combined_df$exchangeable <- params$exchangeable_raters

    all_results[[i]] <- combined_df
  }

  final_results <- do.call(rbind, all_results)

  # Calculate comprehensive summary stats for each group
  calculate_summary <- function(est, se, cov, true_kappa) {
    error <- est - true_kappa
    c(
      bias = mean(error, na.rm = TRUE),
      var = stats::var(est, na.rm = TRUE),
      mse = mean(error^2, na.rm = TRUE),
      avg_se = mean(se, na.rm = TRUE),
      emp_sd = stats::sd(est, na.rm = TRUE),
      coverage = mean(cov, na.rm = TRUE)
    )
  }

  summary_list <- by(final_results, final_results[, c("n_subjects", "missing_method", "exchangeable", "missing_prop")], function(group) {
    true_k <- group$true_kappa[1]

    s_conger <- calculate_summary(group$conger_est, group$conger_se, group$conger_coverage, true_k)
    s_fleiss_raw <- calculate_summary(group$fleiss_raw_est, group$fleiss_raw_se, group$fleiss_raw_coverage, true_k)
    s_bp_raw <- calculate_summary(group$bp_raw_est, group$bp_raw_se, group$bp_raw_coverage, true_k)
    s_fleiss_counts <- calculate_summary(group$fleiss_counts_est, group$fleiss_counts_se, group$fleiss_counts_coverage, true_k)
    s_bp_counts <- calculate_summary(group$bp_counts_est, group$bp_counts_se, group$bp_counts_coverage, true_k)

    data.frame(
      n_subjects = group$n_subjects[1], missing_method = group$missing_method[1],
      exchangeable = group$exchangeable[1], missing_prop = group$missing_prop[1],
      true_kappa = true_k,
      estimator = c("Conger", "Fleiss (raw)", "BP (raw)", "Fleiss (counts)", "BP (counts)"),
      bias = c(s_conger["bias"], s_fleiss_raw["bias"], s_bp_raw["bias"], s_fleiss_counts["bias"], s_bp_counts["bias"]),
      var = c(s_conger["var"], s_fleiss_raw["var"], s_bp_raw["var"], s_fleiss_counts["var"], s_bp_counts["var"]),
      mse = c(s_conger["mse"], s_fleiss_raw["mse"], s_bp_raw["mse"], s_fleiss_counts["mse"], s_bp_counts["mse"]),
      avg_se = c(s_conger["avg_se"], s_fleiss_raw["avg_se"], s_bp_raw["avg_se"], s_fleiss_counts["avg_se"], s_bp_counts["avg_se"]),
      emp_sd = c(s_conger["emp_sd"], s_fleiss_raw["emp_sd"], s_bp_raw["emp_sd"], s_fleiss_counts["emp_sd"], s_bp_counts["emp_sd"]),
      coverage = c(s_conger["coverage"], s_fleiss_raw["coverage"], s_bp_raw["coverage"], s_fleiss_counts["coverage"], s_bp_counts["coverage"])
    )
  })

  summary_df <- do.call(rbind, summary_list)
  rownames(summary_df) <- NULL
  return(summary_df)
}

#' Plot Consistency Study Results
#'
#' Generates plots visualizing the bias from a consistency study.
#' @param summary_df A data frame from `run_consistency_study`.
#' @return A `ggplot` object.
#' @keywords internal
plot_consistency_results <- function(summary_df) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for this plotting function.")
  }
  summary_df$exchangeable_label <- factor(summary_df$exchangeable,
                                          levels = c(TRUE, FALSE),
                                          labels = c("Raters: Exchangeable", "Raters: Non-Exchangeable"))
  ggplot2::ggplot(summary_df, ggplot2::aes(x = n_subjects, y = bias, color = estimator, linetype = as.factor(missing_prop))) +
    ggplot2::geom_line(alpha = 0.8) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray20") +
    ggplot2::facet_grid(exchangeable_label ~ missing_method, labeller = ggplot2::label_value) +
    ggplot2::scale_x_continuous(breaks = unique(summary_df$n_subjects)) +
    ggplot2::labs(title = "Estimator Bias vs. Sample Size", x = "Number of Subjects (n)",
                  y = "Average Bias (Estimate - True Kappa)", color = "Estimator", linetype = "Missing Prop.") +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom", strip.text = ggplot2::element_text(face = "bold"))
}

#' Plot SE and Coverage Results
#'
#' Generates plots for standard error accuracy and CI coverage.
#' @param summary_df A data frame from `run_consistency_study`.
#' @param conf_level The nominal confidence level for the coverage plot.
#' @return A list of two `ggplot` objects.
#' @keywords internal
plot_se_coverage_results <- function(summary_df, conf_level = 0.95) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required for this plotting function.")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("tidyr is required for this plotting function.")

  summary_df$exchangeable_label <- factor(summary_df$exchangeable,
                                          levels = c(TRUE, FALSE),
                                          labels = c("Raters: Exchangeable", "Raters: Non-Exchangeable"))

  se_data <- tidyr::pivot_longer(summary_df,
                                 cols = c("avg_se", "emp_sd"),
                                 names_to = "se_type",
                                 values_to = "se_value",
                                 names_prefix = "se_")
  se_data$se_type <- factor(se_data$se_type, levels = c("avg_se", "emp_sd"), labels = c("Estimated SE", "Empirical SD"))

  p1 <- ggplot2::ggplot(se_data, ggplot2::aes(x = n_subjects, y = se_value, color = estimator, linetype = se_type)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::facet_grid(exchangeable_label ~ missing_method, labeller = ggplot2::label_value) +
    ggplot2::scale_x_continuous(breaks = unique(summary_df$n_subjects)) +
    ggplot2::labs(title = "Standard Error Performance", x = "Number of Subjects (n)", y = "Standard Error",
                  color = "Estimator", linetype = "SE Type") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")

  p2 <- ggplot2::ggplot(summary_df, ggplot2::aes(x = n_subjects, y = coverage, color = estimator, shape = as.factor(missing_prop))) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_hline(yintercept = conf_level, linetype = "dashed", color = "black") +
    ggplot2::facet_grid(exchangeable_label ~ missing_method, labeller = ggplot2::label_value) +
    ggplot2::scale_x_continuous(breaks = unique(summary_df$n_subjects)) +
    ggplot2::scale_y_continuous(limits = c(min(0.85, min(summary_df$coverage, na.rm=T)), 1.0)) +
    ggplot2::labs(title = "Confidence Interval Coverage", x = "Number of Subjects (n)", y = "Coverage Probability",
                  color = "Estimator", shape = "Missing Prop.") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")

  return(list(se_plot = p1, coverage_plot = p2))
}

result <- run_consistency_study(2000)

plot_se_coverage_results(result)
