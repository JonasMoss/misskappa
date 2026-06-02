#' Introduce Missing Data (MCAR Uniform)
#'
#' Sets a proportion of cells in a data matrix to NA, with each cell having
#' an equal and independent probability of being selected.
#'
#' @param data The complete data matrix.
#' @param prop The proportion of cells to set to NA.
#' @return A matrix with missing values.
#' @keywords internal
introduce_missing_mcar_uniform <- function(data, prop) {
  if (prop == 0) return(data)
  n_elements <- length(data)
  n_to_na <- floor(n_elements * prop)
  na_indices <- sample(n_elements, size = n_to_na, replace = FALSE)
  data[na_indices] <- NA
  return(data)
}

#' Introduce Missing Data (MCAR by Rater)
#'
#' Sets cells in a data matrix to NA based on rater-specific probabilities.
#' This models a scenario where some raters are more likely to have missing
#' entries than others, but the missingness is independent of the ratings.
#'
#' @param data The complete data matrix.
#' @param props_by_rater A numeric vector with a missingness probability for
#'   each rater (column).
#' @return A matrix with missing values.
#' @keywords internal
introduce_missing_mcar_by_rater <- function(data, props_by_rater) {
  if (ncol(data) != length(props_by_rater)) {
    stop("Length of 'props_by_rater' must match the number of columns in 'data'.")
  }
  for (j in seq_len(ncol(data))) {
    n_to_na <- floor(nrow(data) * props_by_rater[j])
    if (n_to_na > 0) {
      na_indices <- sample(nrow(data), size = n_to_na, replace = FALSE)
      data[na_indices, j] <- NA
    }
  }
  return(data)
}

#' Introduce Missing Data (MAR by Row Mean)
#'
#' Induces a missing-at-random (MAR) pattern where the probability of a cell
#' being missing depends on the mean of the other observed values for that subject.
#'
#' @param data The complete data matrix.
#' @param prop The target overall proportion of missing data.
#' @param strength A numeric value controlling how strongly the row mean
#'   influences the missingness probability.
#' @return A matrix with missing values.
#' @keywords internal
introduce_missing_mar_by_row_mean <- function(data, prop, strength = 2.0) {
  if (prop == 0) return(data)
  row_means <- rowMeans(data, na.rm = TRUE)

  scaled_means <- scale(row_means)

  objective <- function(b0) {
    probs <- 1 / (1 + exp(-(b0 + strength * scaled_means)))
    mean(probs) - prop
  }

  b0_solution <- tryCatch(
    stats::uniroot(objective, interval = c(-10, 10))$root,
    error = function(e) {
      warning("Could not solve for intercept in MAR generation; using approximation.")
      log(prop / (1 - prop))
    }
  )

  missing_probs <- 1 / (1 + exp(-(b0_solution + strength * scaled_means)))

  missing_mask <- matrix(
    stats::rbinom(length(data), 1, prob = rep(missing_probs, ncol(data))),
    nrow = nrow(data),
    byrow = FALSE
  ) == 1

  data[missing_mask] <- NA
  return(data)
}

#' Run a Single Raw and Counts Data Trial
#'
#' Simulates one dataset, introduces missingness, and computes agreement estimates
#' for both raw and counts-based functions.
#'
#' @inheritParams run_consistency_trial
#' @return A single-row data frame with true kappa, estimates, standard errors,
#'   and CI coverage indicators for all five relevant estimators.
#' @keywords internal
run_consistency_trial <- function(n_subjects,
                                  skill_vec,
                                  n_categories,
                                  missing_func,
                                  missing_params,
                                  seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # 1. Simulate complete data
  true_dist <- rep(1 / n_categories, n_categories)
  sim_data_complete <- simulate_jsm(
    n = n_subjects,
    s = skill_vec,
    model = "fleiss",
    true_dist = true_dist
  )
  true_kappa <- attr(sim_data_complete, "kappa")

  # 2. Introduce missingness
  missing_args <- c(list(data = sim_data_complete), missing_params)
  sim_data_missing <- do.call(missing_func, missing_args)

  # 3. Helper to run an estimator and extract results
  run_and_extract <- function(estimator_fun, ...) {
    tryCatch({
      res <- estimator_fun(...)
      is_covered <- res$ci[1] <= true_kappa && true_kappa <= res$ci[2]
      analytical_se <- attr(res, "full_results")$se
      list(est = res$estimate, se = analytical_se, coverage = is_covered)
    }, error = function(e) list(est = NA, se = NA, coverage = NA))
  }

  # 4. Estimate agreement coefficients
  analysis_opts <- analysis_control(bootstrap_method = "none")

  # Raw data estimators
  conger_res <- run_and_extract(conger_kappa, x = sim_data_missing, analysis_control = analysis_opts)
  fleiss_raw_res <- run_and_extract(fleiss_kappa_raw, x = sim_data_missing, analysis_control = analysis_opts)
  bp_raw_res <- run_and_extract(brennan_prediger_raw, x = sim_data_missing, analysis_control = analysis_opts)

  # Counts data estimators
  counts_data <- to_counts_matrix(sim_data_missing)
  r_total <- ncol(sim_data_complete)

  fleiss_counts_res <- run_and_extract(fleiss_kappa_counts, x = counts_data, r = r_total, analysis_control = analysis_opts)
  bp_counts_res <- run_and_extract(brennan_prediger_counts, x = counts_data, r = r_total, analysis_control = analysis_opts)

  # 5. Return results
  data.frame(
    true_kappa = true_kappa,
    conger_est = conger_res$est, conger_se = conger_res$se, conger_coverage = conger_res$coverage,
    fleiss_raw_est = fleiss_raw_res$est, fleiss_raw_se = fleiss_raw_res$se, fleiss_raw_coverage = fleiss_raw_res$coverage,
    bp_raw_est = bp_raw_res$est, bp_raw_se = bp_raw_res$se, bp_raw_coverage = bp_raw_res$coverage,
    fleiss_counts_est = fleiss_counts_res$est, fleiss_counts_se = fleiss_counts_res$se, fleiss_counts_coverage = fleiss_counts_res$coverage,
    bp_counts_est = bp_counts_res$est, bp_counts_se = bp_counts_res$se, bp_counts_coverage = bp_counts_res$coverage
  )
}
