# This script runs a master simulation study to generate a comprehensive
# table comparing the confidence interval coverage of the misskappa MLE and a
# naive Available-Case (AC) estimator using its standard analytical variance.
#
# The goal is to provide a definitive, publication-ready table that demonstrates
# the robustness of the MLE and the failures of the AC method under realistic
# data conditions.

library(misskappa)
library(tidyr)
library(knitr)
library(irrCAC) # Required for the analytical AC variance function

#' Naive Conger's Kappa with Analytical SE (from irrCAC)
#'
#' A modified version of the conger() function to calculate the naive
#' available-case kappa and its analytical standard error.
#' @return A named vector c(kappa=, se=).
#' @keywords internal
conger_naive_analytical <- function(ratings.mat, weights = "unweighted", conflev = 0.95) {
  n <- nrow(ratings.mat); r <- ncol(ratings.mat); f <- 0
  categ <- sort(unique(na.omit(as.vector(ratings.mat)))); q <- length(categ)

  if (q <= 1) return(c(kappa = NA, se = NA)) # Handle single category case

  if (is.character(weights)) {
    weights.mat <- switch(weights,
                          "quadratic" = irrCAC::quadratic.weights(categ),
                          "ordinal" = irrCAC::ordinal.weights(categ),
                          "linear" = irrCAC::linear.weights(categ),
                          "radical" = irrCAC::radical.weights(categ),
                          "ratio" = irrCAC::ratio.weights(categ),
                          "circular" = irrCAC::circular.weights(categ),
                          "bipolar" = irrCAC::bipolar.weights(categ),
                          irrCAC::identity.weights(categ) # Default
    )
  } else { weights.mat <- weights }

  agree.mat <- matrix(0, nrow = n, ncol = q)
  for (k in 1:q) {
    categ.is.k <- (ratings.mat == categ[k])
    agree.mat[, k] <- (replace(categ.is.k, is.na(categ.is.k), FALSE)) %*% rep(1, r)
  }
  agree.mat.w <- t(weights.mat %*% t(agree.mat))
  classif.mat <- matrix(0, nrow = r, ncol = q)
  for (k in 1:q) {
    with.mis <- (t(ratings.mat) == categ[k])
    without.mis <- replace(with.mis, is.na(with.mis), FALSE)
    classif.mat[, k] <- without.mis %*% rep(1, n)
  }
  ri.vec <- agree.mat %*% rep(1, q)
  sum.q <- (agree.mat * (agree.mat.w - 1)) %*% rep(1, q)
  n2more <- sum(ri.vec >= 2)
  if (n2more == 0) return(c(kappa = NA, se = NA))

  pa <- sum(sum.q[ri.vec >= 2] / ((ri.vec * (ri.vec - 1))[ri.vec >= 2])) / n2more
  ng.vec <- classif.mat %*% rep(1, q)
  pgk.mat <- classif.mat / (ng.vec %*% rep(1, q))
  pgk.mat[is.nan(pgk.mat)] <- 0

  p.mean.k <- (t(pgk.mat) %*% rep(1, r)) / r
  s2kl.mat <- (t(pgk.mat) %*% pgk.mat - r * p.mean.k %*% t(p.mean.k)) / (r-1)
  pe <- sum(weights.mat * (p.mean.k %*% t(p.mean.k) - s2kl.mat / r))

  if (abs(1 - pe) < 1e-9) return(c(kappa = 1, se = 0))
  conger.kappa <- (pa - pe) / (1 - pe)

  # Variance Calculation
  bkl.mat <- (weights.mat + t(weights.mat)) / 2
  lamda.ig.mat <- matrix(0, n, r)
  epsi.ig.mat <- 1 * !is.na(ratings.mat)

  for (k in 1:q) {
    lamda.ig.kmat <- matrix(0, n, r)
    for (l in 1:q) {
      delta.ig.mat <- 1 * (ratings.mat == categ[l]); delta.ig.mat[is.na(delta.ig.mat)] <- 0
      term_in_paren <- (delta.ig.mat - (epsi.ig.mat - rep(1, n) %*% t(ng.vec/n)) * (rep(1, n) %*% t(pgk.mat[,l])))
      lamda.ig.kmat <- lamda.ig.kmat + weights.mat[k, l] * term_in_paren
    }
    lamda.ig.kmat <- lamda.ig.kmat * (rep(1, n) %*% t(n/ng.vec))
    lamda.ig.kmat[is.nan(lamda.ig.kmat) | is.infinite(lamda.ig.kmat)] <- 0
    lamda.ig.mat <- lamda.ig.mat + lamda.ig.kmat * (r * mean(pgk.mat[,k]) - rep(1, n) %*% t(pgk.mat[, k]))
  }
  pe.ivec <- (lamda.ig.mat %*% rep(1, r)) / (r * (r - 1))
  den.ivec <- ri.vec * (ri.vec - 1); den.ivec <- den.ivec - (den.ivec == 0)
  pa.ivec <- sum.q / den.ivec
  pe.r2 <- pe * (ri.vec >= 2)
  conger.ivec <- (n/n2more) * (pa.ivec - pe.r2)/(1-pe)
  conger.ivec.x <- conger.ivec - 2*(1-conger.kappa)*(pe.ivec-pe)/(1-pe)
  var.conger <- ((1-f)/(n*n)) * sum((conger.ivec.x - conger.kappa)^2)

  return(c(kappa = conger.kappa, se = sqrt(var.conger)))
}

#' Helper to construct a Fisher-transformed CI
#' @keywords internal
construct_fisher_ci <- function(kappa, se, conf.level = 0.95) {
  if (is.na(kappa) || is.na(se) || se <= 1e-9 || abs(kappa) >= 1) return(c(NA, NA))
  z_crit <- qnorm(1 - (1 - conf.level) / 2)
  z_val <- 0.5 * log((1 + kappa) / (1 - kappa))
  se_z <- se / (1 - kappa^2)
  ci_z <- z_val + c(-1, 1) * z_crit * se_z
  ci_kappa <- (exp(2 * ci_z) - 1) / (exp(2 * ci_z) + 1)
  return(ci_kappa)
}
#' Introduce Missing Data based on various mechanisms
#'
#' This function introduces missing data into a dataset. It can handle
#' simple independent cell-level MCAR or complex joint MCAR mechanisms.
#'
#' @param data The complete data frame.
#' @param mechanism_def A definition for the missingness mechanism.
#'   - If a numeric vector, it's treated as probabilities for independent
#'     cell-level MCAR for each rater.
#'   - If a multi-dimensional array, it's treated as the joint probability
#'     tensor P(M_1=m1, M_2=m2, ...), where m_i is 0 (missing) or 1 (present).
#'     The dimensions must be 2x2x...x2.
#' @return A data frame with NAs introduced.
#' @keywords internal
introduce_missingness <- function(data, mechanism_def) {
  n_subjects <- nrow(data)
  n_raters <- ncol(data)
  data_missing <- data

  if (is.array(mechanism_def) && length(dim(mechanism_def)) == n_raters) {
    # --- New Mode: Joint Missingness Distribution via Probability Tensor ---
    # Validate the tensor
    if (any(dim(mechanism_def) != 2)) {
      stop("All dimensions of the probability tensor must be 2.")
    }
    if (abs(sum(mechanism_def) - 1) > 1e-9) {
      stop("The sum of probabilities in the tensor must be 1.")
    }

    # Create a mapping from a single index (1 to 2^n_raters) to the patterns
    # This is more efficient than building the full patterns matrix
    all_patterns <- expand.grid(replicate(n_raters, c(1, 0), simplify = FALSE))
    # Note: R fills arrays column-major, so expand.grid's order matches the
    #       flattened array's order. We use c(1,0) to match indices:
    #       index 1 -> m_i=1 (present), index 2 -> m_i=0 (missing)

    # Sample a pattern index for each subject
    # The probability vector is the flattened tensor
    sampled_indices <- sample.int(2^n_raters, size = n_subjects, replace = TRUE, prob = c(mechanism_def))
    chosen_patterns <- all_patterns[sampled_indices, , drop = FALSE]

    # Create a logical mask of where to place NAs (where the pattern is 0)
    na_mask <- (chosen_patterns == 0)
    data_missing[na_mask] <- NA

  } else if (is.vector(mechanism_def) && !is.list(mechanism_def)) {
    # --- Legacy Mode: Independent Cell-Level MCAR ---
    for (j in 1:n_raters) {
      # Note: This is P(missing), while the tensor uses P(present/missing)
      # To be consistent, let's use the same logic
      p_missing <- mechanism_def[j]
      na_mask <- sample(c(TRUE, FALSE), size = n_subjects, replace = TRUE, prob = c(p_missing, 1 - p_missing))
      data_missing[na_mask, j] <- NA
    }
  } else {
    stop("Invalid 'mechanism_def' format. Must be a numeric vector or a multi-dimensional array.")
  }

  return(data_missing)
}


#' Run a Single Simulation Trial for CI Coverage
#' @keywords internal
run_ci_trial <- function(n_subjects, skill_vec, missing_props, kappa_population, seed) {
  set.seed(seed)
  data_complete <- simulate_jsm(n = n_subjects, s = skill_vec, model = "fleiss", true_dist = c(0.5, 0.3, 0.2))
  data_missing <- introduce_missingness(data_complete, missing_props)

  # MLE Confidence Interval
  fit_mle <- tryCatch(conger_kappa(data_missing, weight = "unweighted"), error = function(e) NULL)
  coverage_mle <- NA
  if (!is.null(fit_mle)) {
    ci_mle <- construct_fisher_ci(fit_mle$estimate, fit_mle$se)
    if (all(!is.na(ci_mle))) coverage_mle <- kappa_population >= ci_mle[1] && kappa_population <= ci_mle[2]
  }

  # Available-Case Analytical CI
  fit_ac <- tryCatch(conger_naive_analytical(data_missing), error = function(e) c(kappa=NA, se=NA))
  coverage_ac <- NA
  if (!any(is.na(fit_ac))) {
    ci_ac <- construct_fisher_ci(fit_ac["kappa"], fit_ac["se"])
    if (all(!is.na(ci_ac))) coverage_ac <- kappa_population >= ci_ac[1] && kappa_population <= ci_ac[2]
  }
  data.frame(estimator = c("MLE", "AC"), coverage = c(coverage_mle, coverage_ac))
}

#' Run the Master CI Coverage Study
#' @keywords internal
run_master_ci_study <- function(reps = 1000) {


  p_tensor <- array(0, dim = rep(2, 4))

  p_tensor[1, 1, 1, 1] <- 0.65
  p_tensor[1, 2, 2, 1] <- 0.15
  p_tensor[2, 2, 1, 1] <- 0.10
  p_tensor[1, 1, 1, 2] <- 0.10

  c1 <- 0.7
  c12 <- 0.5
  c123 <- 0.4
  p_unif_pmcar_tensor <- array(0, dim = rep(2, n_raters))

  p_111 <- c123
  p_110 <- c12 - c123
  p_101 <- c12 - c123 # using c12=c13
  p_011 <- c12 - c123 # using c12=c23
  p_100 <- c1 - (p_110 + p_101 + p_111)
  p_010 <- c1 - (p_110 + p_011 + p_111)
  p_001 <- c1 - (p_101 + p_011 + p_111)
  p_000 <- 1 - (p_111+p_110+p_101+p_011+p_100+p_010+p_001)

  # Assign to the tensor (remembering index 1 is present, 2 is missing)
  p_unif_pmcar_tensor[1,1,1] <- p_111
  p_unif_pmcar_tensor[1,1,2] <- p_110
  p_unif_pmcar_tensor[1,2,1] <- p_101
  p_unif_pmcar_tensor[2,1,1] <- p_011
  p_unif_pmcar_tensor[1,2,2] <- p_100
  p_unif_pmcar_tensor[2,1,2] <- p_010
  p_unif_pmcar_tensor[2,2,1] <- p_001
  p_unif_pmcar_tensor[2,2,2] <- p_000



  case_definitions <- list(
    "Case 1 (Exchangeable)" = list(skills = c(0.3, 0.5, 0.7), missing = p_unif_pmcar_tensor)
    #"Case 2 (Moderate Het.)" = list(skills = seq(0.6, 0.9, length.out = 4), missing = seq(0.1, 0.4, length.out = 4)),
    #"Case 3 (High Het. MAR)" = list(skills = c(0.9, 0.9, 0.5, 0.5), missing = c(0.1, 0.1, 0.4, 0.4)),
    #"Case 4 (High Het. MCAR)" = list(skills = c(0.9, 0.9, 0.5, 0.5), missing = rep(0.25, 4))
  )
  param_grid <- tidyr::crossing(case_name = names(case_definitions), n_subjects = c(20, 50, 100, 500, 1000))
  cat("Running Master CI Coverage Study...\n")
  full_results <- list()
  for (i in 1:nrow(param_grid)) {
    params <- param_grid[i, ]; case_params <- case_definitions[[params$case_name]]
    kappa_pop <- attr(simulate_jsm(n = 1, s = case_params$skills, model = "fleiss", true_dist = c(0.5, 0.3, 0.2)), "kappa")
    cat(sprintf("  - Running: %s, n = %d\n", params$case_name, params$n_subjects))
    reps_list <- lapply(1:reps, function(rep_seed) {
      run_ci_trial(n_subjects = params$n_subjects, skill_vec = case_params$skills,
                   missing_props = case_params$missing, kappa_population = kappa_pop, seed = rep_seed * i)
    })
    results_df <- do.call(rbind, reps_list)
    results_df$case_name <- params$case_name; results_df$n_subjects <- params$n_subjects
    full_results[[i]] <- results_df
  }
  final_df <- do.call(rbind, full_results)
  aggregate(coverage ~ case_name + estimator + n_subjects, data = final_df, FUN = mean, na.rm = TRUE)
}

#' Format and Print Results Table
#' @return A formatted kable object.
#' @keywords internal
format_results_to_table <- function(summary_df) {
  summary_df$row_label <- paste(summary_df$case_name, summary_df$estimator)
  summary_df$coverage_char <- sprintf("%.2f", summary_df$coverage)
  table_wide <- tidyr::pivot_wider(
    summary_df[, c("row_label", "n_subjects", "coverage_char")],
    names_from = n_subjects, values_from = coverage_char)
  table_wide <- table_wide[order(table_wide$row_label), ]
  table_wide$Case <- gsub(" .*", "", table_wide$row_label)
  table_wide$Method <- gsub(".* (MLE|AC)", "\\1", table_wide$row_label)
  final_table <- table_wide[, c("Case", "Method", "20", "50", "100", "500", "1000")]
  colnames(final_table) <- c("Case", "Method", "n = 20", "n = 50", "n = 100", "n = 500", "n = 1000")
  kable(final_table, format = "pipe", align = "lccccc",
        caption = "Empirical 95% CI Coverage Probability for MLE and Available-Case (AC) Estimators")
}

# --- Run the study and print the results table ---
#coverage_results <- run_master_ci_study(reps = 10000)
#print(coverage_results)
#format_results_to_table(coverage_results)


coverage_results2 <- run_master_ci_study(reps = 10)
print(coverage_results2)
format_results_to_table(coverage_results2)
