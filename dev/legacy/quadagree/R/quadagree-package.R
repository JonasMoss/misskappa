#' @description
#' Provides functions for the estimating quadratically weighted agreement coefficients
#' Fleiss' kappa and Conger's kappa (multirater Cohen's kappa) with inference. It supports
#' both complete and incomplete data (missing completely at random), and offers
#' confidence intervals via the delta method or a studentized bootstrap.
#'
#' @importFrom Rcpp sourceCpp
#' @useDynLib quadagree, .registration = TRUE
#' @name quadagree-package
"_PACKAGE"

#' Confidence Intervals for Quadratic Kappa Coefficients
#'
#' Estimates the quadratic Fleiss' kappa, Conger's (Cohen's) kappa, and the
#' Brennan--Prediger coefficient for data in either long or aggregated form,
#' along with confidence intervals based on the delta method or
#' the studentized bootstrap.
#'
#' @section Methods:
#' - `kappa_raw()` applies to *long form* data where each row is a rater's scores. Supports data missing completely at random.
#' - `kappa_aggr()` applies to *aggregated form* data where each row contains counts per category. Supports data missing completely at random.
#'
#' Transformations are used to improve confidence interval coverage. Supported options include Fisher's *z*-transform, logit, and arcsin.
#'
#' The bootstrap procedure is a *studentized bootstrap* with variance estimated per replicate, ensuring second-order correctness (Efron, 1987).
#'
#' @param x A matrix or object coercible to a matrix. Each row is a subject, each column a rater (`kappa_raw`) or rating count (`kappa_aggr`).
#' @param values A numeric vector assigning a numeric value to each category (only used in `kappa_aggr`). Defaults to `1:ncol(x)`.
#' @param r Number of raters for `kappa_aggr`. Defaults to the maximum number of ratings per item.  Defaults to the maximal number of raters encountered in a row.
#' @param transform A variance-stabilizing transform to apply before CI construction. One of `"none"` (default), `"fisher"`, `"log"`, `"arcsin"`.
#' @param conf_level Confidence level for the CI. Default is `0.95`.
#' @param alternative Alternative hypothesis for CI. One of `"two.sided"` (default), `"greater"`, `"less"`.
#' @param bootstrap Logical. If `TRUE`, use the studentized bootstrap. Default is `FALSE`.
#' @param n_reps Number of bootstrap replicates (ignored if `bootstrap = FALSE`). Default is `1000`.
#' @param seed Optional random seed passed to bootstrap replicates.
#' @return An object of class `"quadagree_raw"` or `"quadagree_aggr"` with components for Fleiss' and Conger's kappa, including confidence intervals, standard errors, and relevant attributes.
#'
#' @references
#' Efron, B. (1987). Better Bootstrap Confidence Intervals. *JASA*, 82(397), 171–185.
#'
#' van Praag, B. M. S., Dijkstra, T. K., & van Velzen, J. (1985). Least-squares theory based on general distributional assumptions with an application to the incomplete observations problem. *Psychometrika*, 50(1), 25–36.
#'
#' Magnus, J. R., & Neudecker, H. (2019). *Matrix Differential Calculus with Applications in Statistics and Econometrics*. Wiley.
#'
#' Moss, J. & van Oest, R. (work in progress). Inference for quadratically weighted multi-rater kappas with missing raters.
#'
#' @examples
#' kappa_raw(dat.klein2018, transform = "fisher")
#' kappa_aggr(dat.fleiss1971)
#' @name kappa
#' @export
kappa_raw <- function(x,
                      transform = c("none", "fisher", "log", "arcsin"),
                      conf_level = 0.95,
                      alternative = c("two.sided", "greater", "less"),
                      bootstrap = FALSE,
                      n_reps = 1000,
                      seed = NULL) {
  transform <- match.arg(transform)
  alternative <- match.arg(alternative)
  data_name <- deparse(substitute(x))
  x <- as.matrix(x)

  rating_values <- sort(unique(stats::na.omit(c(x))))
  C <- length(rating_values)
  c1 <- (2 / C^2) * (C * sum(rating_values^2) - sum(rating_values)^2)

  res <- KappasRawCpp(
    x = x, c1 = c1, transform_str = transform, conf_level = conf_level,
    alternative_str = alternative, bootstrap = bootstrap, n_reps = n_reps, seed = seed
  )


  attr(res, "data.name") <- data_name
  attr(res, "n") <- nrow(x)
  attr(res, "method") <- if (bootstrap) "Studentized bootstrap CI" else "Delta method CI"
  attr(res, "transform") <- transform
  attr(res, "alternative") <- alternative
  attr(res, "conf_level") <- conf_level
  class(res) <- "quadagree_raw"

  res
}

#' @rdname kappa
#' @export
kappa_aggr <- function(x, values = NULL, r = NULL,
                       transform = c("none", "fisher", "log", "arcsin"),
                       conf_level = 0.95,
                       alternative = c("two.sided", "greater", "less"),
                       bootstrap = FALSE,
                       n_reps = 1000,
                       seed = NULL) {
  transform <- match.arg(transform)
  alternative <- match.arg(alternative)
  data_name <- deparse(substitute(x))
  x <- as.matrix(x)
  if (is.null(values)) values <- seq_len(ncol(x))
  if (is.null(r)) r <- max(rowSums(x, na.rm = TRUE))

  C <- length(values)
  c1 <- (2 / C^2) * (C * sum(values^2) - sum(values)^2)

  res <- KappasAggrCpp(
    x = x, values = values, R = r, c1 = c1, transform_str = transform,
    conf_level = conf_level, alternative_str = alternative,
    bootstrap = bootstrap, n_reps = n_reps, seed = seed
  )

  attr(res, "data.name") <- data_name
  attr(res, "n") <- nrow(x)
  attr(res, "r") <- r
  attr(res, "method") <- if (bootstrap) "Studentized bootstrap CI" else "Delta method CI"
  attr(res, "transform") <- transform
  attr(res, "alternative") <- alternative
  attr(res, "conf_level") <- conf_level
  class(res) <- "quadagree_aggr"

  res
}

#' @export
print.quadagree_raw <- function(x, ...) {
  data_name <- attr(x, "data.name")
  n_eff <- sprintf(" (n_eff = %d)", x[[1]]$n_eff)

  cat("Quadratic Agreement Coefficients (raw data)\n")
  cat("data: ", data_name, n_eff, "\n", sep = "")
  cat("method:", attr(x, "method"))
  if (attr(x, "transform") != "none") cat(sprintf(" (%s transform)", attr(x, "transform")))
  cat("\n\n")

  print_kappa <- function(obj, name) {
    if (is.null(obj)) {
      return()
    }
    ci <- sprintf("[%.3f, %.3f]", obj$interval[1], obj$interval[2])
    ci_label <- sprintf("CI (%.0f%%)", attr(x, "conf_level") * 100)
    cat(sprintf(
      "%-18s: estimate = % .3f, %s %s, std.err = %.4f\n",
      name, obj$estimate, ci_label, ci, obj$std_err
    ))
  }

  print_kappa(x$fleiss, "Fleiss' Kappa")
  print_kappa(x$conger, "Conger's Kappa")
  print_kappa(x$bp, "Brennan-Prediger")

  invisible(x)
}

#' @export
print.quadagree_aggr <- function(x, ...) {
  data_name <- attr(x, "data.name")
  n_eff <- sprintf(" (n_eff = %d)", x[[1]]$n_eff)

  cat("Quadratic Agreement Coefficients (aggregated data)\n")
  cat("data: ", data_name, n_eff, "\n", sep = "")
  cat("method:", attr(x, "method"))
  if (attr(x, "transform") != "none") cat(sprintf(" (%s transform)", attr(x, "transform")))
  cat("\n\n")

  print_kappa <- function(obj, name) {
    if (is.null(obj)) {
      return()
    }
    ci <- sprintf("[%.3f, %.3f]", obj$interval[1], obj$interval[2])
    ci_label <- sprintf("CI (%.0f%%)", attr(x, "conf_level") * 100)
    cat(sprintf(
      "%-18s: estimate = % .3f, %s %s, std.err = %.4f\n",
      name, obj$estimate, ci_label, ci, obj$std_err
    ))
  }

  print_kappa(x$fleiss, "Fleiss' Kappa")
  print_kappa(x$bp, "Brennan-Prediger")

  invisible(x)
}

#' Simulate ratings from a skill-difficulty agreement model
#'
#' @description
#' This function simulates a matrix of categorical ratings under a class of skill-based guessing models
#' described in Moss (2023). Each rater has a skill level between 0 and 1. For each item, the true class
#' is drawn from a true distribution, and raters either guess or respond correctly depending on their skill.
#'
#' Several models are supported:
#' - `"general"`: Fully user-specified `true_dist` and `guessing_dist`
#' - `"cohen-fleiss"`: True labels are marginally distributed like the weighted average of guessing distributions
#' - `"fleiss"`: Guessing distribution equals `true_dist` (classic Fleiss' kappa assumption)
#' - `"bp"`: Guessing is uniform, with user-specified `true_dist`
#' - `"tu"`: True distribution is uniform, with user-specified `guessing_dist`
#'
#' Agreement is governed by the skill vector `s`, and the true agreement is returned as an attribute
#' (`"kappa"`) using the latent agreement definition from Moss (2023).
#'
#' @param n Number of items to simulate.
#' @param s Numeric vector of rater skill levels between 0 and 1.
#' @param model One of `"general"`, `"cohen-fleiss"`, `"fleiss"`, `"bp"`, or `"tu"`.
#' @param true_dist Optional numeric vector of probabilities for the true class distribution.
#' @param guessing_dist Optional guessing distribution: either a vector (shared across raters) or a matrix (rater-specific).
#'
#' @return A matrix of simulated ratings (n rows, J columns) with attributes:
#' - `"n"`: number of items
#' - `"s"`: skill vector
#' - `"true_dist"`: normalized true class distribution
#' - `"guessing_dist"`: normalized guessing distribution
#' - `"kappa"`: true latent agreement value under the model
#'
#' @references
#' Moss, J. (2023). Measuring Agreement Using Guessing Models and Knowledge Coefficients.
#' \emph{Psychometrika}, \doi{10.1007/s11336-023-09887-2} https://arxiv.org/abs/2309.03613
#'
#' @keywords internal
simulate_jsm <- function(n,
                         s,
                         model = c("general", "cohen-fleiss", "fleiss", "bp", "tu"),
                         true_dist = NULL,
                         guessing_dist = NULL) {
  model <- match.arg(model)

  if (!is.numeric(n) || length(n) != 1 || n < 1 || n != as.integer(n)) {
    stop("n must be a single positive integer.")
  }

  if (!is.numeric(s) || any(s < 0 | s > 1)) {
    stop("s must be a numeric vector with values in [0, 1].")
  }

  J <- length(s) # number of raters

  # Normalize and check guessing_dist
  if (!is.null(guessing_dist)) {
    if (is.matrix(guessing_dist)) {
      if (nrow(guessing_dist) != J) {
        stop("If guessing_dist is a matrix, it must have one row per rater.")
      }
      if (any(guessing_dist < 0)) {
        stop("guessing_dist matrix must be non-negative.")
      }
      guessing_dist <- guessing_dist / rowSums(guessing_dist)
    } else {
      if (!is.numeric(guessing_dist)) {
        stop("guessing_dist must be a numeric vector or matrix.")
      }
      if (any(guessing_dist < 0)) {
        stop("guessing_dist must contain non-negative values.")
      }
      guessing_dist <- guessing_dist / sum(guessing_dist)
    }
  }

  # Normalize and check true_dist
  if (!is.null(true_dist)) {
    if (!is.numeric(true_dist) || any(true_dist < 0)) {
      stop("true_dist must be a non-negative numeric vector.")
    }
    true_dist <- true_dist / sum(true_dist)
  }

  # Determine number of categories
  q <- if (!is.null(true_dist)) {
    length(true_dist)
  } else if (!is.null(guessing_dist)) {
    if (is.matrix(guessing_dist)) ncol(guessing_dist) else length(guessing_dist)
  } else {
    stop("Either true_dist or guessing_dist must be supplied.")
  }

  # Model-specific logic
  if (model == "cohen-fleiss") {
    if (is.null(guessing_dist)) stop("guessing_dist must be provided for 'cohen-fleiss' model.")
    if (is.null(true_dist)) {
      if (is.matrix(guessing_dist)) {
        weights <- (1 - s) / sum(1 - s)
        true_dist <- colSums(guessing_dist * weights)
      } else {
        true_dist <- guessing_dist
      }
    }
  } else if (model == "tu") {
    true_dist <- rep(1 / q, q)
  } else if (model == "bp") {
    if (is.null(true_dist)) stop("true_dist must be provided for 'bp' model.")
    guessing_dist <- rep(1 / q, q)
  } else if (model == "fleiss") {
    if (is.null(true_dist)) stop("true_dist must be provided for 'fleiss' model.")
    guessing_dist <- true_dist
  }

  if (is.null(true_dist) || abs(sum(true_dist) - 1) > 1e-8) {
    stop("true_dist must sum to 1.")
  }
  if (is.null(guessing_dist)) {
    stop("guessing_dist could not be inferred.")
  }

  # Sample true labels
  x_star <- sample(q, n, replace = TRUE, prob = true_dist)

  # Simulate responses
  observations <- if (is.matrix(guessing_dist)) {
    sapply(seq_along(s), function(j) {
      z <- stats::rbinom(n, 1, s[j])
      z * x_star + (1 - z) * sample(q, n, replace = TRUE, prob = guessing_dist[j, ])
    })
  } else {
    sapply(seq_along(s), function(j) {
      z <- stats::rbinom(n, 1, s[j])
      z * x_star + (1 - z) * sample(q, n, replace = TRUE, prob = guessing_dist)
    })
  }

  # Latent agreement (true kappa)
  true_jsm <- function(s) {
    ss <- s %*% t(s)
    diag(ss) <- 0
    j <- length(s)
    sum(ss) / (j * (j - 1))
  }

  attr(observations, "n") <- n
  attr(observations, "s") <- s
  attr(observations, "true_dist") <- true_dist
  attr(observations, "guessing_dist") <- guessing_dist
  attr(observations, "kappa") <- true_jsm(s)
  observations
}


#' Population Values of Fleiss' and Conger's Kappa
#'
#' Computes the population values of Fleiss' and Conger's quadratic kappa
#' given the mean vector and covariance matrix of ordinal ratings.
#'
#' @param mu A numeric vector of length \eqn{r} giving the mean ratings per rater.
#' @param sigma A \eqn{r \times r} covariance matrix of the ratings.
#'
#' @return A named list with entries:
#'   \item{fleiss}{Population Fleiss' kappa value.}
#'   \item{conger}{Population Conger's (or Cohen's) kappa value.}
#'
#' @export
population_kappas <- function(mu, sigma) {
  if (!is.numeric(mu) || !is.vector(mu)) stop("`mu` must be a numeric vector.")
  if (!is.matrix(sigma) || !is.numeric(sigma)) stop("`sigma` must be a numeric matrix.")
  if (nrow(sigma) != length(mu) || ncol(sigma) != length(mu)) stop("Dimensions of `sigma` must match length of `mu`.")
  if (!isSymmetric(sigma)) stop("`sigma` must be symmetric.")

  r <- length(mu)
  trace <- sum(diag(sigma))

  # for conger: variance of means is scaled by r^2
  mean_diff_conger <- (mean(mu^2) - mean(mu)^2) * r^2
  top_conger <- sum(sigma) - trace
  bottom_conger <- (r - 1) * trace + mean_diff_conger
  conger <- top_conger / bottom_conger

  # for fleiss: variance of means is scaled by r
  mean_diff_fleiss <- (mean(mu^2) - mean(mu)^2) * r
  top_fleiss <- sum(sigma) - trace - mean_diff_fleiss
  bottom_fleiss <- (r - 1) * (trace + mean_diff_fleiss)
  fleiss <- top_fleiss / bottom_fleiss

  list(fleiss = fleiss, conger = conger)
}

#' Convert Raw Rating Matrix to Counts Format
#'
#' This function takes a matrix of integer ratings (subjects x raters)
#' and converts it into a counts matrix suitable for Fleiss' Kappa analysis.
#' Each row in the output represents a unique subject (by their pattern of ratings),
#' and columns represent the counts of ratings in each category.
#'
#' @param raw_data A matrix or data frame where rows are subjects and columns
#'   are raters. Categorical ratings should be integers (e.g., 1, 2, 3, ...).
#'   `NA` values are handled by being ignored in the counts for that subject.
#'
#' @return A matrix where each row is a subject and each column is a category.
#'   The values are the number of raters who assigned that subject to that category.
#'   The number of columns is determined by the maximum category label in the data.
#'
#' @export
#' @examples
#' raw <- matrix(c(1, 2, 1, NA, 1, 2, 3, 3, 3), nrow = 3, byrow = TRUE)
#' to_counts_matrix(raw)
#' #      [,1] [,2] [,3]
#' # [1,]    2    1    0  (2 ratings of '1', 1 rating of '2')
#' # [2,]    1    1    0  (1 rating of '1', 1 rating of '2', 1 NA)
#' # [3,]    0    0    3  (3 ratings of '3')
to_counts_matrix <- function(raw_data) {
  # Find the number of categories from the data
  # We assume categories are 1, 2, ..., C
  C <- max(raw_data, na.rm = TRUE)
  if (!is.finite(C)) {
    stop("Could not determine the number of categories. Is the data empty or all NA?")
  }

  # Use apply to process each row (subject)
  counts <- t(apply(raw_data, 1, function(row) {
    # 'tabulate' is very fast for this. It counts occurrences of integers.
    # We ignore NAs implicitly as tabulate doesn't see them.
    tabulate(row, nbins = C)
  }))

  # Set column names for clarity
  colnames(counts) <- paste0("cat_", 1:C)

  return(counts)
}

#' Introduce Missing Data into a Matrix
#'
#' @param data The complete data matrix.
#' @param prop The proportion of cells to set to NA.
#'
#' @return A matrix with missing values.
introduce_missingness <- function(data, prop = 0.1) {
  if (prop == 0) return(data)

  # Find the indices of non-NA cells to choose from
  valid_indices <- which(!is.na(data))

  # Determine how many NAs to introduce
  n_to_na <- floor(length(valid_indices) * prop)

  # Sample the indices to set to NA
  na_indices <- sample(valid_indices, size = n_to_na, replace = FALSE)

  # Introduce NAs
  data[na_indices] <- NA

  return(data)
}

introduce_missingness<- function(data, prop) {
  probs <- (1:ncol(data)) / sum((1:ncol(data)))
  # probs: vector of length ncol(data), giving per-variable missingness probs
  stopifnot(length(probs) == ncol(data))

  n <- nrow(data)
  p <- ncol(data)

  for (j in seq_len(p)) {
    missing_j <- rbinom(n, size = 1, prob = probs[j]) == 1
    data[missing_j, j] <- NA
  }

  return(data)
}


#' Run a Single Simulation Trial for Agreement Coefficients
#'
#' This function simulates data, optionally introduces missingness, and then
#' estimates agreement using conger_cpp, fleiss_cpp, and bp_cpp.
#'
#' @param n_subjects Number of subjects to simulate.
#' @param skill_vec A numeric vector of rater skills.
#' @param true_dist The true distribution of categories.
#' @param missing_prop Proportion of data to be made missing.
#' @param seed An optional random seed for reproducibility.
#'
#' @return A tidy data frame with one row containing the true kappa and
#'   the estimates from the three methods.
run_sim_trial <- function(n_subjects, skill_vec, true_dist, missing_prop = 0.1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # 1. Simulate complete data
  # We use the "bp" model (uniform guessing) as a simple, clear case.
  # The true kappa is only a function of skills, so the model choice doesn't
  # affect the target value, only the data generation process.
  sim_data <- simulate_jsm(
    n = n_subjects,
    s = skill_vec,
    model = "fleiss",
    true_dist = true_dist
  )

  true_kappa <- attr(sim_data, "kappa")
  C <- length(true_dist) # Number of categories

  # 2. Introduce missingness
  sim_data_missing <- introduce_missingness(sim_data, prop = missing_prop)

  # The C++ functions expect integer matrices and NA_integer_
  # The `introduce_missingness` function above uses R's default `NA`,
  # which is a logical NA. We must convert it.
  storage.mode(sim_data_missing) <- "integer"
  #sim_data_missing[is.na(sim_data_missing)] <- -2147483648L # Our special NA value

  # 3. Define standard options for C++ functions
  em_opts <- list(tol = 1e-8, max_iter = 1000, prune_tol = 1e-12, start_alpha = 0.5)

  # 4. Estimate with Conger's Kappa (raw data model)
  conger_est <- tryCatch({
    conger_cpp(
      x = sim_data_missing,
      c = C,
      weight_type = "identity",
      values = NULL,
      em_options = em_opts,
      analysis_options = NULL
    )$estimate
  }, error = function(e) NA)

  fleiss_raw_est <- tryCatch({
    fleiss_raw_cpp(
      x = sim_data_missing,
      c = C,
      weight_type = "identity",
      values = NULL,
      em_options = em_opts,
      analysis_options = NULL
    )$estimate
  }, error = function(e) NA)

  bp_raw_est <- tryCatch({
    bp_raw_cpp(
      x = sim_data_missing,
      c = C,
      weight_type = "identity",
      values = NULL,
      em_options = em_opts,
      analysis_options = NULL
    )$estimate
  }, error = function(e) NA)

  # 5. Convert to counts format for Fleiss/BP
  # Important: Convert back to R's NA first for `to_counts_matrix` to work
  sim_data_for_counts <- sim_data_missing
  #sim_data_for_counts[sim_data_for_counts == -2147483648L] <- NA
  counts_matrix <- as.matrix(to_counts_matrix(sim_data_for_counts))

  r_modal <- length(skill_vec)

  # 6. Estimate with Fleiss' Kappa (counts data model)
  fleiss_est <- tryCatch({
    fleiss_cpp(
      x = counts_matrix,
      r = r_modal,
      weight_type = "quadratic",
      values = NULL,
      em_options = em_opts,
      analysis_options = NULL
    )$estimate
  }, error = function(e) NA)

  # 7. Estimate with Brennan-Prediger (counts data model)
  bp_est <- tryCatch({
    bp_cpp(
      x = counts_matrix,
      r = r_modal,
      weight_type = "identity",
      values = NULL,
      em_options = em_opts,
      analysis_options = NULL
    )$estimate
  }, error = function(e) NA)

  # 8. Return tidy data frame
  data.frame(
    true_kappa = true_kappa,
    conger_est = conger_est,
    fleiss_est = fleiss_raw_est,
    bp_est = bp_raw_est,
    n_subjects = n_subjects,
    missing_prop = missing_prop
  )
}
