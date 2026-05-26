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
#' @keywords internal
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
