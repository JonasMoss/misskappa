# --- Judge-Skill Model (JSM) ---

#' Calculate nu for the Judge-Skill Model (Internal)
#' @param s A numeric vector of judge skills in [0, 1].
#' @return The scalar value of nu.
#' @keywords internal
nu_jsm <- function(s) {
  j <- length(s)
  if (j < 2) {
    return(NaN)
  }

  # Efficiently calculates sum of all pairwise products of unique pairs
  # (sum(s)^2 - sum(s^2)) is twice the sum of s_i*s_j for i < j
  n_pairs <- j * (j - 1)

  sum_of_prods <- (sum(s)^2 - sum(s^2))

  # The mean of the products
  sum_of_prods / n_pairs
}

#' Simulate ratings from the Judge-Skill Model (Internal)
#' @param n Number of items.
#' @param s A numeric vector of judge skills.
#' @param true_dist Vector of true category probabilities.
#' @param guess_dist Vector or Matrix of guessing probabilities.
#' @return A matrix of simulated ratings.
#' @keywords internal
sim_jsm <- function(n, s, true_dist, guess_dist) {
  c <- length(true_dist)
  j <- length(s)

  x_star <- sample(seq_len(c), n, replace = TRUE, prob = true_dist)
  ratings <- matrix(NA_integer_, nrow = n, ncol = j)

  is_guess_matrix <- is.matrix(guess_dist)

  for (k in seq_len(j)) {
    is_known <- stats::rbinom(n, 1, prob = s[k]) == 1
    n_guess <- sum(!is_known)

    ratings[is_known, k] <- x_star[is_known]

    if (n_guess > 0) {
      current_guess_dist <- if (is_guess_matrix) guess_dist[k, ] else guess_dist
      ratings[!is_known, k] <- sample(
        seq_len(c), n_guess, replace = TRUE, prob = current_guess_dist
      )
    }
  }
  ratings
}


# --- Item-Difficulty Model (IDM) ---

#' Calculate nu for the Item-Difficulty Model (Internal)
#' @param beta_params A list with shape parameters `a` and `b` for a Beta distribution.
#' @return The scalar value of nu.
#' @keywords internal
nu_idm <- function(beta_params) {
  a <- beta_params$a
  b <- beta_params$b

  # E[delta^2] for delta ~ Beta(a, b)
  (a * (a + 1)) / ((a + b) * (a + b + 1))
}

#' Generate item easiness parameters for IDM (Internal)
#' @param n Number of items.
#' @param beta_params List with shape parameters `a` and `b`.
#' @return A numeric vector of item easiness values (`delta_i`).
#' @keywords internal
gen_idm <- function(n, beta_params) {
  stats::rbeta(n, shape1 = beta_params$a, shape2 = beta_params$b)
}

#' Simulate ratings from the Item-Difficulty Model (Internal)
#' @param n Number of items.
#' @param j Number of judges.
#' @param delta_i A numeric vector of item easiness values.
#' @param true_dist Vector of true category probabilities.
#' @param guess_dist Vector or Matrix of guessing probabilities.
#' @return A matrix of simulated ratings.
#' @keywords internal
sim_idm <- function(n, j, delta_i, true_dist, guess_dist) {
  c <- length(true_dist)
  x_star <- sample(seq_len(c), n, replace = TRUE, prob = true_dist)
  ratings <- matrix(NA_integer_, nrow = n, ncol = j)

  is_guess_matrix <- is.matrix(guess_dist)

  # For this model, all judges either know or guess for a given item
  is_known <- stats::rbinom(n, 1, prob = delta_i) == 1

  # Assign known ratings
  # sweep() is an efficient way to fill rows of a matrix with a vector
  if(any(is_known)) {
    ratings[is_known, ] <- x_star[is_known]
  }

  # Assign guessed ratings
  items_to_guess <- which(!is_known)
  if (length(items_to_guess) > 0) {
    if (is_guess_matrix) {
      # Slower, but necessary for judge-specific guessing
      for (k in seq_len(j)) {
        ratings[items_to_guess, k] <- sample(
          seq_len(c), length(items_to_guess), replace = TRUE, prob = guess_dist[k, ]
        )
      }
    } else {
      # Faster: generate all guesses at once and fill the sub-matrix
      n_total_guesses <- length(items_to_guess) * j
      guessed_values <- sample(
        seq_len(c), n_total_guesses, replace = TRUE, prob = guess_dist
      )
      ratings[items_to_guess, ] <- matrix(guessed_values, ncol = j)
    }
  }
  ratings
}

#' Calculate the True Latent Agreement (nu)
#'
#' @description
#' Calculates the true latent agreement coefficient, `nu`, for a variety of
#' skill-difficulty models. `nu` represents the population probability that a
#' randomly selected pair of judges both know the true category of a randomly
#' selected item.
#'
#' @param model A character string specifying the skill-difficulty model.
#'   Currently supported: `"jsm"` (Judge-Skill Model) and `"idm"`
#'   (Item-Difficulty Model).
#' @param model_params A list of parameters specific to the chosen model:
#'   - For `model = "jsm"`: `list(s = ...)` where `s` is a numeric vector
#'     of judge skill levels.
#'   - For `model = "idm"`: `list(beta_params = list(a = ..., b = ...))` where
#'     `a` and `b` are the shape parameters for the Beta distribution of item
#'     easiness.
#'
#' @return A single numeric value, the true latent agreement `nu`.
#' @export
#' @examples
#' # For 2 judges with skills 0.8 and 0.9
#' calculate_nu(model = "jsm", model_params = list(s = c(0.8, 0.9)))
#' # Expected: 0.8 * 0.9 = 0.72
#'
#' # For an item-difficulty model where easiness ~ Beta(5, 2)
#' calculate_nu(model = "idm", model_params = list(beta_params = list(a = 5, b = 2)))
calculate_nu <- function(model = c("jsm", "idm"), model_params = list()) {

  # --- 1. Argument validation (simplified for now) ---
  model <- match.arg(model)
  # A real implementation would have detailed checks on model_params

  # --- 2. Call the core calculation logic ---
  nu <- switch(
    model,
    "jsm" = nu_jsm(s = model_params$s),
    "idm" = nu_idm(beta_params = model_params$beta_params),
    # ... placeholders for future models
    stop("Model '", model, "' not yet implemented for nu calculation.")
  )

  return(nu)
}


#' Simulate Agreement Data from a Skill-Difficulty Model
#'
#' @description
#' Generates a matrix of simulated categorical ratings based on a specified
#' skill-difficulty model.
#'
#' @param n The number of items (subjects) to simulate.
#' @param j The number of judges (raters).
#' @param c The number of rating categories.
#' @param model A character string specifying the skill-difficulty model.
#'   See `calculate_nu` for details.
#' @param model_params A list of parameters for the chosen model. See `calculate_nu`.
#' @param true_dist A numeric vector of probabilities for the true category
#'   distribution. Must sum to 1. If NULL, a uniform distribution is assumed.
#' @param guessing_dist The distribution used when a judge guesses. Can be:
#'   - A numeric vector of probabilities (shared by all judges).
#'   - A matrix with `j` rows and `c` columns (for judge-specific guessing).
#'   If NULL, a uniform distribution is assumed.
#'
#' @return A `n x j` integer matrix of simulated ratings.
#' @export
#' @examples
#' # Simulate from a JSM with 3 judges
#' ratings_jsm <- simulate_agreement_data(
#'   n = 100, j = 3, c = 4,
#'   model = "jsm",
#'   model_params = list(s = c(0.9, 0.85, 0.7))
#' )
#'
#' # Simulate from an IDM where items are mostly easy (Beta(8, 2))
#' # and judges have judge-specific guessing biases
#' guess_biases <- rbind(c(0.5, 0.2, 0.2, 0.1), c(0.1, 0.5, 0.2, 0.2))
#' ratings_idm <- simulate_agreement_data(
#'   n = 50, j = 2, c = 4,
#'   model = "idm",
#'   model_params = list(beta_params = list(a = 8, b = 2)),
#'   guessing_dist = guess_biases
#' )
simulate_agreement_data <- function(n, j, c,
                                    model = c("jsm", "idm"),
                                    model_params = list(),
                                    true_dist = NULL,
                                    guessing_dist = NULL) {

  # --- 1. Argument validation (simplified for now) ---
  model <- match.arg(model)
  if (is.null(true_dist)) true_dist <- rep(1 / c, c)
  if (is.null(guessing_dist)) guessing_dist <- rep(1 / c, c)
  # ... detailed checks needed for lengths, sums, etc. ...

  # --- 2. Generate random components for the simulation ---
  sim_params <- switch(
    model,
    "jsm" = list(s = model_params$s),
    "idm" = list(delta_i = gen_idm(n, model_params$beta_params)),
    # ... placeholders for future models
    stop("Model '", model, "' not yet implemented for simulation.")
  )

  # --- 3. Call the core simulation logic ---
  ratings <- switch(
    model,
    "jsm" = sim_jsm(n, sim_params$s, true_dist, guessing_dist),
    "idm" = sim_idm(n, j, sim_params$delta_i, true_dist, guessing_dist)
  )

  # --- 4. Return formatted output ---
  # In the future, could add attributes like the true nu, model, etc.
  return(ratings)
}




# =====================================================================
# ==    BLACK-BOX VALIDATION TEST FOR AGREEMENT SIMULATOR            ==
# ==    (Compares against misskappa::kappa)                          ==
# =====================================================================

# We need the misskappa package for its kappa() function.
# Please ensure it is installed: install.packages("misskappa")
library(misskappa)

# --- Test Setup ---
# We'll use a large number of simulated items to ensure the empirical kappa
# is a precise estimate of its true population value.
n_large <- 100000
# We'll use a non-uniform distribution to make the test more rigorous.
# Under this condition, Fleiss' kappa would differ, but Conger's should match.
C <- 4
true_and_guessing_dist <- c(0.6, 0.2, 0.1, 0.1)
tolerance <- 1e-2 # Looser tolerance due to sampling variation in kappa


# =====================================================================
# ==                          TESTING JSM                          ==
# =====================================================================
cat("--- Testing Judge-Skill Model (JSM) vs. Conger's Kappa ---\n")

# 1. Define Population Parameters
jsm_params <- list(s = c(0.9, 0.85))
j_jsm <- length(jsm_params$s)

# 2. Calculate Theoretical Nu
# This is our ground truth for the population parameter.
nu_theoretical_jsm <- calculate_nu(model = "jsm", model_params = jsm_params)
cat("Theoretical nu (ground truth):", nu_theoretical_jsm, "\n")
# Expected: 0.9 * 0.85 = 0.765

# 3. Simulate Data Under the Key Assumption
# We generate a large dataset where the guessing distribution is identical
# to the true category distribution.
cat("Simulating JSM data (n=", n_large, ")... This may take a moment.\n")
ratings_jsm <- simulate_agreement_data(
  n = n_large,
  j = j_jsm,
  c = C,
  model = "jsm",
  model_params = jsm_params,
  true_dist = true_and_guessing_dist,
  guessing_dist = true_and_guessing_dist
)

# 4. Estimate Conger's Kappa from the Simulated Data
cat("Calculating Conger's kappa from simulated data...\n")
kappa_empirical_jsm <- kappa_raw(ratings_jsm)$estimates[1] # Using your package's function
cat("Empirical Conger's kappa:", kappa_empirical_jsm, "\n")

# 5. Check for Closeness
is_close_jsm <- abs(nu_theoretical_jsm - kappa_empirical_jsm) < tolerance
cat("--> Test passed:", is_close_jsm, "\n\n")


# =====================================================================
# ==                          TESTING IDM                          ==
# =====================================================================
cat("--- Testing Item-Difficulty Model (IDM) vs. Conger's Kappa ---\n")

# 1. Define Population Parameters
idm_params <- list(beta_params = list(a = 7, b = 2)) # Items are mostly easy
j_idm <- 3

# 2. Calculate Theoretical Nu
nu_theoretical_idm <- calculate_nu(model = "idm", model_params = idm_params)
cat("Theoretical nu (ground truth):", nu_theoretical_idm, "\n")
# Expected: (7 * 8) / (9 * 10) = 56 / 90 = 0.6222...

# 3. Simulate Data Under the Key Assumption
cat("Simulating IDM data (n=", n_large, ")... This may take a moment.\n")
ratings_idm <- simulate_agreement_data(
  n = n_large,
  j = j_idm,
  c = C,
  model = "idm",
  model_params = idm_params,
  true_dist = true_and_guessing_dist,
  guessing_dist = true_and_guessing_dist
)

# 4. Estimate Conger's Kappa from the Simulated Data
cat("Calculating Conger's kappa from simulated data...\n")
kappa_empirical_idm <- kappa_raw(ratings_idm)$estimates[1]
cat("Empirical Conger's kappa:", kappa_empirical_idm, "\n")

# 5. Check for Closeness
is_close_idm <- abs(nu_theoretical_idm - kappa_empirical_idm) < tolerance
cat("--> Test passed:", is_close_idm, "\n\n")
