#' Simulation helpers for missing-data agreement studies
#'
#' @description
#' `sim` is a single exported object that bundles several simulation
#' generators useful for testing the package's estimators under MCAR / MAR
#' missingness. Exposed as a list of closures so the manual is not flooded
#' by individual function pages.
#'
#' Available members:
#'
#' - `sim$mcar(n, R, C, p, p_missing, seed = NULL)`: simulate an n x R
#'   matrix of integer ratings (categories `1..C`) with marginal
#'   distribution `p` (length C), then drop entries independently per
#'   cell with probability `p_missing` (scalar or length R). Returns the
#'   incomplete matrix with the underlying complete matrix as the
#'   `"complete"` attribute.
#'
#' - `sim$mar_truth(n, R, C, p, pi_truth, seed = NULL)`: like `mcar` but
#'   the per-cell missingness probability depends on the unobserved truth
#'   category; `pi_truth` is a length-C vector of missing probabilities,
#'   one per category. MAR via the underlying truth.
#'
#' - `sim$jsm(n, s, model, true_dist, guessing_dist)`: port of the legacy
#'   `simulate_jsm()` skill-difficulty guessing simulator. See
#'   the legacy docs for details. Returns the simulated rating matrix with
#'   `"kappa"` (true latent agreement) attached.
#'
#' All generators accept a `seed` argument; if non-NULL it is passed to
#' `set.seed()` before sampling.
#'
#' @examples
#' # Simulate MCAR categorical ratings (200 subjects, 3 raters, 4 categories,
#' # 15% missing per cell) and estimate kappa on the incomplete matrix.
#' x <- sim$mcar(n = 200, R = 3, C = 4, p = c(0.4, 0.3, 0.2, 0.1),
#'               p_missing = 0.15, seed = 1)
#' kappa(x, estimator = "ipw")
#'
#' # The underlying complete matrix is kept as an attribute for comparison.
#' dim(attr(x, "complete"))
#'
#' @keywords datasets
#' @export
sim <- local({

  set_seed_if <- function(seed) {
    if (!is.null(seed)) set.seed(seed)
  }

  mcar_impl <- function(n, R, C, p, p_missing, seed = NULL) {
    set_seed_if(seed)
    if (length(p) != C) stop("'p' must have length C.")
    if (length(p_missing) == 1L) p_missing <- rep(p_missing, R)
    if (length(p_missing) != R) stop("'p_missing' must be scalar or length R.")
    x_star <- matrix(sample.int(C, n * R, replace = TRUE, prob = p),
                     nrow = n, ncol = R)
    x <- x_star
    for (j in seq_len(R)) {
      drop <- stats::rbinom(n, 1, p_missing[j]) == 1L
      x[drop, j] <- NA_integer_
    }
    attr(x, "complete") <- x_star
    x
  }

  mar_truth_impl <- function(n, R, C, p, pi_truth, seed = NULL) {
    set_seed_if(seed)
    if (length(p) != C) stop("'p' must have length C.")
    if (length(pi_truth) != C) stop("'pi_truth' must have length C.")
    truth <- sample.int(C, n, replace = TRUE, prob = p)
    x_star <- matrix(NA_integer_, nrow = n, ncol = R)
    for (j in seq_len(R)) {
      x_star[, j] <- truth
    }
    x <- x_star
    for (i in seq_len(n)) {
      pmiss <- pi_truth[truth[i]]
      drop <- stats::rbinom(R, 1, pmiss) == 1L
      x[i, drop] <- NA_integer_
    }
    attr(x, "complete") <- x_star
    attr(x, "truth") <- truth
    x
  }

  jsm_impl <- function(n, s,
                       model = c("general", "cohen-fleiss", "fleiss", "bp", "tu"),
                       true_dist = NULL, guessing_dist = NULL, seed = NULL) {
    set_seed_if(seed)
    model <- match.arg(model)
    J <- length(s)

    if (!is.null(guessing_dist)) {
      if (is.matrix(guessing_dist)) {
        guessing_dist <- guessing_dist / rowSums(guessing_dist)
      } else {
        guessing_dist <- guessing_dist / sum(guessing_dist)
      }
    }
    if (!is.null(true_dist)) true_dist <- true_dist / sum(true_dist)

    q <- if (!is.null(true_dist)) length(true_dist)
         else if (!is.null(guessing_dist))
           if (is.matrix(guessing_dist)) ncol(guessing_dist) else length(guessing_dist)
         else stop("Either true_dist or guessing_dist must be supplied.")

    if (model == "cohen-fleiss") {
      if (is.null(guessing_dist)) stop("guessing_dist required for 'cohen-fleiss'.")
      if (is.null(true_dist)) {
        if (is.matrix(guessing_dist)) {
          w <- (1 - s) / sum(1 - s)
          true_dist <- colSums(guessing_dist * w)
        } else {
          true_dist <- guessing_dist
        }
      }
    } else if (model == "tu") {
      true_dist <- rep(1 / q, q)
    } else if (model == "bp") {
      if (is.null(true_dist)) stop("true_dist required for 'bp'.")
      guessing_dist <- rep(1 / q, q)
    } else if (model == "fleiss") {
      if (is.null(true_dist)) stop("true_dist required for 'fleiss'.")
      guessing_dist <- true_dist
    }

    x_star <- sample.int(q, n, replace = TRUE, prob = true_dist)
    obs <- if (is.matrix(guessing_dist)) {
      sapply(seq_along(s), function(j) {
        z <- stats::rbinom(n, 1, s[j])
        z * x_star + (1 - z) * sample.int(q, n, replace = TRUE, prob = guessing_dist[j, ])
      })
    } else {
      sapply(seq_along(s), function(j) {
        z <- stats::rbinom(n, 1, s[j])
        z * x_star + (1 - z) * sample.int(q, n, replace = TRUE, prob = guessing_dist)
      })
    }

    ss <- s %*% t(s)
    diag(ss) <- 0
    attr(obs, "kappa") <- sum(ss) / (J * (J - 1))
    attr(obs, "truth") <- x_star
    obs
  }

  list(
    mcar = mcar_impl,
    mar_truth = mar_truth_impl,
    jsm = jsm_impl
  )
})
