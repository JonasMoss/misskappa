#' Simulate the true value of a judge skill model.
#'
#' The true values of a judge skill model are easily calculate for the
#'    pair-wise coefficients, but not the g-wise. We used this function to
#'    approximate the true values, used the in the confidence interval
#'    simulation. Since all judges have the same skill, the population value
#'    of Fleiss' kappa and Cohen's kappa are equal, so we may use the much
#'    easier to compute Fleiss' chance agreement.
#'
#' @param n Number of items used in simulation.
#' @param r Number of raters.
#' @param true_dist True distribution.
#' @param s Judge skills.
#' @param model Model, passed to `agreeable::jsm`.
#' @param g The arity of the disagreement.
#' @param disagreement The disagreement function.
#' @return Monte Carlo estimate of the value of the kappa
sim_jsm_true_value <- function(n, r, true_dist, s, model, g, disagreement) {
  if (g == 2) {
    return(attr(agreeable::simulate_jsm(1, s, model, true_dist), "skill"))
  }

  df <- eval(parse(text = paste0("agreer:::", disagreement))[[1]])
  disagreements <- c("nominal", "absolute", "quadratic", "hubert")
  ds <- which(disagreement == disagreements) - 1

  m <- length(true_dist)
  with_item <- agreeable::simulate_jsm(n, s, model, true_dist)
  without_item <- matrix(sample(seq(m), size = n * g, replace = TRUE), nrow = n)

  d <- mean(agreer:::muds_simple_cpp(with_item, g, ds))
  f <- mean(agreer:::muds_simple_cpp(without_item, g, ds))
  1 - d / f
}
