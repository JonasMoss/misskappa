#' Simulate CIs for Cohen's kappa and Fleiss' kappa in a guessing model
#' @param s,true_dist,model Parameters passe to `simulate_jsm`,
#' @param n Number of items rated.
#' @param type Type of agreement coefficient calculated.
#' @param disagreement The disagreement function passed to `cohen_kappa` and
#'    `fleiss_kappa`.
#' @param n_reps Number of repetitions in the simulation.
#' @return A vector of simulated coverages and CI lengths.
sim_ci_guessing_model_g_wise <- function(s,
                                         true_dist,
                                         model = "bp",
                                         n,
                                         type = c("cohen", "fleiss"),
                                         disagreement = "nominal",
                                         g,
                                         true_value,
                                         n_reps = 1000,
                                         seed = 1) {
  type <- match.arg(type)
  set.seed(seed)
  in_ci <- \(ci, true_value) ci[1] <= true_value & ci[2] >= true_value

  ci_sim <- replicate(n_reps, {
    y <- agreeable::simulate_jsm(
      n = n,
      s = s,
      true_dist = true_dist,
      model = model
    )

    arcsine_ci <- agreer::kappa(y, method = "arcsine", disagreement = disagreement, type = type, g = g)$conf.int
    fisher_ci <- agreer::kappa(y, method = "fisher", disagreement = disagreement, type = type, g = g)$conf.int

    c(
      "arcsine_coverage" = in_ci(arcsine_ci, true_value),
      "fisher_coverage" = in_ci(fisher_ci, true_value),
      "arcsine_length" = arcsine_ci[2] - arcsine_ci[1],
      "fisher_length" = fisher_ci[2] - fisher_ci[1]
    )
  })

  return_value <- rowMeans(ci_sim, na.rm = TRUE)
  attr(return_value, "true_value") <- true_value
  return_value
}
