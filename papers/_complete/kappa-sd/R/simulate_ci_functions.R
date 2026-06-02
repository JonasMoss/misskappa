#' Simulate CIs for Cohen's kappa and Fleiss' kappa in a guessing model
#' @param s,true_dist,model Parameters passe to `simulate_jsm`,
#' @param n Number of items rated.
#' @param type Type of agreement coefficient calculated.
#' @param weight The weight passed to `cohen_kappa` and `fleiss_kappa`.
#' @param n_reps Number of repetitions in the simulation.
#' @return A vector of simulated coverages and CI lengths.
sim_ci_guessing_model = function(
  s,
  true_dist,
  model = "bp",
  n,
  type = c("cohen", "fleiss"),
  weight = "unit",
  n_reps = 1000) {

  type = match.arg(type)

  in_ci = function(ci, true_value) ci[1] <= true_value & ci[2] >= true_value

  ci_sim = replicate(n_reps, {

    y = agreeable::simulate_jsm(
      n = n,
      s = s,
      true_dist = true_dist,
      model = model)

    true_value <<- attr(y, "skill")

    if(type == "cohen") {
      basic_ci = agreeable::cohen_kappa(y, method = "basic", weight = weight)$conf.int
      arcsine_ci = agreeable::cohen_kappa(y, method = "arcsine", weight = weight)$conf.int
      fisher_ci = agreeable::cohen_kappa(y, method = "fisher", weight = weight)$conf.int
    } else {
      basic_ci = agreeable::fleiss_kappa(y, method = "basic", weight = weight)$conf.int
      arcsine_ci = agreeable::fleiss_kappa(y, method = "arcsine", weight = weight)$conf.int
      fisher_ci = agreeable::fleiss_kappa(y, method = "fisher", weight = weight)$conf.int
    }

    c("basic_coverage" = in_ci(basic_ci, true_value),
      "arcsine_coverage" = in_ci(arcsine_ci, true_value),
      "fisher_coverage" = in_ci(fisher_ci,true_value),
      "basic_length" = basic_ci[2] - basic_ci[1],
      "arcsine_length" = arcsine_ci[2] - arcsine_ci[1],
      "fisher_length" = fisher_ci[2] - fisher_ci[1])

  })

  return_value = rowMeans(ci_sim, na.rm = TRUE)
  attr(return_value, "true_value") = true_value
  return_value
}


#' Simulate CIs for Cohen's kappa and Fleiss' kappa in a normal model
#' @param rho Common correlation.
#' @param j Number of judges.
#' @param n Number of items rated.
#' @param type Type of agreement coefficient calculated.
#' @param weight The weight passed to `cohen_kappa` and `fleiss_kappa`.
#' @param n_reps Number of repetitions in the simulation.
#' @return A vector of simulated coverages and CI lengths.
sim_ci_normal_model = function(
  rho = 0.7,
  j,
  n,
  type = c("cohen", "fleiss"),
  weight = "unit",
  n_reps = 1000) {

  type = match.arg(type)

  Sigma = matrix(rep(rho, j * j), nrow = j)
  diag(Sigma) = 1
  if(weight == "quadratic") {
    true_value = rho
  } else if (weight == "abs") {
    true_value = 1 - sqrt(1 - rho)
  }

  in_ci = function(ci, true_value) ci[1] <= true_value & ci[2] >= true_value

  ci_sim = replicate(n_reps, {

    y = MASS::mvrnorm(
      n = n,
      mu = rep(0, j),
      Sigma = Sigma)

    if(type == "cohen") {
      basic_ci = agreeable::cohen_kappa(y, method = "basic", weight = weight)$conf.int
      arcsine_ci = agreeable::cohen_kappa(y, method = "arcsine", weight = weight)$conf.int
      fisher_ci = agreeable::cohen_kappa(y, method = "fisher", weight = weight)$conf.int
    } else {
      basic_ci = agreeable::fleiss_kappa(y, method = "basic", weight = weight)$conf.int
      arcsine_ci = agreeable::fleiss_kappa(y, method = "arcsine", weight = weight)$conf.int
      fisher_ci = agreeable::fleiss_kappa(y, method = "fisher", weight = weight)$conf.int
    }

    c("basic_coverage" = in_ci(basic_ci, true_value),
      "arcsine_coverage" = in_ci(arcsine_ci, true_value),
      "fisher_coverage" = in_ci(fisher_ci,true_value),
      "basic_length" = basic_ci[2] - basic_ci[1],
      "arcsine_length" = arcsine_ci[2] - arcsine_ci[1],
      "fisher_length" = fisher_ci[2] - fisher_ci[1])

  })

  return_value = rowMeans(ci_sim, na.rm = TRUE)
  attr(return_value, "true_value") = true_value
  return_value
}