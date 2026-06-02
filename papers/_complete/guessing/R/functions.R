#' Simulate distributions using the Dirichlet distribution
#'
#' @param true_alpha If `0`, the true distribution equals the uniform
#'   distribution when `type = "uniform"` and the mean of the guessing
#'   distributions if `type = "marginal"`. If `>0`, the true distribution is
#'   drawn from a Dirichlet with concentration parameter alpha with base measure
#'   as decided by `type`.
#' @param guessing_alpha If `type = "marginal"`, all guessing distributions are
#'   equal when `guessing_alpha = 0`; if not, they are interspersed according
#'   to a Dirichlet distribution. If `type = "uniform"`, the guessing
#'   distributions will be equal to the true distribution if
#'   `guessing_alpha = 0` and drawn from a Dirichlet with concentration
#'   `guessing_alpha` otherwise.
#' @param s Vector of expected skill-difficulty parameters.
#' @param type The type of study. Details in the description of `true_alpha` and
#'   `guessing_alpha`.
#' @param j,ratings Number of judges and number of ratings.
#' @internal
#' @return A list containing the true distribution and guessing distributions.

sim_distributions <- function(true_alpha = 1,
                              guessing_alpha = 1,
                              s,
                              type = c("uniform", "marginal"),
                              j = NULL,
                              ratings = NULL) {
  if (type == "marginal") {
    eps <- 0.00001
    a <- c(extraDistr::rdirichlet(1, 5 * rep(1, ratings)))
    # We use this mixture to avoid numerical problems.
    basis <- (1 - eps) * a + eps * rep(1, ratings)

    guessing_dist <- matrix(rep((1 - eps) * a + eps * rep(1, ratings), j), nrow = j, byrow = TRUE)
    if (guessing_alpha != 0) {
      a <- extraDistr::rdirichlet(j, guessing_alpha * basis)
      guessing_dist <- (1 - eps) * a + eps * guessing_dist
    }

    true_dist <- colSums(guessing_dist * (1 - s) / sum((1 - s)))
    if (true_alpha != 0) {
      true_dist <- c(extraDistr::rdirichlet(1, true_alpha * true_dist))
    }
  } else {
    basis <- rep(1, ratings)

    true_dist <- if (true_alpha == 0) {
      basis / ratings
    } else {
      c(extraDistr::rdirichlet(1, true_alpha * basis))
    }

    guessing_dist <- if (guessing_alpha == 0) {
      matrix(rep(basis / ratings, j), nrow = j, byrow = TRUE)
    } else {
      extraDistr::rdirichlet(j, guessing_alpha * basis)
    }
  }

  list(true_dist = true_dist, guessing_dist = guessing_dist)
}

#' Sensitivity of agreement coefficients to assumptions.
#'
#' @param true_alpha If `0`, the true distribution equals the uniform
#'   distribution when `type = "uniform"` and the mean of the guessing
#'   distributions if `type = "marginal"`. If `>0`, the true distribution is
#'   drawn from a Dirichlet with concentration parameter alpha with base measure
#'   as decided by `type`.
#' @param guessing_alpha If `type = "marginal"`, all guessing distributions are
#'   equal when `guessing_alpha = 0`; if not, they are interspersed according
#'   to a Dirichlet distribution. If `type = "uniform"`, the guessing
#'   distributions will be equal to the true distribution if
#'   `guessing_alpha = 0` and drawn from a Dirichlet with concentration
#'   `guessing_alpha` otherwise.
#' @param type The type of study. Details in the description of `true_alpha` and
#'   `guessing_alpha`.
#' @param model Choose between the judge skill model and the Gaussian copula
#'   beta model.
#' @param params A list of params. If the model is `"jsm"`,
#'   the params are `shape1` and `shape2`, used to sample `s` from a
#'   Beta distribution. If the model is `gcb`, the params are the
#'   correlaton `"\rho"` and params `shape1` and `shape2`.
#' @param j Function of zero arguments to sample number of judges from.
#'   Defaults to `sample(20, 1) + 2``.
#' @param ratings Function of zero arguments to sample number of categories
#'   from. Defaults to `sample(8, 1) + 2`.
#' @keywords internal
#' @param type The type of sensitivity study.
#' @param j The way to simulate judges.
#' @return A vector of squared errors.
sim_sensitivity <- function(true_alpha = 1,
                            guessing_alpha = 1,
                            type = c("uniform", "marginal"),
                            model = c("jsm", "gcb"),
                            params = NULL,
                            j = NULL,
                            ratings = NULL) {
  type <- match.arg(type)
  model <- match.arg(model)

  j <- if (is.null(j)) sample(20, 1) + 2 else j()
  ratings <- if (is.null(ratings)) sample(8, 1) + 2 else ratings()
  if (is.null(params)) {
    params <- list(
      shape1 = 7, shape2 = 1.5,
      rho = 0.5
    )
  }

  s <- if (model == "jsm") {
    stats::rbeta(j, params$shape1, params$shape2)
  } else if (model == "gcb") {
    mean <- params$shape1 / (params$shape1 + params$shape2)
    rep(mean, j)
  }

  dist <- sim_distributions(true_alpha, guessing_alpha, s, type, j, ratings)
  true_dist <- dist$true_dist
  guessing_dist <- dist$guessing_dist

  if (model == "jsm") {
    knowledge <- agreeable:::true_jsm(s, type = "knowledge", true_dist, guessing_dist)
    se <- function(x) abs(knowledge - x)

    c(
      "cf" = se(agreeable:::true_jsm(s, true_dist, guessing_dist, type = "cf")),
      "fleiss" = se(agreeable:::true_jsm(s, true_dist, guessing_dist, type = "fleiss")),
      "cohen" = se(agreeable:::true_jsm(s, true_dist, guessing_dist, type = "cohen")),
      "bp" = se(agreeable:::true_jsm(s, true_dist, guessing_dist, type = "bp")),
      "cbp" = se(agreeable:::true_jsm(s, true_dist, guessing_dist, type = "cbp"))
    )
  } else if (model == "gcb") {
    ss <- matrix(agreeable:::nbm_second_moment(
      params$rho,
      params$shape1,
      params$shape2,
      params$shape1,
      params$shape2
    ), ncol = j, nrow = j)

    knowledge <- agreeable:::true_knowledge(ss)
    se <- function(x) abs(knowledge - x)

    c(
      "cf" = se(agreeable:::true_cf(s, ss, true_dist, guessing_dist)),
      "fleiss" = se(agreeable:::true_fleiss(s, ss, true_dist, guessing_dist)),
      "cohen" = se(agreeable:::true_cohen(s, ss, true_dist, guessing_dist)),
      "bp" = se(agreeable:::true_bp(s, ss, true_dist, guessing_dist)),
      "cbp" = se(agreeable:::true_cbp(s, ss, true_dist, guessing_dist))
    )
  }
}

#' Simulate batch sensitivities for jsm
#'
#' @keywords internal
#' @param alphas Passed to `simulate_sensitivity`.
#' @param type Type of simulation study.
#' @param j How to sample judges.
#' @return Batch sensitivities.
sens_all_jsm <- function(alphas, n_reps, type = "uniform", j = NULL) {
  true_alpha <- alphas[1]
  guessing_alpha <- alphas[2]
  (rowMeans(replicate(n_reps, sim_sensitivity(
    true_alpha = true_alpha,
    guessing_alpha = guessing_alpha,
    type = type,
    j = j,
    model = "jsm",
    params = list(shape1 = 7, shape2 = 1.5)
  ))))
}

#' Simulate batch sensitivities for gcb.
#'
#' @keywords internal
#' @param alphas Passed to `simulate_sensitivity`.
#' @param type Type of simulation study.
#' @param j How to sample judges.
#' @return Batch sensitivities.
sens_all_gcb <- function(alphas, n_reps, type = "uniform", j = NULL,
                         rho = 0.5, shape1 = 7, shape2 = 1.5) {
  true_alpha <- alphas[1]
  guessing_alpha <- alphas[2]
  (rowMeans(replicate(n_reps, sim_sensitivity(
    true_alpha = true_alpha,
    guessing_alpha = guessing_alpha,
    type = type,
    j = j,
    model = "gcb",
    params = list(rho = rho, shape1 = shape1, shape2 = shape2)
  ))))
}

#' Simulate confidence interval of knowledge coefficients.
#'
#' @param true_alpha,guessing_alpha Passed to `sim_distributions`.
#' @param n Number of items rated.
#' @param n_reps Number of repetitions in the simulation.
#' @return A vector of simulated coverages and CI lengths.
sim_ci <- function(true_alpha, guessing_alpha, type = c("uniform", "marginal"), n, n_reps = 1000) {
  type = match.arg(type)
  models <- c("cf", "fleiss", "cohen", "bp", "cbp")
  in_ci <- function(ci, true_value) ci[1] <= true_value & ci[2] >= true_value
  ci_sim <- replicate(n_reps, {
    j <- rpois(1, 5) + 2
    ratings <- sample(8, 1) + 2
    s <- rbeta(j, 7, 1.5)
    dist <- agreeable:::sim_distributions(true_alpha, guessing_alpha, s, type, j, ratings)
    true_dist <- dist$true_dist
    guessing_dist <- dist$guessing_dist


    # Sometimes the data consists only of agreements; this is uninteresting,
    # so we toss them away.
    all_equal = TRUE
    while(all_equal) {
      y <- agreeable::simulate_jsm(
        n = n,
        s = s,
        true_dist = true_dist,
        guessing_dist = guessing_dist
      )

      all_equal = all(y[1] == y)

    }

    true_value <<- attr(y, "skill")
    cis <- lapply(models, function(model) {
      agreeable::knowledge(y, model = model)$conf.int
    })

    coverage <- sapply(cis, in_ci, true_value = true_value)
    length <- sapply(cis, function(ci) ci[2] - ci[1])
    c(coverage, length)

  })

  return_value <- rowMeans(ci_sim, na.rm = TRUE)
  return_value
}

#' Appropriately format the output of a simulation in the sensitivty analysis.
#'
#' @keywords internal
#' @param arr An array of simulations.
#' @return Formated simulations.
formater = function(arr) {
  arr <- t(arr)
  colnames(arr) <- NULL
  arr[arr < 4e-12] = 0
  arr <- format(c(arr), digits = 2, scientific = TRUE)
  dim(arr) <- c(9, 5)
  arr
}

# This function is passed to `apply.`
sim_helper = function(x, type, n_reps) {
  sim_ci(true_alpha = x[2], guessing_alpha = x[1], n = x[3], type = type, n_reps = n_reps)
}
