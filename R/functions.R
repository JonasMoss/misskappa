### ===========================================================================
###
###   Functions used in the paper.
###
### ===========================================================================


#' Calculate the trace of a matrix.
#' @param mat A square matrix.
#' @return The trace of the matrix.
tr <- \(mat) sum(diag(mat))

#' Test the validity of the Perreault--Leigh model.
#'
#' This test uses a naïve approach, and it is probably possible to improve
#'   it significantly. Moreover, it only tests the validity of a model that
#'   is strictly less general than Gwet's model
#' @param tab A table of judgements.
#' @return A vector containing the statistic and associated p-value obtained
#'   from the chi squared test.
perreault_leigh_test <- \(tab) {
  n <- sum(tab)
  ff_hat <- tab / sum(tab)
  n_cat <- nrow(ff_hat)
  one_mat <- matrix(1, n_cat, n_cat)
  ones <- rep(1, n_cat)

  ff_mod_est <- \(p) {
    s <- p[1]
    t <- p[2:n_cat]
    t <- c(t, 1 - sum(t))
    middle <- (1 - s)^2 * one_mat / n_cat^2
    right <- s * (1 - s) * tcrossprod(t, ones) / n_cat +
      s * (1 - s) * tcrossprod(ones, t) / n_cat
    left <- s^2 * diag(t)
    left + right + middle
  }

  f <- \(p) {
    ff_mod <- ff_mod_est(p)
    n * sum((ff_mod - ff_hat)^2 / ff_mod)
  }

  ui <- rbind(
    c(1, rep(0, n_cat - 1)),
    c(-1, rep(0, n_cat - 1)),
    c(0, rep(1, n_cat - 1)),
    c(0, -rep(1, n_cat - 1)),
    cbind(rep(0, n_cat - 1), diag(n_cat - 1))
  )

  ci <- c(0, -1, 0, -1, rep(0, n_cat - 1))

  result <- constrOptim(
    theta = c(0.5, rep(1 / n_cat, n_cat - 1)),
    f = f,
    grad = NULL,
    ui = ui,
    ci = ci
  )

  stat <- result$value
  c(
    stat = result$value,
    p_value = pchisq(stat, n_cat^2 - n_cat - 2, lower.tail = FALSE)
  )
}

#' Test validity of the marginal guessing distribution model.
#'
#' @param tab A table of judgements.
#' @return A vector containing the statistic and associated p-value obtained
#'   from the chi squared test.
mgdm_test <- \(tab) {
  n <- sum(tab)
  ff_hat <- tab / sum(tab)
  f_hat <- (rowSums(ff_hat) + colSums(ff_hat)) / 2
  n_cat <- nrow(ff_hat)
  one_mat <- matrix(1, n_cat, n_cat)
  ones <- rep(1, n_cat)

  f <- \(nu) {
    right <- (1 - nu) * tcrossprod(f_hat, f_hat)
    left <- nu * diag(f_hat)
    ff_mod <- left + right
    sum((ff_mod - ff_hat)^2 / ff_mod)
  }

  stat <- n * optimize(f, c(0, 1))$objective

  c(
    stat = stat,
    p_value = pchisq(stat, n_cat^2 - n_cat - 1, lower.tail = FALSE)
  )
}
