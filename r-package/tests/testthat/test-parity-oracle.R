# Static parity oracle. The expected estimates and diagonal variances below
# were generated from the irrCAC-based `irrcacsmoke` oracle (complete-data
# Conger, Fleiss, and Brennan-Prediger coefficients, asymptotic Gamma / n
# variance convention) on 2026-06-13 and frozen here, so the regression runs
# without any external package. To refresh, install that oracle locally and
# recompute the constants from the same inputs.

test_that("available raw complete-data matches the irrCAC oracle", {
  x <- matrix(
    c(
      0, 0, 1,
      0, 1, 1,
      1, 1, 1,
      0, 0, 0,
      1, 0, 1,
      0, 1, 0,
      1, 1, 0,
      0, 0, 1,
      1, 1, 1,
      0, 0, 0
    ),
    ncol = 3,
    byrow = TRUE
  )
  storage.mode(x) <- "integer"

  fit <- estimate_kappa_raw(x, method = "available", weight = "identity")

  oracle_estimate <- c(
    Conger = 0.21052631578947364,
    Fleiss = 0.19999999999999996,
    `Brennan-Prediger` = 0.19999999999999996
  )
  oracle_variance <- c(
    Conger = 0.039752802694884175,
    Fleiss = 0.042666666666666679,
    `Brennan-Prediger` = 0.042666666666666679
  )

  nms <- names(oracle_estimate)
  expect_equal(unname(fit$estimates[nms]), unname(oracle_estimate),
               tolerance = 1e-9)
  expect_equal(unname(diag(fit$vcov)[nms]), unname(oracle_variance),
               tolerance = 1e-9)
})

test_that("counts complete-data matches the irrCAC oracle", {
  fit <- kappa_counts(dat.fleiss1971, weight = "nominal")

  oracle_estimate <- c(
    Fleiss = 0.43024452006014097,
    `Brennan-Prediger` = 0.44444444444444442
  )
  oracle_variance <- c(
    Fleiss = 0.0028396071239620256,
    `Brennan-Prediger` = 0.0029372427983539095
  )

  nms <- names(oracle_estimate)
  expect_equal(unname(fit$estimates[nms]), unname(oracle_estimate),
               tolerance = 1e-9)
  expect_equal(unname(diag(fit$vcov)[nms]), unname(oracle_variance),
               tolerance = 1e-9)
})
