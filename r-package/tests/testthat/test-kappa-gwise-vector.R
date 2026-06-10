test_that("kappa() g = 2 is unchanged by passing g explicitly", {
  x <- matrix(c(
    1, 1, 2, 1,
    1, 2, 2, 2,
    2, 1, 1, 1,
    2, 3, 3, 3,
    1, 1, 1, 2,
    3, 3, 2, 3
  ), nrow = 6, byrow = TRUE)

  default_fit <- kappa(x, estimator = "ipw", weight = "nominal")
  g2_fit <- kappa(x, estimator = "ipw", weight = "nominal", g = 2L)
  expect_equal(g2_fit$estimates, default_fit$estimates)
  expect_equal(g2_fit$vcov, default_fit$vcov)
})

test_that("kappa() routes g > 2 categorical kernels and labels Conger/Fleiss", {
  x <- matrix(c(
    1, 1, 2, 1, 1,
    1, 2, 3, 2, 2,
    2, 1, 1, 1, 1,
    2, 3, 4, 4, 5,
    1, 1, 1, 2, 1,
    3, 3, 2, 3, 3
  ), nrow = 6, byrow = TRUE)

  fit <- kappa(x, estimator = "pairwise", weight = "nominal", g = 3L)
  expect_s3_class(fit, "misskappa_estimate")
  expect_named(fit$estimates, c("Conger", "Fleiss"))
  expect_equal(fit$g, 3L)
  expect_true(all(is.finite(fit$estimates)))

  # estimator = "pairwise" with g > 2 is the complete-data g-wise estimator.
  direct <- kappa_gwise(x, distance = "nominal", method = "complete", g = 3L)
  expect_equal(unname(fit$estimates), unname(direct$estimates))
  expect_equal(unname(fit$vcov), unname(direct$vcov))

  # IPW on complete data matches the complete-data point estimates.
  fit_ipw <- kappa(x, estimator = "ipw", weight = "nominal", g = 3L)
  expect_equal(unname(fit_ipw$estimates), unname(fit$estimates),
               tolerance = 1e-8)

  # cat_fiml runs and yields finite estimates.
  fit_fiml <- kappa(x, estimator = "cat_fiml", weight = "nominal", g = 3L)
  expect_named(fit_fiml$estimates, c("Conger", "Fleiss"))
  expect_true(all(is.finite(fit_fiml$estimates)))

  # hubert kernel (g > 2 only).
  fit_hub <- kappa(x, estimator = "ipw", weight = "hubert", g = 3L)
  expect_true(all(is.finite(fit_hub$estimates)))
})

test_that("kappa() g > 2 with NAs: ipw/cat_fiml work, pairwise errors", {
  complete_support <- as.matrix(expand.grid(
    r1 = 1:2, r2 = 1:2, r3 = 1:2, r4 = 1:2, r5 = 1:2
  ))
  missing_rows <- matrix(c(
    1, NA, 2, 1, 1,
    2, 1, NA, 2, 2,
    1, 2, 1, NA, 1,
    2, 2, 2, 1, NA
  ), ncol = 5, byrow = TRUE)
  x <- rbind(complete_support, missing_rows)

  fit_ipw <- kappa(x, estimator = "ipw", weight = "nominal", g = 3L)
  expect_true(all(is.finite(fit_ipw$estimates)))

  fit_fiml <- kappa(x, estimator = "cat_fiml", weight = "nominal", g = 3L)
  expect_true(all(is.finite(fit_fiml$estimates)))

  expect_error(
    kappa(x, estimator = "pairwise", weight = "nominal", g = 3L),
    "complete"
  )
})

test_that("kappa() g > 2 continuous (linear) supports pairwise/ipw only", {
  x <- matrix(c(
    1, 1, 2, 1, 1,
    1, 2, 3, 2, 2,
    2, 1, 1, 1, 1,
    2, 3, 4, 4, 5,
    1, 1, 1, 2, 1,
    3, 3, 2, 3, 3
  ), nrow = 6, byrow = TRUE)

  fit_pair <- kappa(x, estimator = "pairwise", weight = "linear", g = 3L)
  expect_named(fit_pair$estimates, c("Conger", "Fleiss"))
  expect_true(all(is.finite(fit_pair$estimates)))

  fit_ipw <- kappa(x, estimator = "ipw", weight = "linear", g = 3L)
  expect_true(all(is.finite(fit_ipw$estimates)))

  expect_error(
    kappa(x, estimator = "cat_fiml", weight = "linear", g = 3L),
    "continuous g-wise"
  )
  expect_error(
    kappa(x, estimator = "nt_fiml", weight = "linear", g = 3L),
    "quadratic"
  )
})

test_that("kappa() quadratic weighting is g-invariant (cheap path)", {
  x <- matrix(c(
    1, 1, 2, 1, 1,
    1, 2, 3, 2, 2,
    2, 1, 1, 1, 1,
    2, 3, 4, 4, 5,
    1, 1, 1, 2, 1,
    3, 3, 2, 3, 3
  ), nrow = 6, byrow = TRUE)

  base <- kappa(x, estimator = "pairwise", weight = "quadratic")
  for (gg in c(2L, 3L, 5L)) {
    fit <- kappa(x, estimator = "pairwise", weight = "quadratic", g = gg)
    expect_equal(fit$estimates, base$estimates)
    expect_equal(fit$vcov, base$vcov)
  }

  # The same g-invariance holds for the categorical IPW quadratic path.
  base_ipw <- kappa(x, estimator = "ipw", weight = "quadratic")
  fit_ipw <- kappa(x, estimator = "ipw", weight = "quadratic", g = 4L)
  expect_equal(fit_ipw$estimates, base_ipw$estimates)
})

test_that("kappa() errors clearly on unsupported g > 2 combinations", {
  x <- matrix(c(
    1, 1, 2, 1, 1,
    1, 2, 3, 2, 2,
    2, 1, 1, 1, 1,
    2, 3, 4, 4, 5
  ), nrow = 4, byrow = TRUE)

  expect_error(kappa(x, estimator = "ipw", weight = "ordinal", g = 3L),
               "should be one of")
  expect_error(kappa(x, estimator = "ipw", weight = "hubert", g = 2L),
               "g > 2")
})

test_that("kappa() detects 3-D arrays as vector-valued ratings", {
  set.seed(1)
  n <- 8L
  R <- 3L
  p <- 2L
  arr <- array(sample(0:1, n * R * p, replace = TRUE), dim = c(n, R, p))
  storage.mode(arr) <- "double"

  fit <- kappa(arr, estimator = "pairwise")
  expect_s3_class(fit, "misskappa_estimate")
  expect_named(fit$estimates, c("Conger", "Fleiss"))

  # default weight = "nominal" maps to the Hamming component loss.
  direct <- kappa_vector(arr, method = "pairwise", loss = "hamming")
  expect_equal(unname(fit$estimates), unname(direct$estimates))

  # weight -> loss mapping for the squared component loss.
  fit_sq <- kappa(arr, estimator = "ipw", weight = "quadratic")
  direct_sq <- kappa_vector(arr, method = "ipw", loss = "squared")
  expect_equal(unname(fit_sq$estimates), unname(direct_sq$estimates))

  # rms / feature_weights reachable through ...
  fit_rms <- kappa(arr, estimator = "ipw", loss = "rms",
                   feature_weights = c(1, 2))
  direct_rms <- kappa_vector(arr, method = "ipw", loss = "rms",
                             feature_weights = c(1, 2))
  expect_equal(unname(fit_rms$estimates), unname(direct_rms$estimates))

  expect_error(kappa(arr, estimator = "cat_fiml"), "pairwise")
  expect_error(kappa(arr, estimator = "nt_fiml"), "pairwise")
})
