test_that("agreement() g = 2 is unchanged by passing g explicitly", {
  x <- matrix(c(
    1, 1, 2, 1,
    1, 2, 2, 2,
    2, 1, 1, 1,
    2, 3, 3, 3,
    1, 1, 1, 2,
    3, 3, 2, 3
  ), nrow = 6, byrow = TRUE)

  default_fit <- agreement(x, estimator = "ipw", weight = "nominal")
  g2_fit <- agreement(x, estimator = "ipw", weight = "nominal", g = 2L)
  expect_equal(g2_fit$estimates, default_fit$estimates)
  expect_equal(g2_fit$vcov, default_fit$vcov)
})

test_that("agreement() routes g > 2 categorical kernels and labels Conger/Fleiss", {
  x <- matrix(c(
    1, 1, 2, 1, 1,
    1, 2, 3, 2, 2,
    2, 1, 1, 1, 1,
    2, 3, 4, 4, 5,
    1, 1, 1, 2, 1,
    3, 3, 2, 3, 3
  ), nrow = 6, byrow = TRUE)

  fit <- agreement(x, estimator = "ipw", weight = "nominal", g = 3L)
  expect_s3_class(fit, "misskappa_estimate")
  expect_named(fit$estimates, c("Conger", "Fleiss"))
  expect_equal(fit$g, 3L)
  expect_true(all(is.finite(fit$estimates)))

  # On complete data IPW is the complete-data g-wise estimator, which is why
  # agreement() no longer offers "pairwise" for g > 2.
  direct <- research_kappa_gwise(x, distance = "nominal", method = "complete", g = 3L)
  expect_equal(unname(fit$estimates), unname(direct$estimates),
               tolerance = 1e-12)
  expect_equal(unname(fit$vcov), unname(direct$vcov), tolerance = 1e-12)
  expect_error(agreement(x, estimator = "pairwise", weight = "nominal", g = 3L),
               "quadratic estimator")

  # cat_fiml runs and yields finite estimates.
  fit_fiml <- agreement(x, estimator = "cat_fiml", weight = "nominal", g = 3L)
  expect_named(fit_fiml$estimates, c("Conger", "Fleiss"))
  expect_true(all(is.finite(fit_fiml$estimates)))

  # hubert kernel (g > 2 only).
  fit_hub <- agreement(x, estimator = "ipw", weight = "hubert", g = 3L)
  expect_true(all(is.finite(fit_hub$estimates)))
})

test_that("agreement() g > 2 with NAs: ipw/cat_fiml work", {
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

  fit_ipw <- agreement(x, estimator = "ipw", weight = "nominal", g = 3L)
  expect_true(all(is.finite(fit_ipw$estimates)))

  fit_fiml <- agreement(x, estimator = "cat_fiml", weight = "nominal", g = 3L)
  expect_true(all(is.finite(fit_fiml$estimates)))
})

test_that("agreement() g > 2 continuous (linear) supports ipw only", {
  x <- matrix(c(
    1, 1, 2, 1, 1,
    1, 2, 3, 2, 2,
    2, 1, 1, 1, 1,
    2, 3, 4, 4, 5,
    1, 1, 1, 2, 1,
    3, 3, 2, 3, 3
  ), nrow = 6, byrow = TRUE)

  fit_ipw <- agreement(x, estimator = "ipw", weight = "linear", g = 3L)
  expect_named(fit_ipw$estimates, c("Conger", "Fleiss"))
  expect_true(all(is.finite(fit_ipw$estimates)))
  direct <- research_kappa_gwise(x, distance = "absolute", method = "complete", g = 3L)
  expect_equal(unname(fit_ipw$estimates), unname(direct$estimates),
               tolerance = 1e-10)
  expect_error(agreement(x, estimator = "pairwise", weight = "linear", g = 3L),
               "quadratic estimator")

  expect_error(
    agreement(x, estimator = "cat_fiml", weight = "linear", g = 3L),
    "continuous g-wise"
  )
  expect_error(
    agreement(x, estimator = "nt_fiml", weight = "linear", g = 3L),
    "quadratic"
  )
})

test_that("agreement() quadratic weighting is g-invariant (cheap path)", {
  x <- matrix(c(
    1, 1, 2, 1, 1,
    1, 2, 3, 2, 2,
    2, 1, 1, 1, 1,
    2, 3, 4, 4, 5,
    1, 1, 1, 2, 1,
    3, 3, 2, 3, 3
  ), nrow = 6, byrow = TRUE)

  base <- agreement(x, estimator = "pairwise", weight = "quadratic")
  for (gg in c(2L, 3L, 5L)) {
    fit <- agreement(x, estimator = "pairwise", weight = "quadratic", g = gg)
    expect_equal(fit$estimates, base$estimates)
    expect_equal(fit$vcov, base$vcov)
  }

  # The same g-invariance holds for the categorical IPW quadratic path.
  base_ipw <- agreement(x, estimator = "ipw", weight = "quadratic")
  fit_ipw <- agreement(x, estimator = "ipw", weight = "quadratic", g = 4L)
  expect_equal(fit_ipw$estimates, base_ipw$estimates)
})

test_that("agreement() errors clearly on unsupported g > 2 combinations", {
  x <- matrix(c(
    1, 1, 2, 1, 1,
    1, 2, 3, 2, 2,
    2, 1, 1, 1, 1,
    2, 3, 4, 4, 5
  ), nrow = 4, byrow = TRUE)

  expect_error(agreement(x, estimator = "ipw", weight = "ordinal", g = 3L),
               "should be one of")
  expect_error(agreement(x, estimator = "ipw", weight = "hubert", g = 2L),
               "g > 2")
})

test_that("agreement() detects 3-D arrays as vector-valued ratings", {
  set.seed(1)
  n <- 8L
  R <- 3L
  p <- 2L
  arr <- array(sample(0:1, n * R * p, replace = TRUE), dim = c(n, R, p))
  storage.mode(arr) <- "double"

  fit <- agreement(arr, estimator = "pairwise")
  expect_s3_class(fit, "misskappa_estimate")
  expect_named(fit$estimates, c("Conger", "Fleiss"))

  # default weight = "nominal" maps to the Hamming component loss.
  direct <- .kappa_vector(arr, method = "pairwise", loss = "hamming")
  expect_equal(unname(fit$estimates), unname(direct$estimates))

  # weight -> loss mapping for the squared component loss.
  fit_sq <- agreement(arr, estimator = "ipw", weight = "quadratic")
  direct_sq <- .kappa_vector(arr, method = "ipw", loss = "squared")
  expect_equal(unname(fit_sq$estimates), unname(direct_sq$estimates))

  # rms / feature_weights reachable through ...
  fit_rms <- agreement(arr, estimator = "ipw", loss = "rms",
                   feature_weights = c(1, 2))
  direct_rms <- .kappa_vector(arr, method = "ipw", loss = "rms",
                             feature_weights = c(1, 2))
  expect_equal(unname(fit_rms$estimates), unname(direct_rms$estimates))

  expect_error(agreement(arr, estimator = "cat_fiml"), "pairwise")
  # nt_fiml / nt_fiml_cs are supported on arrays (full-W quadratic FIML,
  # default identity metric); Brennan-Prediger needs a declared `ranges`.
  fit_ml <- agreement(arr, estimator = "nt_fiml")
  expect_equal(names(coef(fit_ml)), c("Conger", "Fleiss"))
  expect_equal(fit_ml$method, "nt_fiml")
  rng <- cbind(lo = apply(arr, 3L, min, na.rm = TRUE),
               hi = apply(arr, 3L, max, na.rm = TRUE))
  fit_bp <- agreement(arr, estimator = "nt_fiml", ranges = rng)
  expect_equal(names(coef(fit_bp)), c("Conger", "Fleiss", "Brennan-Prediger"))
  expect_equal(coef(fit_bp)[1:2], coef(fit_ml))
})
