test_that("kappa() routes the four methods and returns the S3 shape", {
  x <- matrix(
    c(0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1),
    nrow = 10, byrow = TRUE
  )
  storage.mode(x) <- "integer"

  for (m in c("available", "ipw", "fiml", "gwet")) {
    fit <- kappa(x, method = m, weight = "identity")
    expect_s3_class(fit, "misskappa_estimate")
    expect_named(fit$estimates, c("Conger", "Fleiss", "Brennan-Prediger"))
    expect_equal(dim(fit$vcov), c(3L, 3L))
    expect_true(all(is.finite(fit$estimates)))
  }
})

test_that("available-case on a Cohen-style 2-rater example matches the textbook value", {
  x <- matrix(
    c(0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1),
    nrow = 10, byrow = TRUE
  )
  fit <- kappa(x, method = "available", weight = "identity")
  # Cohen's kappa for this contingency table is 0.4.
  expect_equal(unname(fit$estimates["Conger"]), 0.4, tolerance = 1e-9)
  expect_equal(unname(fit$estimates["Brennan-Prediger"]), 0.4, tolerance = 1e-9)
})

test_that("S3 methods on misskappa_estimate behave", {
  x <- matrix(c(0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1), nrow = 6, byrow = TRUE)
  fit <- kappa(x, method = "available", weight = "identity")

  est <- coef(fit)
  expect_named(est, c("Conger", "Fleiss", "Brennan-Prediger"))

  V <- vcov(fit)
  expect_equal(dim(V), c(3L, 3L))
  expect_true(isSymmetric(V))

  ci <- confint(fit)
  expect_equal(dim(ci), c(3L, 2L))
  expect_true(all(ci[, 1] <= ci[, 2]))

  df <- as.data.frame(fit)
  expect_equal(nrow(df), 3L)
  expect_named(df, c("coefficient", "estimate", "se"))

  expect_output(print(fit), "method=available")
})

test_that("IPW rejects a rater with zero observations cleanly", {
  x <- matrix(c(0, 0, NA, 1, 1, NA, 0, 1, NA, 1, 0, NA, 0, 0, NA),
              nrow = 5, byrow = TRUE)
  expect_error(kappa(x, method = "ipw", weight = "identity"),
               regexp = "singular")
})

test_that("FIML matches available on complete data, identity weight", {
  x <- matrix(
    c(0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1),
    nrow = 10, byrow = TRUE
  )
  fit_av <- kappa(x, method = "available", weight = "identity")
  fit_ml <- kappa(x, method = "fiml", weight = "identity")
  expect_equal(unname(fit_ml$estimates), unname(fit_av$estimates),
               tolerance = 1e-6)
})

test_that("sim$mcar produces an n x R matrix with the right missingness rate", {
  set.seed(7)
  x <- sim$mcar(n = 500, R = 4, C = 3, p = c(0.4, 0.4, 0.2),
                p_missing = 0.2, seed = 7)
  expect_equal(dim(x), c(500L, 4L))
  rate <- mean(is.na(x))
  expect_true(abs(rate - 0.2) < 0.05)
  comp <- attr(x, "complete")
  expect_equal(dim(comp), c(500L, 4L))
  expect_true(!any(is.na(comp)))
})

test_that("sim$jsm port still produces sensible output", {
  set.seed(1)
  x <- sim$jsm(n = 100, s = c(0.7, 0.5, 0.8),
               model = "fleiss", true_dist = c(0.5, 0.3, 0.2), seed = 1)
  expect_equal(dim(x), c(100L, 3L))
  expect_true(!is.null(attr(x, "kappa")))
  expect_true(attr(x, "kappa") > 0 && attr(x, "kappa") < 1)
})

test_that("Bundled datasets load", {
  expect_true(is.numeric(as.matrix(dat.fleiss1971)))
  expect_true(is.numeric(as.matrix(dat.gwet2014)))
  expect_true(is.numeric(as.matrix(dat.klein2018)))
  expect_true(is.numeric(as.matrix(dat.zapf2016)))
})

test_that("kappa_continuous() handles available/ipw/gwet and matches legacy", {
  x <- matrix(c(
    1.0, 1.2, 0.9,  3.5, 3.4, 3.6,  2.0, 2.1, 1.9,  4.0, 4.2, 3.8,
    1.5, 1.4, 1.6,  3.0, 2.9, 3.1,  2.5, 2.6, 2.4,  4.5, 4.4, 4.6,
    1.0, NA,  1.1,  3.5, 3.6, NA
  ), nrow = 10, byrow = TRUE)

  fit_av <- kappa_continuous(x, method = "available", weight = "quadratic")
  expect_s3_class(fit_av, "misskappa_estimate")
  expect_named(fit_av$estimates, c("Conger", "Fleiss"))
  expect_equal(unname(fit_av$estimates["Conger"]), 0.9893902893012867,
               tolerance = 1e-9)
  expect_equal(unname(fit_av$estimates["Fleiss"]), 0.9893337228643995,
               tolerance = 1e-9)

  fit_ipw <- kappa_continuous(x, method = "ipw", weight = "quadratic")
  expect_equal(unname(fit_ipw$estimates["Conger"]), 0.9889592202554336,
               tolerance = 1e-9)

  fit_gw <- kappa_continuous(x, method = "gwet", weight = "quadratic")
  expect_equal(unname(fit_gw$estimates["Conger"]), 0.9893896826901533,
               tolerance = 1e-9)
})

test_that("kappa_continuous() with identity loss on perfect agreement -> 1", {
  x <- matrix(c(1.0, 1.0, 2.5, 2.5, 3.3, 3.3), nrow = 3, byrow = TRUE)
  fit <- kappa_continuous(x, method = "available", weight = "quadratic")
  expect_equal(unname(fit$estimates["Conger"]), 1.0, tolerance = 1e-9)
})
