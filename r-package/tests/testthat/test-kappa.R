test_that("kappa() routes the four methods and returns the S3 shape", {
  x <- matrix(
    c(0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1),
    nrow = 10, byrow = TRUE
  )
  storage.mode(x) <- "integer"

  for (m in c("available", "ipw", "fiml", "gwet")) {
    fit <- estimate_kappa_raw(x, method = m, weight = "identity")
    expect_s3_class(fit, "misskappa_estimate")
    expect_named(fit$estimates, c("Conger", "Fleiss", "Brennan-Prediger"))
    expect_equal(dim(fit$vcov), c(3L, 3L))
    expect_true(all(is.finite(fit$estimates)))
  }
})

test_that("alpha() routes numeric available and normal FIML paths", {
  set.seed(11)
  n <- 120L
  p <- 4L
  L <- chol(0.25 + 0.75 * diag(p))
  x <- matrix(rnorm(n * p), n, p) %*% L

  fit_av <- alpha(x, estimator = "pairwise")
  fit_ml <- alpha(x, estimator = "nt_fiml")

  expect_s3_class(fit_av, "misskappa_estimate")
  expect_s3_class(fit_ml, "misskappa_estimate")
  expect_named(fit_av$estimates, "alpha")
  expect_equal(dim(fit_av$vcov), c(1L, 1L))
  expect_equal(dim(fit_av$psi), c(n, 1L))
  expect_equal(unname(fit_ml$estimates), unname(fit_av$estimates),
               tolerance = 1e-8)
})

test_that("alpha_cat_fiml() exposes saturated categorical FIML", {
  x <- matrix(c(
    1, 1, 2,
    2, 2, 2,
    3, 2, 3,
    2, 3, 3,
    1, 2, 1,
    3, 3, 2
  ), nrow = 6, byrow = TRUE)
  fit_av <- alpha(x, estimator = "pairwise")
  fit_ml <- alpha_cat_fiml(x)

  expect_s3_class(fit_av, "misskappa_estimate")
  expect_s3_class(fit_ml, "misskappa_estimate")
  expect_named(fit_av$estimates, "alpha")
  expect_equal(dim(fit_av$vcov), c(1L, 1L))
  expect_equal(dim(fit_av$psi), c(6L, 1L))
  expect_equal(unname(fit_ml$estimates), unname(fit_av$estimates),
               tolerance = 1e-8)
})

test_that("available-case on a Cohen-style 2-rater example matches the textbook value", {
  x <- matrix(
    c(0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1),
    nrow = 10, byrow = TRUE
  )
  fit <- estimate_kappa_raw(x, method = "available", weight = "identity")
  # Cohen's kappa for this contingency table is 0.4.
  expect_equal(unname(fit$estimates["Conger"]), 0.4, tolerance = 1e-9)
  expect_equal(unname(fit$estimates["Brennan-Prediger"]), 0.4, tolerance = 1e-9)
})

test_that("S3 methods on misskappa_estimate behave", {
  x <- matrix(c(0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1), nrow = 6, byrow = TRUE)
  fit <- estimate_kappa_raw(x, method = "available", weight = "identity")

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

  expect_output(print(fit), "estimator=available")
})

test_that("IPW rejects a rater with zero observations cleanly", {
  x <- matrix(c(0, 0, NA, 1, 1, NA, 0, 1, NA, 1, 0, NA, 0, 0, NA),
              nrow = 5, byrow = TRUE)
  expect_error(kappa(x, estimator = "ipw"),
               regexp = "every rater must be observed")
})

test_that("FIML matches available on complete data, identity weight", {
  x <- matrix(
    c(0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1),
    nrow = 10, byrow = TRUE
  )
  fit_av <- estimate_kappa_raw(x, method = "available", weight = "identity")
  fit_ml <- estimate_kappa_raw(x, method = "fiml", weight = "identity")
  expect_equal(unname(fit_ml$estimates), unname(fit_av$estimates),
               tolerance = 1e-6)
})

test_that("multi-weight FIML helper matches repeated FIML fits", {
  complete_support <- as.matrix(expand.grid(
    r1 = 0:1, r2 = 0:1, r3 = 0:1, r4 = 0:1
  ))
  missing_rows <- matrix(c(
    0, NA, 0, 1,
    1, 0, NA, 0,
    1, 1, 1, NA
  ), ncol = 4, byrow = TRUE)
  x <- rbind(complete_support, missing_rows)

  multi <- estimate_kappa_fiml_multi(
    x,
    weights = c("identity", "linear", "quadratic"),
    values = c(0, 1)
  )
  expect_named(multi, c("identity", "linear", "quadratic"))

  for (w in names(multi)) {
    single <- estimate_kappa_raw(x, method = "fiml", weight = w,
                                 values = c(0, 1))
    expect_equal(unname(multi[[w]]$estimates), unname(single$estimates),
                 tolerance = 1e-10)
    expect_equal(unname(multi[[w]]$vcov), unname(single$vcov),
                 tolerance = 1e-10)
    expect_equal(unname(multi[[w]]$psi), unname(single$psi),
                 tolerance = 1e-10)
    expect_equal(multi[[w]]$method, "fiml")
    expect_equal(multi[[w]]$weight, w)
  }
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

test_that("kappa_quadratic() matches legacy on a small raw fixture", {
  x <- matrix(c(
    1, 2, 1,  2, 2, 3,  3, 3, 3,  4, 4, 5,
    5, 5, 5,  1, 1, 2,  3, 3, 4,  4, 5, 4,
    2, 3, 2,  5, 4, 5,  1, 2, NA, 3, NA, 4
  ), nrow = 12, byrow = TRUE)
  fit <- kappa_quadratic(x, values = c(1, 2, 3, 4, 5))
  expect_named(fit$estimates, c("Conger", "Fleiss", "Brennan-Prediger"))
  expect_equal(unname(fit$estimates["Conger"]), 0.8374174174174175,
               tolerance = 1e-9)
  expect_equal(unname(fit$estimates["Fleiss"]), 0.8344881466279811,
               tolerance = 1e-9)
})

test_that("kappa_quadratic() exposes empirical covariance and influence rows", {
  x <- matrix(c(
    1, 2, 1,
    2, 2, 3,
    3, 3, 3,
    4, 4, 5,
    5, 5, 5,
    NA, NA, NA
  ), nrow = 6, byrow = TRUE)
  fit <- kappa_quadratic(x, values = c(1, 2, 3, 4, 5))
  expect_true(all(is.finite(fit$vcov)))
  expect_equal(dim(fit$psi), c(6L, 3L))
  expect_equal(unname(crossprod(fit$psi) / 36),
               unname(vcov(fit)), tolerance = 1e-10)
})

test_that("kappa_counts(method='fiml') matches fleiss_cuzick on Fleiss 1971", {
  fit_av <- kappa_counts(dat.fleiss1971, estimator = "fleiss_cuzick", weight = "nominal")
  fit_ml <- kappa_counts(dat.fleiss1971, estimator = "cat_fiml", weight = "nominal")
  expect_equal(unname(fit_av$estimates), unname(fit_ml$estimates),
               tolerance = 1e-9)
})

test_that("kappa_counts(method='fiml') differs from fleiss_cuzick with partial counts", {
  # r_total = 4 but only one row sums to 4; others sum to 3 (partial counts).
  y <- matrix(c(
    3, 1, 0,  0, 2, 1,  2, 0, 1,  1, 2, 0,
    0, 0, 3,  2, 1, 0,  3, 0, 0,  0, 1, 2
  ), nrow = 8, byrow = TRUE)
  fit_av <- kappa_counts(y, estimator = "fleiss_cuzick", weight = "nominal")
  fit_ml <- kappa_counts(y, estimator = "cat_fiml", weight = "nominal", r_total = 4)
  expect_false(isTRUE(all.equal(unname(fit_av$estimates),
                                unname(fit_ml$estimates),
                                tolerance = 1e-9)))
  # Frozen against legacy.
  expect_equal(unname(fit_ml$estimates["Fleiss"]), 0.2725053215314407,
               tolerance = 1e-7)
})

test_that("kappa_counts() reproduces Fleiss 1971", {
  fit_id <- kappa_counts(dat.fleiss1971, weight = "nominal")
  expect_s3_class(fit_id, "misskappa_estimate")
  expect_named(fit_id$estimates, c("Fleiss", "Brennan-Prediger"))
  expect_equal(unname(fit_id$estimates["Fleiss"]), 0.4302445200601408,
               tolerance = 1e-9)
  expect_equal(unname(fit_id$estimates["Brennan-Prediger"]), 4 / 9,
               tolerance = 1e-9)

  fit_q <- kappa_counts(dat.fleiss1971, weight = "quadratic",
                        values = c(1, 2, 3, 4, 5))
  expect_equal(unname(fit_q$estimates["Fleiss"]), 0.2840722495894910,
               tolerance = 1e-9)
})

test_that("kappa_counts() accepts historical pairwise alias", {
  fit_fc <- kappa_counts(dat.fleiss1971, estimator = "fleiss_cuzick")
  fit_old <- kappa_counts(dat.fleiss1971, estimator = "pairwise")
  expect_identical(fit_old$method, "fleiss_cuzick")
  expect_equal(unname(fit_old$estimates), unname(fit_fc$estimates),
               tolerance = 1e-12)
})

test_that("kappa_continuous() with identity loss on perfect agreement -> 1", {
  x <- matrix(c(1.0, 1.0, 2.5, 2.5, 3.3, 3.3), nrow = 3, byrow = TRUE)
  fit <- kappa_continuous(x, method = "available", weight = "quadratic")
  expect_equal(unname(fit$estimates["Conger"]), 1.0, tolerance = 1e-9)
})

test_that("kappa_vector() internal wrapper handles component-separable losses", {
  x <- array(c(
    0, 1,  0, 1,  1, 1,
    1, 0,  1, 0,  1, 1,
    0, 0,  1, 0,  0, 0,
    1, 1,  1, 0,  0, 1,
    0, 1,  0, 0,  0, 1
  ), dim = c(5L, 3L, 2L))

  fit_pair <- kappa_vector(x, method = "pairwise", loss = "hamming",
                           feature_weights = c(1, 2))
  fit_ipw <- kappa_vector(x, method = "ipw", loss = "hamming",
                          feature_weights = c(1, 2))
  expect_s3_class(fit_pair, "misskappa_estimate")
  expect_named(fit_pair$estimates, c("Conger", "Fleiss"))
  expect_equal(unname(fit_pair$estimates), unname(fit_ipw$estimates),
               tolerance = 1e-9)
  expect_false("kappa_vector" %in% getNamespaceExports("misskappa"))

  x[1, 1, 2] <- NA_real_
  x[2, 3, 1] <- NA_real_
  fit_rms <- kappa_vector(x, method = "ipw", loss = "rms",
                          feature_weights = c(1, 2))
  expect_true(all(is.finite(fit_rms$estimates)))
  expect_equal(dim(fit_rms$psi), c(5L, 2L))
  expect_equal(unname(crossprod(fit_rms$psi) / 25), unname(vcov(fit_rms)),
               tolerance = 1e-10)
})

test_that("kappa_gwise() estimates Frechet and Hubert complete-data coefficients", {
  x <- matrix(c(
    1, 1, 2, 1, 1,
    1, 2, 3, 2, 2,
    2, 1, 1, 1, 1,
    2, 3, 4, 4, 5
  ), nrow = 4, byrow = TRUE)

  fit_nom <- kappa_gwise(x, distance = "nominal")
  expect_s3_class(fit_nom, "misskappa_estimate")
  expect_named(fit_nom$estimates, c("Conger", "Fleiss"))
  expect_equal(unname(fit_nom$estimates["Conger"]), 0.23353293413173815,
               tolerance = 1e-9)
  expect_equal(unname(fit_nom$estimates["Fleiss"]), 0.21735788407113500,
               tolerance = 1e-9)
  expect_equal(dim(fit_nom$psi), c(4L, 2L))

  fit_abs <- kappa_gwise(x, distance = "absolute")
  expect_equal(unname(fit_abs$estimates["Conger"]), 0.45877378435517835,
               tolerance = 1e-9)
  expect_equal(unname(fit_abs$estimates["Fleiss"]), 0.44821878125323300,
               tolerance = 1e-9)

  fit_hubert <- kappa_gwise(matrix(c(1, 1, 1, 2, 2, 2), nrow = 2, byrow = TRUE),
                            distance = "hubert")
  expect_equal(unname(fit_hubert$estimates), c(1, 1), tolerance = 1e-9)
})
