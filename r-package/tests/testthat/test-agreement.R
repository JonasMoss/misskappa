test_that("agreement() routes the four methods and returns the S3 shape", {
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
  expect_named(df, c("estimator", "coefficient", "estimate", "se", "lower",
                     "upper"))
  expect_equal(df$lower, unname(ci[, 1]))

  expect_output(print(fit), "estimator = available")
  expect_output(print(fit), "n = ")
})

test_that("IPW rejects a rater with zero observations cleanly", {
  x <- matrix(c(0, 0, NA, 1, 1, NA, 0, 1, NA, 1, 0, NA, 0, 0, NA),
              nrow = 5, byrow = TRUE)
  expect_error(agreement(x, estimator = "ipw"),
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

test_that("agreement_counts(method='fiml') matches fleiss_cuzick on Fleiss 1971", {
  fit_av <- agreement_counts(dat.fleiss1971, estimator = "fleiss_cuzick", weight = "nominal")
  fit_ml <- agreement_counts(dat.fleiss1971, estimator = "cat_fiml", weight = "nominal")
  expect_equal(unname(fit_av$estimates), unname(fit_ml$estimates),
               tolerance = 1e-9)
})

test_that("agreement_counts(method='fiml') differs from fleiss_cuzick with partial counts", {
  # r_total = 4 but only one row sums to 4; others sum to 3 (partial counts).
  y <- matrix(c(
    3, 1, 0,  0, 2, 1,  2, 0, 1,  1, 2, 0,
    0, 0, 3,  2, 1, 0,  3, 0, 0,  0, 1, 2
  ), nrow = 8, byrow = TRUE)
  fit_av <- agreement_counts(y, estimator = "fleiss_cuzick", weight = "nominal")
  fit_ml <- agreement_counts(y, estimator = "cat_fiml", weight = "nominal", r_total = 4)
  expect_false(isTRUE(all.equal(unname(fit_av$estimates),
                                unname(fit_ml$estimates),
                                tolerance = 1e-9)))
  # Frozen against legacy.
  expect_equal(unname(fit_ml$estimates["Fleiss"]), 0.2725053215314407,
               tolerance = 1e-7)
})

test_that("ratings_to_counts() aggregates and discards rater identity", {
  X <- rbind(
    c(0, 1, 1, NA),
    c(2, 2, NA, 2),
    c(0, 0, 0, 0)
  )
  cnt <- ratings_to_counts(X)
  expect_equal(colnames(cnt), c("0", "1", "2"))
  expect_equal(unname(cnt[1, ]), c(1L, 2L, 0L))   # one 0, two 1s, one missing
  expect_equal(unname(cnt[2, ]), c(0L, 0L, 3L))   # three 2s, one missing
  expect_equal(rowSums(cnt), c(3, 3, 4))          # row total = # observed
  expect_true(is.integer(cnt))
  # Permuting raters within each subject leaves the counts unchanged.
  Xperm <- t(apply(X, 1L, sample))
  expect_equal(ratings_to_counts(Xperm, categories = c(0, 1, 2)), cnt)
})

test_that("ratings_to_counts() honours an explicit category set and errors off-set", {
  X <- rbind(c(1, 1, 2), c(2, 2, 1))
  cnt <- ratings_to_counts(X, categories = 1:4)
  expect_equal(colnames(cnt), as.character(1:4))
  expect_equal(unname(colSums(cnt)), c(3L, 3L, 0L, 0L))  # cats 3,4 unobserved
  expect_error(ratings_to_counts(X, categories = c(1, 3)), "absent from")
  expect_error(ratings_to_counts(matrix(NA_real_, 2, 2)), "no observed")
})

test_that("ratings_to_counts() -> cat_fiml equals exchangeable-restricted raw fit", {
  # On counts-summarised data, agreement_counts(cat_fiml) is the exchangeable
  # count-FIML; check the bridge reproduces a hand-built counts matrix and
  # that fleiss_cuzick is reachable the same way.
  set.seed(11)
  X <- matrix(sample(c(0:2, NA), 60, replace = TRUE, prob = c(.3, .3, .3, .1)),
              nrow = 15)
  cnt <- ratings_to_counts(X)
  ref <- t(apply(X, 1L, function(r) tabulate(match(r, 0:2), nbins = 3L)))
  expect_equal(unname(cnt), ref)
  expect_s3_class(agreement_counts(cnt, estimator = "cat_fiml", r_total = 4),
                  "misskappa_estimate")
  expect_s3_class(agreement_counts(cnt, estimator = "fleiss_cuzick"),
                  "misskappa_estimate")
})

test_that("agreement_counts(nt_fiml) matches raw exchangeable normal FIML", {
  set.seed(12)
  X <- matrix(sample(c(1:5, NA), 240, replace = TRUE,
                     prob = c(rep(0.18, 5), 0.10)), nrow = 60)
  cnt <- ratings_to_counts(X, categories = 1:5)

  # Oracle: the exchangeable (equal-mean, compound-symmetry) raw normal FIML.
  # mean = "equal" is not a public agreement() option (mean-free surface), so the
  # internal backend is exercised directly.
  raw <- kappa_quadratic_fiml(
    X, values = 1:5,
    em_options = list(covariance = "variance_components", mean = "equal"))
  by_counts <- agreement_counts(cnt, estimator = "nt_fiml_cs", values = 1:5)

  # The direct-on-counts fit (BFGS over the sufficient statistics) and the raw
  # exchangeable normal FIML (C++ EM) target the same MLE; they agree to
  # optimiser tolerance, not bit-for-bit.
  expect_named(coef(by_counts), c("Fleiss", "Brennan-Prediger"))
  expect_lt(max(abs(coef(by_counts) -
                    coef(raw)[c("Fleiss", "Brennan-Prediger")])), 1e-5)
  expect_lt(max(abs(vcov(by_counts) -
                    vcov(raw)[c("Fleiss", "Brennan-Prediger"),
                              c("Fleiss", "Brennan-Prediger")])), 1e-5)
  expect_equal(unname(crossprod(by_counts$psi) / nrow(cnt)^2),
               unname(vcov(by_counts)), tolerance = 1e-10)
  expect_equal(by_counts$method, "nt_fiml_cs")
  expect_null(by_counts$mean)
  expect_null(by_counts$covariance)
})

test_that("agreement_counts(nt_fiml_cs) is quadratic only and defaults to it", {
  expect_error(
    agreement_counts(dat.fleiss1971, estimator = "nt_fiml_cs", weight = "nominal"),
    "requires weight = \"quadratic\""
  )
  fit <- agreement_counts(dat.fleiss1971, estimator = "nt_fiml_cs")
  expect_equal(fit$weight, "quadratic")
})

test_that("agreement_counts() reproduces Fleiss 1971", {
  fit_id <- agreement_counts(dat.fleiss1971, weight = "nominal")
  expect_s3_class(fit_id, "misskappa_estimate")
  expect_named(fit_id$estimates, c("Fleiss", "Brennan-Prediger"))
  expect_equal(unname(fit_id$estimates["Fleiss"]), 0.4302445200601408,
               tolerance = 1e-9)
  expect_equal(unname(fit_id$estimates["Brennan-Prediger"]), 4 / 9,
               tolerance = 1e-9)

  fit_q <- agreement_counts(dat.fleiss1971, weight = "quadratic",
                        values = c(1, 2, 3, 4, 5))
  expect_equal(unname(fit_q$estimates["Fleiss"]), 0.2840722495894910,
               tolerance = 1e-9)
})

test_that("counts CS nt_fiml coincides with fleiss_cuzick on complete data", {
  # On complete balanced counts (every row sums to r_total) the
  # compound-symmetry FIML and the Fleiss-Cuzick moment estimator collapse to
  # the same point and robust SE.
  vals <- 1:5
  fc <- agreement_counts(dat.fleiss1971, estimator = "fleiss_cuzick",
                     weight = "quadratic", values = vals)
  f  <- agreement_counts(dat.fleiss1971, estimator = "nt_fiml_cs", values = vals)
  expect_named(coef(f), c("Fleiss", "Brennan-Prediger"))
  expect_equal(unname(coef(f)), unname(coef(fc)), tolerance = 1e-5)
  expect_equal(unname(diag(vcov(f))), unname(diag(vcov(fc))), tolerance = 1e-5)
  # robust influence-function identity
  expect_equal(unname(crossprod(f$psi) / nrow(dat.fleiss1971)^2),
               unname(vcov(f)), tolerance = 1e-10)
})

test_that("counts pairwise and listwise are not valid estimators", {
  expect_error(agreement_counts(dat.fleiss1971, estimator = "pairwise"),
               "should be one of")
  expect_error(agreement_counts(dat.fleiss1971, estimator = "listwise"),
               "should be one of")
})

test_that("nt_fiml / nt_fiml_cs select the covariance structure by name", {
  set.seed(31)
  X <- matrix(rnorm(120 * 4), 120, 4)
  X[cbind(sample(120, 30), sample(4, 30, TRUE))] <- NA

  sat <- agreement(X, estimator = "nt_fiml", weight = "quadratic")
  cs  <- agreement(X, estimator = "nt_fiml_cs", weight = "quadratic")
  expect_named(coef(cs), c("Conger", "Fleiss"))   # no declared scale, no BP
  expect_equal(sat$method, "nt_fiml")
  expect_equal(cs$method, "nt_fiml_cs")
  expect_null(cs$covariance)            # name carries the structure
  expect_null(cs$mean)
  # nt_fiml_cs matches the internal compound-symmetry backend
  ref <- kappa_quadratic_fiml(
    X, em_options = list(covariance = "variance_components"))
  expect_equal(unname(coef(cs)), unname(coef(ref)[1:2]), tolerance = 1e-6)
  # saturated and CS are genuinely different fits under missingness
  expect_gt(max(abs(coef(sat) - coef(cs))), 1e-4)
})

test_that("pairwise_cs is the moment companion to nt_fiml_cs", {
  set.seed(131)
  n <- 4000L; R <- 4L; mu <- c(0, .3, .6, 1)
  b <- rnorm(n, 0, 1)
  X <- matrix(mu, n, R, byrow = TRUE) + b + matrix(rnorm(n * R, 0, sqrt(.5)), n, R)

  fit <- .kappa_pairwise_cs(X)          # internal since 2.0.0
  expect_equal(fit$method, "pairwise_cs")
  expect_equal(fit$weight, "quadratic")
  expect_equal(names(coef(fit)), c("Conger", "Fleiss", "Brennan-Prediger"))
  expect_equal(dim(vcov(fit)), c(3L, 3L))
  # influence contract
  expect_equal(unname(crossprod(fit$psi) / n^2), unname(vcov(fit)),
               tolerance = 1e-10)
  # shares the probability limit with nt_fiml_cs (both consistent for the same
  # CS estimand); agree on large complete data
  cs <- kappa_quadratic_fiml(
    X, em_options = list(covariance = "variance_components"))
  expect_equal(unname(coef(fit)), unname(coef(cs)), tolerance = 0.01)
})

test_that("internal pairwise_cs allows connected incomplete co-observation", {
  set.seed(132)
  n <- 300L; R <- 4L
  b <- rnorm(n, 0, 1)
  X <- matrix(0, n, R) + b + matrix(rnorm(n * R, 0, sqrt(.5)), n, R)
  # connected but NOT complete: group A sees raters 1-3, group B sees 2-4, so the
  # pair (1, 4) is never co-observed yet the co-observation graph is connected.
  half <- n / 2
  X[seq_len(half), 4L] <- NA
  X[(half + 1L):n, 1L] <- NA
  fit <- expect_silent(.kappa_pairwise_cs(X))
  expect_true(all(is.finite(coef(fit))))
  expect_true(all(is.finite(diag(vcov(fit)))))
  expect_error(agreement(X, estimator = "pairwise_cs"), "should be one of")
})

test_that("em_options are validated per estimator", {
  set.seed(33)
  X <- matrix(rnorm(120 * 4), 120, 4)
  X[cbind(sample(120, 25), sample(4, 25, TRUE))] <- NA
  expect_error(agreement(X, estimator = "nt_fiml",
                     em_options = list(covariance = "variance_components")),
               "unknown em_options")
  expect_error(agreement(X, estimator = "pairwise", em_options = list(tol = 1e-6)),
               "does not apply")
  expect_error(agreement(X, estimator = "nt_fiml", em_options = list(flatten = 1)),
               "unknown em_options")
  expect_error(agreement_counts(dat.fleiss1971, estimator = "cat_fiml",
                            em_options = list(flatten = 1)),
               "unknown em_options")
  fit <- agreement(X, estimator = "nt_fiml", em_options = list(tol = 1e-9))
  expect_s3_class(fit, "misskappa_estimate")
  expect_error(agreement(X, estimator = "nt_fiml", W = diag(2)), "not used by")
  expect_error(agreement(X, estimator = "ipw", foo = 1), "not used by")
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
