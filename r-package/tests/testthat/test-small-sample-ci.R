test_that("complete categorical kappa carries exact delete-one summaries", {
  x <- matrix(c(
    1, 1, 2,
    2, 2, 2,
    3, 2, 3,
    2, 3, 3,
    1, 2, 1,
    3, 3, 2,
    1, 3, 2,
    2, 1, 3
  ), nrow = 8, byrow = TRUE)

  for (weight in c("nominal", "linear", "quadratic")) {
    fit <- agreement(x, estimator = "ipw", weight = weight)
    expect_equal(dim(fit$delete_one$leave_one_out), c(nrow(x), 3L))
    expect_equal(unname(fit$delete_one$estimate), unname(coef(fit)),
                 tolerance = 1e-10)

    ci <- confint(fit, method = "jackknife", transform = "fisher")
    expect_equal(dim(ci), c(3L, 2L))
    expect_true(all(is.finite(ci)))
    expect_true(all(ci[, 1L] <= ci[, 2L]))
    expect_identical(attr(ci, "method"), "jackknife")
  }
})

test_that("counts and continuous kappa expose the same interval option", {
  counts <- agreement_counts(dat.fleiss1971, weight = "quadratic")
  expect_equal(dim(counts$delete_one$leave_one_out),
               c(nrow(dat.fleiss1971), 2L))
  expect_true(all(is.finite(confint(
    counts, method = "jackknife", transform = "log_complement"
  ))))

  ratings <- matrix(c(
    1.0, 1.2, 1.1,
    2.0, 2.1, 2.0,
    3.0, 2.8, 3.2,
    2.0, 2.2, 2.4,
    1.0, 1.5, 1.2,
    3.0, 2.9, 2.7
  ), nrow = 6, byrow = TRUE)
  continuous <- agreement(ratings, estimator = "pairwise")
  expect_equal(dim(continuous$delete_one$leave_one_out), c(6L, 2L))
  expect_true(all(is.finite(confint(
    continuous, method = "jackknife", transform = "fisher"
  ))))
})

test_that("perfect kappa uses the boundary fallback", {
  counts <- matrix(c(
    4, 0, 0,
    0, 4, 0,
    0, 0, 4,
    4, 0, 0,
    0, 4, 0,
    0, 0, 4
  ), nrow = 6, byrow = TRUE)
  fit <- agreement_counts(counts)
  ci <- confint(fit, method = "jackknife", transform = "fisher")

  expect_true(all(attr(ci, "boundary_fallback")))
  expect_equal(unname(ci[, 2L]), c(1, 1))
  expect_true(all(is.finite(ci[, 1L])))
  expect_true(all(ci[, 1L] < 1))
})

test_that("jackknife intervals remain complete-data only", {
  x <- matrix(c(
    1, 1, 2,
    2, NA, 2,
    3, 2, 3,
    2, 3, 3,
    1, 2, 1,
    3, 3, 2
  ), nrow = 6, byrow = TRUE)
  fit <- agreement(x, estimator = "ipw", weight = "nominal")
  expect_null(fit$delete_one)
  expect_error(confint(fit, method = "jackknife"),
               "unavailable for this fit")
})

test_that("Wald remains the unchanged default and accepts log-complement", {
  fit <- agreement_counts(dat.fleiss1971)
  ci <- confint(fit)
  z <- qnorm(0.975)
  expected <- cbind(coef(fit) - z * sqrt(diag(vcov(fit))),
                    coef(fit) + z * sqrt(diag(vcov(fit))))
  expect_equal(as.numeric(ci), as.numeric(expected), tolerance = 1e-12)

  transformed <- confint(fit, transform = "log_complement")
  expect_true(all(transformed[, 2L] < 1))
})
