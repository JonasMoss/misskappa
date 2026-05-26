test_that("kappa(type=) dispatches to the right implementation", {
  x_raw <- matrix(c(1, 1, NA,
                    2, 2, 2,
                    1, 2, 1),
                  nrow = 3, byrow = TRUE)

  res <- kappa(x_raw, type = "raw", method = "available")
  expect_equal(res$estimates, kappa_raw(x_raw, method = "available")$estimates)
  expect_equal(res$vcov, kappa_raw(x_raw, method = "available")$vcov)
  expect_equal(attr(res, "type", exact = TRUE), "raw")

  x_cont <- matrix(c(0.1, 0.1, NA,
                     0.0, 0.2, 0.1,
                     0.9, 1.0, 0.8),
                   nrow = 3, byrow = TRUE)

  res <- kappa(x_cont, type = "continuous", method = "available", weight = "identity")
  expect_equal(
    res$estimates,
    kappa_continuous(x_cont, method = "available", weight = "identity")$estimates
  )
  expect_equal(
    res$vcov,
    kappa_continuous(x_cont, method = "available", weight = "identity")$vcov
  )
  expect_equal(attr(res, "type", exact = TRUE), "continuous")

  x_counts <- matrix(c(2, 1, 0,
                       1, 1, 1,
                       0, 3, 0),
                     nrow = 3, byrow = TRUE)

  res <- kappa(x_counts, type = "counts", method = "available", r = 3)
  expect_equal(res$estimates, kappa_counts(x_counts, method = "available", r = 3)$estimates)
  expect_equal(res$vcov, kappa_counts(x_counts, method = "available", r = 3)$vcov)
  expect_equal(attr(res, "type", exact = TRUE), "counts")
})

test_that("kappa(type='auto') uses expected heuristics", {
  x_cont <- matrix(c(0.1, 0.1, NA,
                     0.0, 0.2, 0.1,
                     0.9, 1.0, 0.8),
                   nrow = 3, byrow = TRUE)
  res <- kappa(x_cont, type = "auto", method = "available", weight = "identity")
  expect_equal(res$estimates, kappa_continuous(x_cont, method = "available", weight = "identity")$estimates)
  expect_equal(attr(res, "type", exact = TRUE), "continuous")

  x_raw <- matrix(c(1, 1, NA,
                    2, 2, 2,
                    1, 2, 1),
                  nrow = 3, byrow = TRUE)
  res <- kappa(x_raw, type = "auto", method = "available", weight = "quadratic")
  expect_equal(res$estimates, kappa_raw(x_raw, method = "available", weight = "quadratic")$estimates)
  expect_equal(attr(res, "type", exact = TRUE), "raw")

  x_counts <- matrix(c(2, 1, 0,
                       1, 1, 1,
                       0, 3, 0),
                     nrow = 3, byrow = TRUE)
  res <- kappa(x_counts, type = "auto", method = "available", r = 3)
  expect_equal(res$estimates, kappa_counts(x_counts, method = "available", r = 3)$estimates)
  expect_equal(attr(res, "type", exact = TRUE), "counts")
})
