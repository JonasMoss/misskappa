test_that("available raw complete-data matches irrcacsmoke", {
  skip_if_not_installed("irrcacsmoke")

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
  oracle <- rbind(
    Conger = irrcacsmoke::conger.kappa.raw(x),
    Fleiss = irrcacsmoke::fleiss.kappa.raw(x),
    `Brennan-Prediger` = irrcacsmoke::bp.coeff.raw(x)
  )

  expect_true(is.numeric(oracle))
  expect_true(all(is.finite(oracle[, "estimate"])))
  expect_true(all(is.finite(oracle[, "variance"])))
  expect_true(all(oracle[, "variance"] >= 0))
  expect_equal(unname(fit$estimates[rownames(oracle)]),
               unname(oracle[, "estimate"]),
               tolerance = 1e-9)
  expect_equal(unname(diag(fit$vcov)[rownames(oracle)]),
               unname(oracle[, "variance"]),
               tolerance = 1e-9)
})

test_that("counts complete-data matches irrcacsmoke", {
  skip_if_not_installed("irrcacsmoke")

  fit <- kappa_counts(dat.fleiss1971, weight = "nominal")
  oracle <- rbind(
    Fleiss = irrcacsmoke::fleiss.kappa.dist(dat.fleiss1971),
    `Brennan-Prediger` = irrcacsmoke::bp.coeff.dist(dat.fleiss1971)
  )

  expect_true(is.numeric(oracle))
  expect_true(all(is.finite(oracle[, "estimate"])))
  expect_true(all(is.finite(oracle[, "variance"])))
  expect_true(all(oracle[, "variance"] >= 0))
  expect_equal(unname(fit$estimates[rownames(oracle)]),
               unname(oracle[, "estimate"]),
               tolerance = 1e-9)
  expect_equal(unname(diag(fit$vcov)[rownames(oracle)]),
               unname(oracle[, "variance"]),
               tolerance = 1e-9)
})
