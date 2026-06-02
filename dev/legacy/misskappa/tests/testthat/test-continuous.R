library(testthat)
library(misskappa)

complete_data <- dat.zapf2016
res_quad <- kappa_continuous(complete_data, method = "quadratic")
res_np <- kappa_continuous(complete_data, method = "available", weight = "quadratic")
res_ipw <- kappa_continuous(complete_data, method = "ipw", weight = "quadratic")

test_that("Continuous methods are identical on complete data with quadratic weights", {
  # Compare np vs quadratic
  expect_equal(
    res_np$estimates,
    res_quad$estimates,
    tolerance = 1e-9,
    info = "Estimates from 'np' and 'quadratic' should match on complete data."
  )
  expect_equal(
    res_np$vcov,
    res_quad$vcov,
    tolerance = 1e-9,
    info = "V-Cov from 'np' and 'quadratic' should match on complete data."
  )

  # Compare ipw vs quadratic (should also be identical as weights become 1)
  expect_equal(
    res_ipw$estimates,
    res_quad$estimates,
    tolerance = 1e-9,
    info = "Estimates from 'ipw' and 'quadratic' should match on complete data."
  )
  expect_equal(
    res_ipw$vcov,
    res_quad$vcov,
    tolerance = 1e-9,
    info = "V-Cov from 'ipw' and 'quadratic' should match on complete data."
  )
})

incomplete_data <- dat.klein2018

test_that("Continuous methods run without error on incomplete data", {
  # Just check that they run and return the correct structure
  expect_no_error(kappa_continuous(incomplete_data, method = "quadratic"))
  expect_no_error(kappa_continuous(incomplete_data, method = "available"))
  expect_no_error(kappa_continuous(incomplete_data, method = "ipw"))

  res <- kappa_continuous(incomplete_data, method = "ipw")
  expect_named(res, c("estimates", "vcov"))
  expect_length(res$estimates, 2)
  expect_equal(dim(res$vcov), c(2, 2))
})

test_that("'np' and 'ipw' methods differ on incomplete data", {
  res_np <- kappa_continuous(incomplete_data, method = "available", weight = "quadratic")
  res_ipw <- kappa_continuous(incomplete_data, method = "ipw", weight = "quadratic")

  # The estimates should NOT be equal
  expect_false(
    isTRUE(all.equal(res_np$estimates, res_ipw$estimates)),
    info = "'np' and 'ipw' estimates should differ on incomplete data."
  )
  expect_false(
    isTRUE(all.equal(res_np$vcov, res_ipw$vcov)),
    info = "'np' and 'ipw' V-Cov matrices should differ on incomplete data."
  )
})

test_that("Different weights produce different results", {
  # Using complete data for simplicity
  res_quad <- kappa_continuous(complete_data, method = "available", weight = "quadratic")
  res_linear <- kappa_continuous(complete_data, method = "available", weight = "linear")

  # Quadratic vs Linear
  expect_false(
    isTRUE(all.equal(res_quad$estimates, res_linear$estimates)),
    info = "Quadratic and linear weights should produce different estimates."
  )
})

test_that("'np' and 'quadratic' methods differ on incomplete continuous data", {
  incomplete_data <- dat.klein2018

  res_np <- kappa_continuous(incomplete_data, method = "available", weight = "quadratic")
  res_quad <- kappa_continuous(incomplete_data, method = "quadratic")

  expect_false(
    isTRUE(all.equal(res_np$estimates, res_quad$estimates)),
    info = "'np' and 'quadratic' estimates should differ on incomplete continuous data."
  )
})

test_that("All supported weights run for continuous methods without error", {
  # Define the weights applicable to the continuous 'np'/'ipw' methods
  continuous_weights <- c("linear", "quadratic", "radical", "ratio")

  # Using incomplete data to ensure IPW is also tested meaningfully
  data <- dat.klein2018

  for (w in continuous_weights) {
    # Test the 'np' method
    expect_no_error(
      kappa_continuous(data, method = "available", weight = w)
    )

    # Test the 'ipw' method
    expect_no_error(
      kappa_continuous(data, method = "ipw", weight = w)
    )
  }

  # Also confirm that unsupported weights correctly throw an error
  unsupported_weights <- c("ordinal", "circular", "bipolar")

  for (w in unsupported_weights) {
    expect_error(
      kappa_continuous(data, method = "available", weight = w)
    )
  }
})
