# This test file verifies that for complete data, different estimation methods
# produce identical results where they are theoretically expected to.
# Specifically:
# 1. For quadratic weights, 'np', 'ipw', and 'quadratic' methods should match.
# 2. For all other weights, 'ml', 'np', and 'ipw' should match.

# Load the necessary library and data
library(testthat)
library(misskappa)

# Use a small, complete dataset for testing
# dat.zapf2016 has no missing values. We take a subset for speed.
complete_data <- dat.zapf2016[c(1:5), ]

# Define the set of weights to be tested
# We exclude "quadratic" from this main loop as it's tested separately.
weights_to_test <- c(
  "identity", "unweighted", "linear", "ordinal",
  "radical", "ratio", "circular", "bipolar"
)


test_that("Methods are equivalent for quadratic weights on complete data", {
  # For complete data with quadratic weights, the nonparametric, IPW, and
  # parametric quadratic methods should all yield identical results.

  # Calculate results from the three methods
  res_quad <- kappa_raw(complete_data, method = "quadratic")
  res_np <- kappa_raw(complete_data, method = "available", weight = "quadratic")
  res_ipw <- kappa_raw(complete_data, method = "ipw", weight = "quadratic")

  # Compare estimates
  expect_equal(
    res_quad$estimates,
    res_np$estimates,
    tolerance = 1e-7,
    info = "Estimates from 'quadratic' and 'np' methods should match."
  )
  expect_equal(
    res_np$estimates,
    res_ipw$estimates,
    tolerance = 1e-7,
    info = "Estimates from 'np' and 'ipw' methods should match."
  )

  # Compare variance-covariance matrices
  expect_equal(
    res_quad$vcov,
    res_np$vcov,
    tolerance = 1e-7,
    info = "V-Cov matrices from 'quadratic' and 'np' methods should match."
  )
  expect_equal(
    res_np$vcov,
    res_ipw$vcov,
    tolerance = 1e-7,
    info = "V-Cov matrices from 'np' and 'ipw' methods should match."
  )
})


# Loop through the other weight types to test ml, np, and ipw
for (w in weights_to_test) {
  test_that(paste("Methods (ml, np, ipw) are equivalent for '", w, "' weights on complete data"), {
    # On complete data, the Maximum Likelihood, Nonparametric, and IPW
    # estimators should produce identical point estimates and V-Cov matrices
    # because the missing data mechanisms (EM algorithm, inverse weights)
    # have no effect.

    # Calculate results from the three methods
    # EM options are set to high precision for robust comparison.
    res_ml <- kappa_raw(complete_data, method = "ml", weight = w, em_options = list(tol = 1e-12))
    res_np <- kappa_raw(complete_data, method = "available", weight = w)
    res_ipw <- kappa_raw(complete_data, method = "ipw", weight = w)

    # Compare estimates
    expect_equal(
      res_ml$estimates,
      res_np$estimates,
      tolerance = 1e-7,
      info = paste("Estimates from 'ml' and 'np' should match for weight:", w)
    )
    expect_equal(
      res_np$estimates,
      res_ipw$estimates,
      tolerance = 1e-7,
      info = paste("Estimates from 'np' and 'ipw' should match for weight:", w)
    )

    # Compare variance-covariance matrices
    expect_equal(
      res_ml$vcov,
      res_np$vcov,
      tolerance = 1e-7,
      info = paste("V-Cov matrices from 'ml' and 'np' should match for weight:", w)
    )
    expect_equal(
      res_np$vcov,
      res_ipw$vcov,
      tolerance = 1e-7,
      info = paste("V-Cov matrices from 'np' and 'ipw' should match for weight:", w)
    )
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ADDENDUM: INCOMPLETE DATA
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

test_that("Methods for incomplete raw data produce distinct results", {
  # On incomplete data, the different methods ('ml', 'np', 'ipw', 'quadratic')
  # represent distinct statistical approaches and are expected to produce
  # different results.

  incomplete_data <- dat.klein2018

  # Check they run without error
  res_ml <- kappa_raw(incomplete_data, method = "ml", weight = "quadratic")
  res_np <- kappa_raw(incomplete_data, method = "available", weight = "quadratic")
  res_ipw <- kappa_raw(incomplete_data, method = "ipw", weight = "quadratic")
  res_quad <- kappa_raw(incomplete_data, method = "quadratic")

  # Use expect_false(isTRUE(all.equal(...))) for robust inequality check
  # We check the estimates vector, which is sufficient.
  expect_false(
    isTRUE(all.equal(res_ml$estimates, res_np$estimates)),
    info = "ml and np should differ on incomplete raw data."
  )
  expect_false(
    isTRUE(all.equal(res_np$estimates, res_ipw$estimates)),
    info = "np and ipw should differ on incomplete raw data."
  )
  expect_false(
    isTRUE(all.equal(res_ipw$estimates, res_quad$estimates)),
    info = "ipw and quadratic should differ on incomplete raw data."
  )
})
