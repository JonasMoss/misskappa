library(testthat)
library(misskappa)

# Use a complete dataset for straightforward comparison
raw_data <- dat.zapf2016
counts_data <- to_counts_matrix(raw_data)
r_total <- ncol(raw_data) # Total number of raters

# Define the set of weights to be tested
weights_to_test <- c(
  "identity", "unweighted", "linear", "quadratic", "ordinal",
  "radical", "ratio", "circular", "bipolar"
)


test_that("Counts 'np' method matches raw 'np' method", {
  # The `kappa_counts(..., method = 'np')` is a wrapper around the raw data
  # version. This test ensures the wrapper works correctly and that the
  # underlying theory (that np is invariant to rater identity) holds.

  # Calculate results directly from the raw data
  res_raw_np <- kappa_raw(raw_data, method = "available", weight = "quadratic")

  # Calculate results from the counts data using the np wrapper
  res_counts_np <- kappa_counts(counts_data, r = r_total, method = "available", weight = "quadratic")

  # We compare only the estimates they have in common: Fleiss and BP
  common_estimates <- c("Fleiss", "Brennan-Prediger")

  # Compare estimates
  expect_equal(
    res_raw_np$estimates[common_estimates],
    res_counts_np$estimates,
    tolerance = 1e-9,
    info = "Estimates from counts 'np' and raw 'np' should match."
  )

  # Compare variance-covariance matrices
  expect_equal(
    res_raw_np$vcov[common_estimates, common_estimates],
    res_counts_np$vcov,
    tolerance = 1e-9,
    info = "V-Cov matrices from counts 'np' and raw 'np' should match."
  )
})


# Loop through the weight types to test ml and quadratic methods
for (w in weights_to_test) {
  test_that(paste("Counts methods match raw methods for '", w, "' weights"), {
    # 1. Calculate from raw data using 'ml' (our reference for Fleiss/BP)
    res_raw_ml <- kappa_raw(raw_data, method = "ml", weight = w, em_options = list(tol = 1e-12))

    # We compare against the Fleiss and BP coefficients from the raw calculation
    raw_fleiss_bp <- res_raw_ml$estimates[c("Fleiss", "Brennan-Prediger")]
    raw_vcov_fleiss_bp <- res_raw_ml$vcov[c("Fleiss", "Brennan-Prediger"), c("Fleiss", "Brennan-Prediger")]

    # 2. Calculate from counts data using 'ml'
    res_counts_ml <- kappa_counts(counts_data, r = r_total, method = "ml", weight = w, em_options = list(tol = 1e-12))

    # Compare raw 'ml' with counts 'ml'
    expect_equal(
      raw_fleiss_bp,
      res_counts_ml$estimates,
      tolerance = 1e-7,
      info = paste("Estimates from raw 'ml' and counts 'ml' should match for weight:", w)
    )
    expect_equal(
      raw_vcov_fleiss_bp,
      res_counts_ml$vcov,
      tolerance = 1e-7,
      info = paste("V-Cov from raw 'ml' and counts 'ml' should match for weight:", w)
    )

    # 3. If the weight is quadratic, also test the 'quadratic' method
    if (w == "quadratic") {
      res_counts_quad <- kappa_counts(counts_data, r = r_total, method = "quadratic", weight = "quadratic")

      # Compare counts 'ml' with counts 'quadratic'
      expect_equal(
        res_counts_ml$estimates,
        res_counts_quad$estimates,
        tolerance = 1e-7,
        info = "Estimates from counts 'ml' and counts 'quadratic' should match for quadratic weights."
      )
      expect_equal(
        res_counts_ml$vcov,
        res_counts_quad$vcov,
        tolerance = 1e-7,
        info = "V-Cov from counts 'ml' and counts 'quadratic' should match for quadratic weights."
      )
    }
  })
}

test_that("Counts 'np' matches raw 'np' for INCOMPLETE data", {
  # This is a critical test of the new C++ backend for counts data.
  # It verifies that analyzing incomplete data in raw or counts form
  # yields identical results for Fleiss/BP kappas.

  raw_data_incomplete <- dat.klein2018
  counts_data_incomplete <- to_counts_matrix(raw_data_incomplete)
  r_total_incomplete <- ncol(raw_data_incomplete) # r is not constant, but max r

  # Using "available" method (unweighted)
  res_raw <- kappa_raw(raw_data_incomplete, method = "available", weight = "quadratic")
  res_counts <- kappa_counts(
    counts_data_incomplete,
    r = r_total_incomplete,
    method = "available",
    weight = "quadratic"
  )

  common_ests <- c("Fleiss", "Brennan-Prediger")

  # Compare estimates
  expect_equal(
    res_raw$estimates[common_ests],
    res_counts$estimates,
    tolerance = 1e-9,
    info = "Fleiss/BP estimates from raw 'np' and counts 'np' should match on incomplete data."
  )

  # Compare variance-covariance matrices
  expect_equal(
    res_raw$vcov[common_ests, common_ests],
    res_counts$vcov,
    tolerance = 1e-9,
    info = "Fleiss/BP V-Cov from raw 'np' and counts 'np' should match on incomplete data."
  )
})
