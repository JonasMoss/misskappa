# The internal joint covariance engine, joint_vcov(). It is unexported; the
# public agreement_test() exercises the common paths, so here we drive the
# validation branches directly (the package namespace makes them reachable
# from testthat).

make_pair <- function(n = 150L, seed = 11L) {
  set.seed(seed)
  list(a = agreement(matrix(sample.int(4, n * 3, TRUE), n, 3), estimator = "pairwise"),
       b = agreement(matrix(sample.int(4, n * 3, TRUE), n, 3), estimator = "pairwise"))
}

test_that("joint_vcov validates inputs and defaults the block names", {
  fits <- make_pair()

  # Unnamed fits fall back to fit1./fit2. prefixes.
  p <- length(coef(fits$a))
  V <- joint_vcov(fits$a, fits$b)
  expect_true(all(grepl("^fit[12]\\.", rownames(V))))
  expect_equal(dim(V), c(2L * p, 2L * p))

  expect_error(joint_vcov(fits$a), "at least two")
  expect_error(joint_vcov(fits$a, 42), "misskappa_estimate")

  # A fit stripped of its influence functions is rejected.
  bad <- fits$b
  bad$psi <- NULL
  expect_error(joint_vcov(fits$a, bad), "influence functions")

  # Mismatched subject counts are rejected.
  set.seed(99)
  small <- agreement(matrix(sample.int(4, 80 * 3, TRUE), 80, 3), estimator = "pairwise")
  expect_error(joint_vcov(fits$a, small), "same number of subjects")
})
