# Input-validation guards for the continuous / quadratic-FIML / vector-quadratic
# estimator backends. These are internal functions, reachable bare from the
# package namespace (as the other tests already do).

test_that("alpha_continuous validates input and drops all-missing rows", {
  expect_error(alpha_continuous(list(1, 2)), "matrix or data frame")
  expect_error(alpha_continuous(matrix("a", 5, 3)), "must be numeric")
  expect_error(alpha_continuous(matrix(1, 5, 1)), "at least two items")

  set.seed(1)
  X <- matrix(rnorm(100), 25, 4)
  X_col <- X; X_col[, 1] <- NA_real_
  expect_error(alpha_continuous(X_col), "every item must be observed")

  X_row <- X; X_row[1, ] <- NA_real_         # dropped, then estimated on the rest
  expect_s3_class(alpha_continuous(X_row), "misskappa_estimate")
})

test_that("kappa_quadratic_fiml validates input", {
  expect_error(kappa_quadratic_fiml(list(1, 2)), "matrix or data frame")
  expect_error(kappa_quadratic_fiml(matrix("a", 5, 3)), "must be numeric")
  expect_error(kappa_quadratic_fiml(matrix(1, 5, 1)), "at least two raters")

  set.seed(2)
  X <- matrix(rnorm(100), 25, 4)
  X[, 1] <- NA_real_
  expect_error(kappa_quadratic_fiml(X), "every rater must be observed")
})

test_that("kappa_vector_quadratic validates the array, W, and observation", {
  expect_error(kappa_vector_quadratic(matrix(1, 2, 2)),
               "subjects-by-raters-by-features")
  expect_error(kappa_vector_quadratic(array(1, c(2, 1, 1))), "two raters")
  expect_error(kappa_vector_quadratic(array("a", c(2, 2, 1))), "must be numeric")
  expect_error(kappa_vector_quadratic(array(1, c(5, 2, 2)), W = matrix(1, 3, 3)),
               "features-by-features")

  set.seed(3)
  arr <- array(rnorm(5 * 2 * 2), c(5, 2, 2))
  arr[, 1, 1] <- NA_real_                    # one rater-feature cell never seen
  expect_error(kappa_vector_quadratic(arr, method = "pairwise"),
               "rater-feature cell")
  expect_error(kappa_vector_quadratic(arr, method = "nt_fiml"),
               "rater-feature cell")
})
