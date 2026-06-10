pairwise_only_3rater_2cat <- function(reps = 3L) {
  rows <- list()
  k <- 1L
  for (rep in seq_len(reps)) {
    for (a in 0:1) {
      for (b in 0:1) {
        rows[[k]] <- c(a, b, NA_integer_); k <- k + 1L
        rows[[k]] <- c(a, NA_integer_, b); k <- k + 1L
        rows[[k]] <- c(NA_integer_, a, b); k <- k + 1L
      }
    }
  }
  do.call(rbind, rows)
}

chain_missing_pair_3rater_2cat <- function() {
  matrix(c(
    0, 0, NA,
    0, 0, NA,
    0, 0, NA,
    1, 1, NA,
    1, 1, NA,
    0, 1, NA,
    NA, 0, 0,
    NA, 0, 0,
    NA, 1, 1,
    NA, 1, 1,
    NA, 1, 0,
    NA, 0, 1
  ), ncol = 3, byrow = TRUE)
}

test_that("cat_fiml allows benign nuisance non-identification", {
  x <- pairwise_only_3rater_2cat()
  fit <- misskappa::kappa(x, estimator = "cat_fiml", weight = "nominal")

  expect_s3_class(fit, "misskappa_estimate")
  expect_named(fit$estimates, c("Conger", "Fleiss", "Brennan-Prediger"))
  expect_true(all(is.finite(fit$estimates)))
  expect_true(all(is.finite(fit$vcov)))
  expect_equal(crossprod(fit$psi) / nrow(x)^2, fit$vcov, tolerance = 1e-10)
})

test_that("raw FIML backend errors on coefficient non-identification", {
  expect_error(
    misskappa:::estimate_kappa_raw(chain_missing_pair_3rater_2cat(),
                                   method = "fiml", weight = "identity"),
    "not identified"
  )
})

sparse_benign_3rater_3cat <- function() {
  # All rater pairs co-observed, but n = 12 cannot support the saturated 3^3
  # joint: the former hard gradient gate fired here; now it is a warning plus
  # a null_frac diagnostic.
  x <- matrix(c(
    1, 1, 1,
    2, 2, NA,
    3, 3, 3,
    1, 1, 2,
    NA, 3, 2,
    3, 2, 3,
    1, 2, 1,
    2, 2, 3,
    3, 3, 3,
    1, 1, 1,
    2, 1, 2,
    3, 3, 2
  ), ncol = 3, byrow = TRUE)
  storage.mode(x) <- "integer"
  x
}

test_that("sparse benign cat_fiml succeeds with a null_frac warning", {
  x <- sparse_benign_3rater_3cat()
  expect_warning(
    fit <- misskappa::kappa(x, estimator = "cat_fiml", weight = "quadratic"),
    "not uniquely identified"
  )
  expect_true(all(is.finite(fit$estimates)))
  expect_true(all(is.finite(fit$vcov)))
  expect_named(fit$null_frac, c("Conger", "Fleiss", "Brennan-Prediger"))
  expect_true(all(fit$null_frac >= 0 & fit$null_frac <= 1))
  expect_gt(max(fit$null_frac), 0.01)
})

test_that("flatten quiets the warning and selects a nearby posterior mode", {
  x <- sparse_benign_3rater_3cat()
  strict <- suppressWarnings(
    misskappa::kappa(x, estimator = "cat_fiml", weight = "quadratic")
  )
  expect_no_warning(
    flat <- misskappa::kappa(x, estimator = "cat_fiml", weight = "quadratic",
                             em_options = list(flatten = 0.1))
  )
  expect_true(all(is.finite(flat$estimates)))
  # Total pseudo-mass 0.1 against n = 12 shrinks toward uniform with weight
  # ~ 0.1 / 12.1; the estimates stay close to the strict-ML face.
  expect_lt(max(abs(flat$estimates - strict$estimates)), 0.1)
  # The flattened mode is start-independent.
  flat2 <- misskappa::kappa(x, estimator = "cat_fiml", weight = "quadratic",
                            em_options = list(flatten = 0.1, start_alpha = 1))
  expect_lt(max(abs(flat$estimates - flat2$estimates)), 1e-5)
})

test_that("complete-data cat_fiml carries a zero null_frac diagnostic", {
  set.seed(7)
  x <- matrix(sample(1:3, 60, TRUE), ncol = 3)
  fit <- misskappa::kappa(x, estimator = "cat_fiml", weight = "nominal")
  expect_true(all(fit$null_frac < 1e-6))
})
