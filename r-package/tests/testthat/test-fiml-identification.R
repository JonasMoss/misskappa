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
    "not identified from the Louis information"
  )
})
