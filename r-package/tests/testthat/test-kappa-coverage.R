# Targeted coverage for the kappa()/alpha() dispatch and the internal estimator
# wrappers: argument validation, the values= scoring path, the cat_fiml route,
# the g-wise and vector dispatch branches, and the confint() transforms.

cat_mat <- function(n = 20L, R = 3L, k = 3L, seed = 1L) {
  set.seed(seed)
  matrix(sample.int(k, n * R, TRUE), n, R)
}

# --- public kappa() dispatch validation ---

test_that("kappa() validates g and the estimator/weight combination", {
  x <- cat_mat()
  expect_error(kappa(x, g = 1), "integer >= 2")
  expect_error(kappa(x, g = 2.5), "integer >= 2")
  expect_error(kappa(x, estimator = "pairwise", weight = "linear"),
               "requires weight")
  expect_error(kappa(x, estimator = "nt_fiml", weight = "nominal"),
               "requires weight")
})

test_that("kappa() routes bad raw input through the raw/score validators", {
  expect_error(kappa(list(1, 2), estimator = "ipw"), "matrix or data frame")
  expect_error(kappa(matrix("a", 4, 3), estimator = "ipw"), "must be numeric")
  expect_error(kappa(list(1, 2), estimator = "pairwise"), "matrix or data frame")
  expect_error(kappa(matrix("a", 4, 3), estimator = "pairwise"), "must be numeric")
})

test_that("kappa() maps integer codes through values=", {
  x <- cat_mat(k = 3L)
  fit <- kappa(x, estimator = "pairwise", values = c(0, 1, 2))
  expect_s3_class(fit, "misskappa_estimate")
  expect_error(kappa(x, estimator = "pairwise", values = c(0, 1)),
               "unique observed categories")
})

# --- g-wise and vector dispatch ---

test_that("g-wise dispatch honours max_chance_tuples and rejects bad weights", {
  x <- cat_mat(R = 4L)
  fit <- kappa(x, estimator = "ipw", weight = "nominal", g = 3,
               max_chance_tuples = 1e5)
  expect_s3_class(fit, "misskappa_estimate")
  expect_true(!is.null(fit$g) && fit$g == 3L)
  expect_output(print(fit), "g=3")          # print header for g > 2
  expect_error(kappa(x, estimator = "ipw", weight = "ordinal", g = 3),
               "should be one of")
})

test_that("vector dispatch maps weights to component losses", {
  arr <- array(rbinom(2 * 3 * 4, 1, 0.5), dim = c(4, 3, 2))
  expect_s3_class(kappa(arr, estimator = "pairwise", weight = "linear"),
                  "misskappa_estimate")
  expect_error(kappa(arr, estimator = "pairwise", weight = "ordinal"),
               "should be one of")
  expect_error(kappa(arr, estimator = "cat_fiml"), "support estimator")
})

# --- confint transforms ---

test_that("confint() supports the Fisher transform and parm subsetting", {
  fit <- kappa(cat_mat(), estimator = "ipw")
  z <- confint(fit, transform = "fisher")
  expect_true(all(z > -1 & z < 1))
  one <- confint(fit, parm = "Conger")
  expect_equal(rownames(one), "Conger")

  # Boundary estimate (perfect agreement) -> NA Fisher limits with a warning.
  perfect <- cbind(c(1, 2, 1, 2, 1, 2), c(1, 2, 1, 2, 1, 2))
  pf <- kappa(perfect, estimator = "pairwise", values = c(0, 1))
  expect_warning(confint(pf, transform = "fisher"), "Fisher transform requires")
})

# --- alpha() dispatch and scoring ---

test_that("alpha() validates input and maps values=", {
  expect_error(alpha(1:5), "matrix or data frame")
  expect_error(alpha(matrix("a", 5, 3)), "must be numeric")
  expect_error(alpha(matrix(1, 5, 1)), "at least two items")

  x <- cat_mat(R = 4L, k = 3L)
  fit <- alpha(x, estimator = "pairwise", values = c(0, 1, 2))
  expect_s3_class(fit, "misskappa_estimate")
  expect_error(alpha(x, estimator = "pairwise", values = c(0, 1)),
               "unique observed categories")
  expect_error(alpha(matrix(NA_real_, 5, 3), estimator = "pairwise",
                     values = c(1, 2, 3)),
               "every item must be observed")
})

test_that("alpha(estimator = 'cat_fiml') runs and validates", {
  x <- cat_mat(R = 4L, k = 3L)
  fit <- alpha(x, estimator = "cat_fiml")
  expect_s3_class(fit, "misskappa_estimate")
  expect_equal(fit$method, "cat_fiml")
  expect_error(alpha_cat_fiml(list(1, 2)), "matrix or data frame")
  expect_error(alpha_cat_fiml(matrix("a", 5, 3)), "must be numeric")
})

# --- public and internal weighting schemes exercise their glue branches ---

test_that("public weighting schemes run through the raw and counts glue", {
  weights <- c("nominal", "linear", "quadratic")
  x <- cat_mat(R = 3L, k = 4L)
  for (w in weights) {
    expect_s3_class(kappa(x, estimator = "ipw", weight = w), "misskappa_estimate")
    expect_s3_class(kappa_counts(dat.fleiss1971, estimator = "pairwise", weight = w),
                    "misskappa_estimate")
    expect_s3_class(kappa_counts(dat.fleiss1971, estimator = "cat_fiml", weight = w),
                    "misskappa_estimate")
  }

  # values= length is validated in the C++ glue (raw and counts variants).
  expect_error(kappa(x, estimator = "ipw", weight = "linear", values = c(1, 2)),
               "Length of 'values'")
  expect_error(
    kappa_counts(dat.fleiss1971, weight = "linear", values = c(1, 2)),
               "Length of 'values'")
})

test_that("legacy weighting schemes are internal-only", {
  legacy_weights <- c("ordinal", "radical", "ratio", "circular", "bipolar")
  x <- cat_mat(R = 3L, k = 4L)
  counts <- as.matrix(dat.fleiss1971)
  storage.mode(counts) <- "integer"

  for (w in legacy_weights) {
    expect_error(kappa(x, estimator = "ipw", weight = w), "should be one of")
    expect_error(kappa_counts(counts, estimator = "pairwise", weight = w),
                 "should be one of")

    expect_s3_class(estimate_kappa_raw(x, method = "ipw", weight = w),
                    "misskappa_estimate")
    expect_type(rcpp_kappa_counts(counts, weight_type = w, values = NULL),
                "list")
    expect_type(rcpp_kappa_fiml_counts(
      counts, weight_type = w, values = NULL, r_total = 6L,
      em_options = list()), "list")
  }
})

# --- internal estimator wrappers: validation guards ---

test_that("internal estimator wrappers validate their inputs", {
  expect_error(kappa_quadratic(list(1), 1:3), "matrix or data frame")
  expect_error(kappa_quadratic(matrix(1, 4, 3), "a"), "'values' must be numeric")

  expect_error(kappa_counts(list(1)), "matrix or data frame")
  expect_error(kappa_counts(matrix("a", 4, 3)), "must be numeric")

  expect_error(kappa_continuous(list(1)), "matrix or data frame")
  expect_error(kappa_continuous(matrix("a", 4, 2)), "must be numeric")

  expect_error(kappa_vector(matrix(1, 2, 2)), "subjects-by-raters-by-features")
  expect_error(kappa_vector(array(1, c(2, 1, 1))), "two raters")
  expect_error(kappa_vector(array("a", c(2, 2, 1))), "must be numeric")
  expect_error(kappa_vector(array(1, c(3, 2, 2)), feature_weights = c(-1, 1)),
               "non-negative")

  expect_error(kappa_gwise(list(1)), "matrix or data frame")
  expect_error(kappa_gwise(matrix("a", 4, 3)), "must be numeric")
  expect_error(kappa_gwise(matrix(1:12, 4, 3), g = 99), "between 2 and ncol")
  expect_error(kappa_gwise(matrix(1:12, 4, 3), max_chance_tuples = -5),
               "positive integer")
  expect_error(kappa_gwise(matrix(1:12, 4, 3), distance = "absolute",
                           method = "fiml"),
               "categorical g-wise")
})
