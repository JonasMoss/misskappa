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

# ---- co-observation identifiability guard ---------------------------------
# A saturated coefficient (alpha; Conger/Cohen/BP kappa) needs every column pair
# co-observed by >= 1 subject -- a COMPLETE co-observation graph, not merely a
# connected one. See R/check_identifiable.R and the alpha-missing identifiability
# note.

# 3-column NA mask whose (1,3) pair is never jointly observed: the first block
# sees columns 1-2, the second sees columns 2-3. Connected but incomplete.
.incomplete_mask <- function(n = 80L) {
  blk <- seq_len(n) <= n %/% 2L
  mask <- matrix(FALSE, n, 3L)
  mask[blk, 3L]  <- TRUE                      # block 1: columns 1, 2 observed
  mask[!blk, 1L] <- TRUE                      # block 2: columns 2, 3 observed
  mask
}

test_that(".check_pattern_identifiable flags an incomplete co-observation graph", {
  obs <- !.incomplete_mask(60L)
  expect_error(
    .check_pattern_identifiable(obs, unit = "item", coefficient = "coefficient alpha"),
    "never jointly observed")

  full <- matrix(TRUE, 30L, 3L)              # complete graph passes
  N <- .check_pattern_identifiable(full, unit = "item")
  expect_equal(unname(N[1L, 3L]), 30)
})

test_that(".check_pattern_identifiable warns on a single-overlap pair", {
  obs <- matrix(TRUE, 20L, 3L)
  obs[-1L, 3L] <- FALSE                       # column 3 seen by subject 1 only
  expect_warning(.check_pattern_identifiable(obs, unit = "item"),
                 "only one subject")
})

test_that(".check_pattern_identifiable enforces the exchangeable relaxation", {
  obs <- matrix(TRUE, 10L, 4L)
  obs[1L, 2:4] <- FALSE                       # subject 1 has one observed rater
  expect_error(
    .check_pattern_identifiable(obs, unit = "rater", require = "each_subject_2"),
    "at least two observed")
  expect_silent(
    .check_pattern_identifiable(matrix(TRUE, 6L, 3L), unit = "rater",
                                require = "each_subject_2"))
})

test_that(".check_pattern_identifiable supports g-wise tuple arity", {
  obs <- rbind(
    c(TRUE, TRUE, FALSE),
    c(TRUE, TRUE, FALSE),
    c(TRUE, FALSE, TRUE),
    c(TRUE, FALSE, TRUE),
    c(FALSE, TRUE, TRUE),
    c(FALSE, TRUE, TRUE)
  )                                           # all pairs appear, no 1-2-3 row

  expect_silent(.check_pattern_identifiable(obs, unit = "rater", arity = 2L))
  expect_error(
    .check_pattern_identifiable(obs, unit = "rater", coefficient = "3-wise kappa",
                                arity = 3L),
    "3-tuple")
})

test_that("alpha() errors on a never-co-observed item pair for every estimator", {
  set.seed(11)
  mask <- .incomplete_mask(80L)
  Xc <- matrix(rnorm(80L * 3L), 80L, 3L);                 Xc[mask] <- NA_real_
  Xk <- matrix(sample(1:3, 80L * 3L, replace = TRUE), 80L, 3L); Xk[mask] <- NA_integer_
  expect_error(alpha(Xc, estimator = "pairwise"), "never jointly observed")
  expect_error(alpha(Xc, estimator = "nt_fiml"),  "never jointly observed")
  expect_error(alpha(Xk, estimator = "cat_fiml"), "never jointly observed")
})

test_that("kappa() errors on a never-co-observed rater pair for every estimator", {
  set.seed(12)
  mask <- .incomplete_mask(80L)
  Xc <- matrix(rnorm(80L * 3L), 80L, 3L);                 Xc[mask] <- NA_real_
  Xk <- matrix(sample(1:3, 80L * 3L, replace = TRUE), 80L, 3L); Xk[mask] <- NA_integer_
  expect_error(kappa(Xk, estimator = "ipw"),      "never jointly observed")
  expect_error(kappa(Xk, estimator = "cat_fiml"), "never jointly observed")
  expect_error(kappa(Xc, estimator = "pairwise", weight = "quadratic"),
               "never jointly observed")
  expect_error(kappa(Xc, estimator = "nt_fiml", weight = "quadratic"),
               "never jointly observed")
})

test_that("kappa() g-wise paths require every requested rater tuple", {
  x <- matrix(c(
    0, 0, NA,
    0, 1, NA,
    0, NA, 0,
    1, NA, 0,
    NA, 0, 0,
    NA, 0, 1,
    NA, 1, 1,
    1, 0, NA
  ), ncol = 3, byrow = TRUE)

  expect_silent(kappa(x, estimator = "ipw", weight = "nominal", g = 2L))
  expect_error(kappa(x, estimator = "ipw", weight = "nominal", g = 3L),
               "3-tuple")
  expect_error(kappa(x, estimator = "cat_fiml", weight = "nominal", g = 3L),
               "3-tuple")
})

test_that("kappa_counts is the exchangeable representation and is not gated", {
  # Complete counts (Fleiss 1971): rater identity is discarded, so the
  # complete-graph condition is vacuous and no guard fires.
  expect_s3_class(kappa_counts(dat.fleiss1971, estimator = "fleiss_cuzick"),
                  "misskappa_estimate")
})
