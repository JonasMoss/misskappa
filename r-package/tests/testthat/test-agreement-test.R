test_that("kappa_test handles independent, paired, one-sample, and G-way", {
  set.seed(1)
  x1 <- matrix(sample.int(4, 200 * 3, TRUE), 200, 3)
  x2 <- matrix(sample.int(4, 200 * 3, TRUE), 200, 3)
  k1 <- kappa(x1, estimator = "ipw")
  k2 <- kappa(x2, estimator = "ipw")

  ind <- kappa_test(a = k1, b = k2, coef = "Conger", paired = FALSE)
  expect_s3_class(ind, "htest")
  expect_named(ind$statistic, "X-squared")
  expect_equal(unname(ind$parameter), 1)
  expect_named(ind$estimate, c("a", "b"))
  expect_true(is.finite(ind$p.value))

  one <- kappa_test(k1, coef = "Conger", value = 0)
  expect_s3_class(one, "htest")
  expect_equal(unname(one$parameter), 1)

  x3 <- matrix(sample.int(4, 200 * 3, TRUE), 200, 3)
  g3 <- kappa_test(k1, k2, kappa(x3, estimator = "ipw"),
                   coef = "Conger", paired = FALSE)
  expect_equal(unname(g3$parameter), 2)
})

test_that("alpha_test(paired) reproduces the internal joint_vcov contrast", {
  set.seed(2)
  n <- 300L
  f <- rnorm(n)
  X <- sapply(1:10, function(j) 0.6 * f + rnorm(n, sd = 0.8))
  fa <- alpha(X[, 1:5], estimator = "pairwise")     # same n subjects, two item sets
  fb <- alpha(X[, 6:10], estimator = "pairwise")

  kt <- alpha_test(fa, fb, paired = TRUE)
  V  <- joint_vcov(a = fa, b = fb)
  d  <- coef(fa)[["alpha"]] - coef(fb)[["alpha"]]
  v  <- V[1, 1] + V[2, 2] - 2 * V[1, 2]
  expect_equal(unname(kt$statistic), d^2 / v, tolerance = 1e-10)

  # paired uses the cross-term; independent drops it -> different statistic
  ki <- alpha_test(fa, fb, paired = FALSE)
  expect_false(isTRUE(all.equal(unname(kt$statistic), unname(ki$statistic))))
})

test_that("test front doors validate inputs", {
  set.seed(3)
  x <- matrix(sample.int(4, 120 * 3, TRUE), 120, 3)
  k <- kappa(x, estimator = "ipw")
  expect_error(kappa_test(k, k, coef = "Nope"), "not found")
})

test_that("equality tests reject non-fits, empty calls, and broken pairing", {
  set.seed(4)
  k1 <- kappa(matrix(sample.int(4, 200 * 3, TRUE), 200, 3), estimator = "pairwise")
  k2 <- kappa(matrix(sample.int(4, 200 * 3, TRUE), 200, 3), estimator = "pairwise")

  # non-'misskappa_estimate' inputs and empty calls
  expect_error(kappa_test(1, 2), "misskappa_estimate")
  expect_error(kappa_test(), "at least one fit")
  expect_error(alpha_test(), "at least one fit")

  # paired requires influence functions
  k2_nopsi <- k2
  k2_nopsi$psi <- NULL
  expect_error(kappa_test(k1, k2_nopsi, paired = TRUE), "influence functions")

  # paired fits must share the subject count
  small <- kappa(matrix(sample.int(4, 120 * 3, TRUE), 120, 3), estimator = "pairwise")
  expect_error(kappa_test(k1, small, paired = TRUE), "same number of subjects")

  # paired fits must be row-aligned when subject ids are present
  k1_id <- k1
  k2_id <- k2
  rownames(k1_id$psi) <- as.character(seq_len(nrow(k1_id$psi)))
  rownames(k2_id$psi) <- as.character(nrow(k2_id$psi) + seq_len(nrow(k2_id$psi)))
  expect_error(kappa_test(k1_id, k2_id, paired = TRUE), "row-aligned")
})
