test_that("agreement_test handles independent, paired, one-sample, and G-way", {
  set.seed(1)
  x1 <- matrix(sample.int(4, 200 * 3, TRUE), 200, 3)
  x2 <- matrix(sample.int(4, 200 * 3, TRUE), 200, 3)
  k1 <- agreement(x1, estimator = "ipw")
  k2 <- agreement(x2, estimator = "ipw")

  ind <- agreement_test(a = k1, b = k2, coef = "Conger", paired = FALSE)
  expect_s3_class(ind, "htest")
  expect_named(ind$statistic, "X-squared")
  expect_equal(unname(ind$parameter), 1)
  expect_named(ind$estimate, c("a", "b"))
  expect_true(is.finite(ind$p.value))

  one <- agreement_test(k1, coef = "Conger", value = 0)
  expect_s3_class(one, "htest")
  expect_equal(unname(one$parameter), 1)

  x3 <- matrix(sample.int(4, 200 * 3, TRUE), 200, 3)
  g3 <- agreement_test(k1, k2, agreement(x3, estimator = "ipw"),
                   coef = "Conger", paired = FALSE)
  expect_equal(unname(g3$parameter), 2)
})

test_that("agreement_test(paired) reproduces the internal joint_vcov contrast", {
  set.seed(2)
  n <- 300L
  f <- rnorm(n)
  X <- sapply(1:6, function(j) 0.6 * f + rnorm(n, sd = 0.8))
  fa <- agreement(X[, 1:3], estimator = "pairwise")   # same n subjects, two rater sets
  fb <- agreement(X[, 4:6], estimator = "pairwise")

  kt <- agreement_test(fa, fb, paired = TRUE)
  V  <- research_joint_vcov(a = fa, b = fb)
  d  <- coef(fa)[["Conger"]] - coef(fb)[["Conger"]]
  v  <- V["a.Conger", "a.Conger"] + V["b.Conger", "b.Conger"] -
    2 * V["a.Conger", "b.Conger"]
  expect_equal(unname(kt$statistic), d^2 / v, tolerance = 1e-10)

  # paired uses the cross-term; independent drops it -> different statistic
  ki <- agreement_test(fa, fb, paired = FALSE)
  expect_false(isTRUE(all.equal(unname(kt$statistic), unname(ki$statistic))))
})

test_that("test front doors validate inputs", {
  set.seed(3)
  x <- matrix(sample.int(4, 120 * 3, TRUE), 120, 3)
  k <- agreement(x, estimator = "ipw")
  expect_error(agreement_test(k, k, coef = "Nope"), "not found")
})

test_that("equality tests reject non-fits, empty calls, and broken pairing", {
  set.seed(4)
  k1 <- agreement(matrix(sample.int(4, 200 * 3, TRUE), 200, 3), estimator = "pairwise")
  k2 <- agreement(matrix(sample.int(4, 200 * 3, TRUE), 200, 3), estimator = "pairwise")

  # non-'misskappa_estimate' inputs and empty calls
  expect_error(agreement_test(1, 2), "misskappa_estimate")
  expect_error(agreement_test(), "at least one fit")

  # paired requires influence functions
  k2_nopsi <- k2
  k2_nopsi$psi <- NULL
  expect_error(agreement_test(k1, k2_nopsi, paired = TRUE), "influence functions")

  # without subject ids, paired fits must share the subject count
  small <- agreement(matrix(sample.int(4, 120 * 3, TRUE), 120, 3), estimator = "pairwise")
  rownames(small$psi) <- NULL
  k1_noid <- k1
  rownames(k1_noid$psi) <- NULL
  expect_error(agreement_test(k1_noid, small, paired = TRUE),
               "same number of subjects")

  # disjoint subject ids share no influence, so paired equals independent
  k1_id <- k1
  k2_id <- k2
  rownames(k1_id$psi) <- as.character(seq_len(nrow(k1_id$psi)))
  rownames(k2_id$psi) <- as.character(nrow(k2_id$psi) + seq_len(nrow(k2_id$psi)))
  expect_equal(unname(agreement_test(k1_id, k2_id, paired = TRUE)$statistic),
               unname(agreement_test(k1_id, k2_id)$statistic))
})

test_that("agreement_test defaults: independent samples, first coefficient", {
  set.seed(5)
  k1 <- agreement(matrix(sample.int(4, 150 * 3, TRUE), 150, 3), estimator = "ipw")
  k2 <- agreement(matrix(sample.int(4, 150 * 3, TRUE), 150, 3), estimator = "ipw")
  def <- agreement_test(k1, k2)
  ind <- agreement_test(k1, k2, coef = "Conger", paired = FALSE)
  expect_equal(def$statistic, ind$statistic)
  expect_match(def$method, "Independent-sample")
  expect_equal(def$null.value, c("difference in Conger" = 0))
  expect_equal(def$alternative, "two.sided")

  c1 <- agreement_counts(dat.fleiss1971)
  one <- agreement_test(c1, value = 0.4)
  expect_equal(one$null.value, c(Fleiss = 0.4))
  expect_error(agreement_test(k1, k2, k1, value = 0.1), "must be 0")
})
