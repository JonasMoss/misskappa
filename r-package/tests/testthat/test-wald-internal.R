# The internal Wald engine: joint_vcov(), wald_test(), normalise_contrast().
# These are unexported -- the public kappa_test()/alpha_test() exercise the
# common paths, so here we drive the contrast-handling and validation branches
# directly (the package namespace makes them reachable from testthat).

make_pair <- function(n = 150L, seed = 11L) {
  set.seed(seed)
  list(a = kappa(matrix(sample.int(4, n * 3, TRUE), n, 3), estimator = "pairwise"),
       b = kappa(matrix(sample.int(4, n * 3, TRUE), n, 3), estimator = "pairwise"))
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
  small <- kappa(matrix(sample.int(4, 80 * 3, TRUE), 80, 3), estimator = "pairwise")
  expect_error(joint_vcov(fits$a, small), "same number of subjects")
})

test_that("wald_test validates fits, null values, and contrast type", {
  fits <- make_pair()
  cn <- names(coef(fits$a))

  expect_error(wald_test(contrast = cn[1]), "at least one fit")
  expect_error(wald_test(42, contrast = cn[1]), "misskappa_estimate")
  expect_error(wald_test(fits$a, contrast = cn[1], value = Inf), "finite")
  expect_error(wald_test(fits$a, contrast = cn[1], value = c(0, 0)),
               "scalar or have one entry")
  expect_error(wald_test(fits$a, contrast = TRUE), "numeric or character")
})

test_that("normalise_contrast accepts unnamed vectors and matrix contrasts", {
  fits <- make_pair()
  cn <- names(coef(fits$a))
  p <- length(cn)

  # Unnamed numeric vector of the right length -> one contrast row.
  w_vec <- wald_test(fits$a, contrast = rep(1, p) / p)
  expect_s3_class(w_vec, "htest")
  expect_equal(unname(w_vec$parameter), 1)

  # Matrix with (a subset of) column names -> padded to full width.
  Lc <- matrix(c(1, -1), nrow = 1, dimnames = list("diff", cn[1:2]))
  w_named <- wald_test(fits$a, contrast = Lc)
  expect_s3_class(w_named, "htest")
  expect_equal(rownames(w_named$estimate), NULL)
  expect_named(w_named$estimate, "diff")

  # Unnamed matrix of full width -> column names assigned, rows auto-named.
  Lu <- diag(p)[1:2, , drop = FALSE]
  w_unnamed <- wald_test(fits$a, contrast = Lu)
  expect_s3_class(w_unnamed, "htest")
  expect_equal(unname(w_unnamed$parameter), 2)
  expect_match(names(w_unnamed$estimate)[1], "^contrast")

  # Unknown coefficient names are rejected.
  expect_error(wald_test(fits$a, contrast = c(nope = 1)), "Unknown coefficient")
})
