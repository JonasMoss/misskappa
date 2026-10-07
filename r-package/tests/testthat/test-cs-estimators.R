# Compound-symmetry estimators on a connected-but-incomplete co-observation
# graph (the synthetic planned-incomplete dat.tenhove2025 design), plus the
# influence-function contract for the CS / counts-CS estimators.

test_that("dat.tenhove2025 is a connected-but-incomplete design", {
  x <- dat.tenhove2025
  expect_equal(dim(x), c(60L, 6L))
  # Every subject has exactly three of the six raters.
  expect_true(all(rowSums(is.finite(x)) == 3L))
  co <- crossprod(is.finite(x))
  diag(co) <- 0
  expect_true(any(co[upper.tri(co)] == 0))         # some pair never co-observed
  expect_true(all(colSums(is.finite(x)) > 0))      # every rater appears
})

test_that("saturated coefficients are not identified on an incomplete graph", {
  for (est in c("pairwise", "nt_fiml")) {
    expect_error(kappa(dat.tenhove2025, estimator = est), "not identified")
  }
})

test_that("compound-symmetry estimators fit the incomplete design", {
  x <- dat.tenhove2025
  fits <- list(
    nt_fiml_cs = kappa(x, estimator = "nt_fiml_cs", values = 1:5),
    pairwise_cs = .kappa_pairwise_cs(x, values = 1:5)  # internal since 2.0.0
  )
  for (fit in fits) {
    expect_s3_class(fit, "misskappa_estimate")
    expect_named(fit$estimates, c("Conger", "Fleiss", "Brennan-Prediger"))
    expect_true(all(is.finite(fit$estimates)))
    expect_true(all(is.finite(diag(fit$vcov))))
    # Influence-function contract: vcov == crossprod(psi) / n^2, and the
    # registered stats::influence() accessor returns the same matrix.
    expect_equal(unname(crossprod(fit$psi) / nrow(x)^2), unname(fit$vcov),
                 tolerance = 1e-8)
    expect_equal(stats::influence(fit), fit$psi)
  }
})

test_that("nt_fiml_cs and pairwise_cs share their MCAR limit", {
  x <- dat.tenhove2025
  a <- kappa(x, estimator = "nt_fiml_cs", values = 1:5)$estimates
  b <- .kappa_pairwise_cs(x, values = 1:5)$estimates
  expect_lt(max(abs(a - b)), 0.05)
})

test_that("counts nt_fiml_cs (exchangeable CS) carries the IF identity", {
  fit <- kappa_counts(dat.fleiss1971, estimator = "nt_fiml_cs")
  expect_named(fit$estimates, c("Fleiss", "Brennan-Prediger"))
  expect_true(all(is.finite(fit$estimates)))
  expect_equal(unname(crossprod(fit$psi) / nrow(dat.fleiss1971)^2),
               unname(fit$vcov), tolerance = 1e-8)
})
