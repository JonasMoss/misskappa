# Tests for the normal-FIML quadratic (Conger / Fleiss) kappa.

make_ratings <- function(n, R, rho = 0.5, miss = 0.15, seed = 1L,
                         mu = NULL) {
  set.seed(seed)
  S <- matrix(rho, R, R); diag(S) <- 1
  X <- matrix(stats::rnorm(n * R), n, R) %*% chol(S)
  if (!is.null(mu)) X <- sweep(X, 2L, mu, "+")
  if (miss > 0) X[matrix(stats::runif(n * R) < miss, n, R)] <- NA
  colnames(X) <- paste0("r", seq_len(R))
  X
}

# Reference Conger / Fleiss from a mean vector and covariance matrix.
ref_kappa <- function(mu, S) {
  R <- length(mu)
  t1 <- sum(S); t2 <- sum(diag(S)); t3 <- sum((mu - mean(mu))^2)
  c(Conger = (t1 - t2) / ((R - 1) * t2 + R * t3),
    Fleiss = (t1 - t2 - t3) / ((R - 1) * (t2 + t3)))
}

test_that("complete-data point matches the MLE-covariance plug-in", {
  X <- make_ratings(400L, 4L, miss = 0, seed = 11L, mu = c(0, 0.4, -0.3, 0.2))
  fit <- kappa_quadratic_fiml(X)
  n <- nrow(X)
  Smle <- stats::cov(X) * (n - 1) / n
  manual <- ref_kappa(colMeans(X), Smle)
  expect_equal(coef(fit)[c("Conger", "Fleiss")], manual, tolerance = 1e-8)
  expect_true(fit$moments$converged)
})

test_that("complete-data point matches kappa_quadratic(empirical)", {
  set.seed(21L)
  # scored 1..5 categorical ratings so kappa_quadratic's `values` is natural
  X <- matrix(sample(1:5, 350L * 4L, replace = TRUE), 350L, 4L)
  storage.mode(X) <- "double"
  fit <- kappa_quadratic_fiml(X)
  ref <- kappa_quadratic(X, values = 1:5)
  expect_equal(unname(coef(fit)["Conger"]), unname(coef(ref)["Conger"]),
               tolerance = 1e-8)
  expect_equal(unname(coef(fit)["Fleiss"]), unname(coef(ref)["Fleiss"]),
               tolerance = 1e-8)
})

test_that("sandwich vcov satisfies the crossprod(psi)/n^2 influence contract", {
  X <- make_ratings(500L, 5L, miss = 0.2, seed = 2L)
  fit <- kappa_quadratic_fiml(X)
  psi <- fit$psi
  expect_true(is.matrix(psi))
  expect_equal(dim(psi), c(nrow(X), 2L))
  expect_equal(colnames(psi), c("Conger", "Fleiss"))
  n <- nrow(X)
  expect_equal(unname(crossprod(psi) / (n * n)), unname(vcov(fit)),
               tolerance = 1e-10)
})

test_that("analytic gradient G matches a numeric derivative wrt theta", {
  X <- make_ratings(300L, 4L, miss = 0.2, seed = 4L)
  em <- misskappa:::.amc_em(X, misskappa:::.amc_patterns(is.finite(X)))
  p <- ncol(X); pstar <- p * (p + 1L) / 2L
  grad <- misskappa:::.kqf_grad(em$mu, em$Sigma)

  kap_of_theta <- function(theta) {
    pr <- misskappa:::.amc_unpack(theta, p, pstar)
    ref_kappa(pr$mu, pr$Sigma)
  }
  theta <- c(em$mu, misskappa:::.amc_vech(em$Sigma))
  h <- 1e-6
  Gnum <- matrix(0, 2L, length(theta))
  for (k in seq_along(theta)) {
    tp <- theta; tp[k] <- tp[k] + h
    tm <- theta; tm[k] <- tm[k] - h
    Gnum[, k] <- (kap_of_theta(tp) - kap_of_theta(tm)) / (2 * h)
  }
  # .amc_unpack places off-diagonal vech entries symmetrically, so the numeric
  # derivative already carries the weight-2 off-diagonal convention that the
  # analytic gradient uses (matching .amc_score_matrix). Compare G directly.
  expect_equal(unname(grad$G), Gnum, tolerance = 1e-6)
})

test_that("standard errors are stable across finite-difference steps", {
  X <- make_ratings(450L, 4L, rho = 0.4, miss = 0.2, seed = 8L)
  hs <- c(1e-4, 1e-5, 1e-6)
  ses <- vapply(hs, function(h) {
    fit <- kappa_quadratic_fiml(X, em_options = list(fd_h = h))
    sqrt(vcov(fit)[1, 1])
  }, numeric(1))
  rel_span <- diff(range(ses)) / mean(ses)
  expect_lt(rel_span, 1e-4)
})

test_that("kappa_continuous(method = 'fiml') routes here and validates weight", {
  X <- make_ratings(300L, 4L, miss = 0.15, seed = 9L)
  fit <- kappa_continuous(X, method = "fiml", weight = "quadratic")
  expect_equal(fit$method, "quadratic-fiml")
  expect_equal(names(coef(fit)), c("Conger", "Fleiss"))
  expect_error(
    kappa_continuous(X, method = "fiml", weight = "linear"),
    "quadratic"
  )
})

test_that("fits plug into joint_vcov / wald_test", {
  X <- make_ratings(350L, 5L, miss = 0.2, seed = 5L)
  fitA <- kappa_quadratic_fiml(X)
  fitB <- kappa_quadratic_fiml(X[, 1:4])
  V <- joint_vcov(full = fitA, sub = fitB)
  expect_equal(dim(V), c(4L, 4L))
  expect_true(isSymmetric(V))
  # Conger == Fleiss within a single fit (same subjects) is a valid contrast.
  wt <- wald_test(fitA, contrast = c("Conger" = 1, "Fleiss" = -1))
  expect_s3_class(wt, "htest")
  expect_true(is.finite(wt$p.value))
})

test_that("saturated EM covariance matches magmaan's oracle", {
  skip_if_not_installed("magmaan")
  X <- make_ratings(500L, 5L, miss = 0.2, seed = 6L)
  fit <- kappa_quadratic_fiml(X)
  mask <- is.finite(X); storage.mode(mask) <- "logical"
  mag <- magmaan::magmaan_core$estimate_saturated_em_moments(
    list(X = X, mask = mask))
  expect_equal(max(abs(fit$moments$Sigma - mag$cov[[1L]])), 0, tolerance = 1e-5)
})

test_that("sandwich vcov matches magmaan's saturated-moment sandwich ACOV", {
  # magmaan exposes the saturated-EM sandwich ingredients H, J and
  # acov = H^-1 J H^-1 for theta = (mu, vech Sigma) in the same column-major
  # ordering. The functional vcov G acov G^T must reproduce the per-subject
  # influence-function vcov, validated through magmaan's independent FIML path.
  skip_if_not_installed("magmaan")
  X <- make_ratings(500L, 4L, miss = 0.2, seed = 12L,
                    mu = c(-0.4, -0.1, 0.2, 0.5))
  fit <- kappa_quadratic_fiml(X)
  G <- misskappa:::.kqf_grad(fit$moments$mu, fit$moments$Sigma)$G
  mask <- is.finite(X); storage.mode(mask) <- "logical"
  acov <- as.matrix(magmaan::magmaan_core$estimate_saturated_em_moments(
    list(X = X, mask = mask))$acov)
  V_oracle <- G %*% acov %*% t(G)
  expect_equal(unname(vcov(fit)), unname(V_oracle), tolerance = 1e-6)
})
