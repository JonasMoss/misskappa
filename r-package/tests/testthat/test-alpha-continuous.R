# Tests for the normal-FIML continuous-item coefficient alpha.

make_data <- function(n, p, rho = 0.4, miss = 0.15, seed = 1L) {
  set.seed(seed)
  S <- matrix(rho, p, p); diag(S) <- 1
  X <- matrix(stats::rnorm(n * p), n, p) %*% chol(S)
  if (miss > 0) X[matrix(stats::runif(n * p) < miss, n, p)] <- NA
  colnames(X) <- paste0("y", seq_len(p))
  X
}

test_that("complete-data alpha equals the MLE-covariance plug-in formula", {
  X <- make_data(300L, 5L, miss = 0)
  fit <- alpha_continuous(X)
  Smle <- stats::cov(X) * (nrow(X) - 1) / nrow(X)
  p <- ncol(X)
  manual <- (p / (p - 1)) * (1 - sum(diag(Smle)) / sum(Smle))
  expect_equal(unname(coef(fit)["alpha"]), manual, tolerance = 1e-8)
  expect_true(fit$moments$converged)
})

test_that("sandwich vcov satisfies the crossprod(psi)/n^2 influence contract", {
  X <- make_data(400L, 6L, miss = 0.2, seed = 2L)
  fit <- alpha_continuous(X)
  psi <- fit$psi
  expect_true(is.matrix(psi))
  expect_equal(dim(psi), c(nrow(X), 1L))
  expect_equal(colnames(psi), "alpha")
  n <- nrow(X)
  expect_equal(unname(crossprod(psi) / (n * n)), unname(vcov(fit)),
               tolerance = 1e-10)
})

test_that("standard errors ignore legacy fd_h option", {
  X <- make_data(450L, 5L, rho = 0.35, miss = 0.2, seed = 8L)
  hs <- c(1e-4, 1e-5, 1e-6)
  ses <- vapply(hs, function(h) {
    fit <- alpha_continuous(X, em_options = list(fd_h = h))
    sqrt(vcov(fit)[1, 1])
  }, numeric(1))
  rel_span <- diff(range(ses)) / mean(ses)
  expect_equal(rel_span, 0, tolerance = 1e-12)
})

test_that("analytic casewise scores match numerical differentiation", {
  X <- make_data(200L, 5L, miss = 0.2, seed = 4L)
  p <- ncol(X); pstar <- p * (p + 1L) / 2L
  vp <- misskappa:::.amc_vech_pos(p)
  pat <- misskappa:::.amc_patterns(is.finite(X))
  em <- misskappa:::.amc_em(X, pat)
  theta <- c(em$mu, misskappa:::.amc_vech(em$Sigma))
  S <- misskappa:::.amc_score_matrix(em$mu, em$Sigma, X, pat, vp, p, pstar)

  loglik_i <- function(th, xi) {
    pr <- misskappa:::.amc_unpack(th, p, pstar)
    o <- which(is.finite(xi))
    e <- xi[o] - pr$mu[o]
    So <- pr$Sigma[o, o, drop = FALSE]
    -0.5 * (length(o) * log(2 * pi) +
              as.numeric(determinant(So, TRUE)$modulus) +
              drop(t(e) %*% solve(So, e)))
  }
  h <- 1e-6
  for (r in c(1L, 17L, 88L, 150L)) {
    num <- vapply(seq_along(theta), function(k) {
      tp <- theta; tp[k] <- tp[k] + h
      tm <- theta; tm[k] <- tm[k] - h
      (loglik_i(tp, X[r, ]) - loglik_i(tm, X[r, ])) / (2 * h)
    }, numeric(1))
    expect_equal(S[r, ], num, tolerance = 1e-6,
                 info = sprintf("row %d", r))
  }
})

test_that("alpha_continuous fits plug into joint_vcov / wald_test", {
  X1 <- make_data(350L, 5L, miss = 0.2, seed = 5L)
  # A second item subset on the SAME subjects (drop last item).
  fitA <- alpha_continuous(X1)
  fitB <- alpha_continuous(X1[, 1:4])
  V <- joint_vcov(full = fitA, sub = fitB)
  expect_equal(dim(V), c(2L, 2L))
  expect_true(isSymmetric(V))
  expect_equal(unname(V[1, 1]), unname(vcov(fitA)[1, 1]), tolerance = 1e-10)
  expect_equal(unname(V[2, 2]), unname(vcov(fitB)[1, 1]), tolerance = 1e-10)

  wt <- wald_test(full = fitA, sub = fitB,
                  contrast = c("full.alpha" = 1, "sub.alpha" = -1))
  expect_s3_class(wt, "htest")
  expect_true(is.finite(wt$p.value))
})

test_that("saturated EM moments match lavaan's h1 estimator", {
  skip_if_not_installed("lavaan")
  X <- make_data(500L, 6L, miss = 0.2, seed = 6L)
  fit <- alpha_continuous(X)
  Mp <- lavaan:::lav_data_mi_patterns(X)
  lav <- lavaan:::lav_mvn_mi_h1_est_moments(
    y = X, mp = Mp, tol = 1e-12, max_iter = 10000L)
  expect_equal(unname(fit$moments$mu), unname(lav$Mu), tolerance = 1e-6)
  expect_equal(unname(fit$moments$Sigma), unname(lav$Sigma), tolerance = 1e-6)
})

test_that("sandwich SE matches coefficientalpha(varphi = 0)", {
  skip_if_not_installed("coefficientalpha")
  X <- make_data(500L, 6L, miss = 0.2, seed = 7L)
  fit <- alpha_continuous(X)
  invisible(utils::capture.output(
    ca <- coefficientalpha::alpha(X, varphi = 0, se = TRUE, test = FALSE,
                                  silent = TRUE)))
  expect_equal(unname(coef(fit)["alpha"]), as.numeric(ca$alpha),
               tolerance = 1e-4)
  expect_equal(sqrt(vcov(fit)[1, 1]), as.numeric(ca$se), tolerance = 5e-3 *
                 as.numeric(ca$se))
})
