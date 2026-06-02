test_that("influence() returns the per-subject IF for categorical raw fits", {
  set.seed(1)
  n <- 200L
  R <- 4L
  C <- 5L
  x <- matrix(sample.int(C, n * R, replace = TRUE), n, R)
  x[matrix(stats::runif(n * R) < 0.2, n, R)] <- NA
  storage.mode(x) <- "integer"

  for (m in c("available", "ipw", "gwet")) {
    fit <- kappa(x, method = m, weight = "identity")
    psi <- stats::influence(fit)
    expect_true(is.matrix(psi))
    expect_equal(dim(psi), c(n, 3L))
    expect_equal(colnames(psi), c("Conger", "Fleiss", "Brennan-Prediger"))
  }
})

test_that("vcov = (1 / n^2) crossprod(psi) to floating-point noise", {
  set.seed(2)
  n <- 300L
  R <- 5L
  C <- 4L
  x <- matrix(sample.int(C, n * R, replace = TRUE), n, R)
  x[matrix(stats::runif(n * R) < 0.25, n, R)] <- NA
  storage.mode(x) <- "integer"

  for (m in c("available", "ipw", "gwet")) {
    for (w in c("identity", "linear", "quadratic")) {
      fit <- kappa(x, method = m, weight = w)
      psi <- stats::influence(fit)
      v_reconstructed <- crossprod(psi) / (n * n)
      v_orig <- vcov(fit)
      expect_equal(
        v_reconstructed,
        v_orig,
        tolerance = 1e-10,
        ignore_attr = TRUE,
        info = sprintf("method=%s, weight=%s", m, w)
      )
    }
  }
})

test_that("influence() returns NULL for estimators that do not expose IFs", {
  x <- matrix(
    c(0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1),
    nrow = 10, byrow = TRUE
  )
  storage.mode(x) <- "integer"

  # FIML uses EM + Louis info; no IF exposed in this version.
  fit_fiml <- kappa(x, method = "fiml", weight = "identity")
  expect_null(stats::influence(fit_fiml))

  # Counts, continuous, and quadratic also do not expose IFs (yet).
  counts <- matrix(c(5, 5, 0, 8, 2, 0, 0, 3, 7, 4, 1, 5), nrow = 4, byrow = TRUE)
  storage.mode(counts) <- "integer"
  fit_counts <- kappa_counts(counts, method = "available")
  expect_null(stats::influence(fit_counts))

  xc <- matrix(stats::rnorm(40), nrow = 10)
  fit_cont <- kappa_continuous(xc, method = "available", weight = "quadratic")
  expect_null(stats::influence(fit_cont))

  fit_quad <- kappa_quadratic(xc, values = c(1, 2, 3, 4, 5)[seq_len(min(ncol(xc), 5))])
  expect_null(stats::influence(fit_quad))
})

test_that("joint_vcov() assembles a block matrix and recovers per-fit vcov on the diagonal blocks", {
  set.seed(3)
  n <- 250L
  R <- 4L
  C <- 5L
  x <- matrix(sample.int(C, n * R, replace = TRUE), n, R)
  x[matrix(stats::runif(n * R) < 0.3, n, R)] <- NA
  storage.mode(x) <- "integer"

  ac  <- kappa(x, method = "available")
  ipw <- kappa(x, method = "ipw")
  V   <- joint_vcov(ac = ac, ipw = ipw)

  expect_equal(dim(V), c(6L, 6L))
  expect_true(isSymmetric(V))

  block_ac  <- V[1:3, 1:3]
  block_ipw <- V[4:6, 4:6]
  expect_equal(unname(block_ac),  unname(vcov(ac)),  tolerance = 1e-10)
  expect_equal(unname(block_ipw), unname(vcov(ipw)), tolerance = 1e-10)
  expect_equal(
    rownames(V),
    c("ac.Conger", "ac.Fleiss", "ac.Brennan-Prediger",
      "ipw.Conger", "ipw.Fleiss", "ipw.Brennan-Prediger")
  )
})

test_that("joint_vcov() errors on mismatched n or non-IF fits", {
  set.seed(4)
  x1 <- matrix(sample.int(4L, 200, replace = TRUE), 50, 4)
  x2 <- matrix(sample.int(4L, 240, replace = TRUE), 60, 4)
  storage.mode(x1) <- "integer"
  storage.mode(x2) <- "integer"

  fit1 <- kappa(x1, method = "available")
  fit2 <- kappa(x2, method = "available")
  expect_error(joint_vcov(fit1, fit2), "same number of subjects")

  fit_fiml <- kappa(x1, method = "fiml")
  expect_error(joint_vcov(fit1, fit_fiml), "do not expose influence")

  expect_error(joint_vcov(fit1), "at least two")
})

test_that("Hausman contrast kappa_AC - kappa_IPW agrees between joint_vcov and direct computation", {
  set.seed(5)
  n <- 400L
  R <- 4L
  C <- 5L
  x <- matrix(sample.int(C, n * R, replace = TRUE), n, R)
  x[matrix(stats::runif(n * R) < 0.3, n, R)] <- NA
  storage.mode(x) <- "integer"

  ac  <- kappa(x, method = "available")
  ipw <- kappa(x, method = "ipw")
  V   <- joint_vcov(ac = ac, ipw = ipw)

  c_vec <- c(1, 0, 0, -1, 0, 0)
  v_contrast <- as.numeric(crossprod(c_vec, V %*% c_vec))

  # Directly: Var(kappa_AC.Conger - kappa_IPW.Conger) via the n x 6
  # stacked IF.
  psi <- cbind(stats::influence(ac), stats::influence(ipw))
  diff_psi <- psi[, 1] - psi[, 4]   # ac.Conger - ipw.Conger
  v_direct <- sum(diff_psi^2) / (n * n)

  expect_equal(v_contrast, v_direct, tolerance = 1e-10)
})

test_that("wald_test() tests single-fit and joint linear contrasts", {
  set.seed(6)
  n <- 350L
  R <- 4L
  C <- 5L
  x <- matrix(sample.int(C, n * R, replace = TRUE), n, R)
  x[matrix(stats::runif(n * R) < 0.25, n, R)] <- NA
  storage.mode(x) <- "integer"

  ac <- kappa(x, method = "available")
  single <- wald_test(ac, contrast = "Conger", value = 0)
  expect_s3_class(single, "htest")
  expect_named(single$statistic, "X-squared")
  expect_equal(unname(single$estimate), unname(coef(ac)["Conger"]))
  expect_equal(unname(single$parameter), 1)
  expect_true(is.finite(single$p.value))

  ipw <- kappa(x, method = "ipw")
  joint <- wald_test(
    ac = ac,
    ipw = ipw,
    contrast = c("ac.Conger" = 1, "ipw.Conger" = -1)
  )
  V <- joint_vcov(ac = ac, ipw = ipw)
  c_vec <- c(1, 0, 0, -1, 0, 0)
  delta <- coef(ac)["Conger"] - coef(ipw)["Conger"]
  stat <- as.numeric(delta * delta / crossprod(c_vec, V %*% c_vec))

  expect_equal(unname(joint$statistic), stat, tolerance = 1e-10)
  expect_equal(unname(joint$parameter), 1)
  expect_true(is.finite(joint$p.value))
})

test_that("wald_test() validates contrasts and influence-function support", {
  x <- matrix(sample.int(4L, 200, replace = TRUE), 50, 4)
  storage.mode(x) <- "integer"
  fit <- kappa(x, method = "available")

  expect_error(wald_test(fit, contrast = "Nope"), "Unknown coefficient")
  expect_error(wald_test(fit, contrast = c(1, 0)), "length")

  fit_fiml <- kappa(x, method = "fiml")
  expect_error(wald_test(fit, fit_fiml, contrast = c("fit1.Conger" = 1)),
               "do not expose influence")
})
