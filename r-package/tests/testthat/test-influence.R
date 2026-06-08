test_that("influence() returns the per-subject IF for categorical raw fits", {
  set.seed(1)
  n <- 200L
  R <- 4L
  C <- 5L
  x <- matrix(sample.int(C, n * R, replace = TRUE), n, R)
  x[matrix(stats::runif(n * R) < 0.2, n, R)] <- NA
  storage.mode(x) <- "integer"

  for (m in c("available", "ipw", "gwet")) {
    fit <- estimate_kappa_raw(x, method =m, weight = "identity")
    psi <- fit$psi
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
      fit <- estimate_kappa_raw(x, method =m, weight = w)
      psi <- fit$psi
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

test_that("influence() returns per-subject IFs for FIML fits", {
  x <- matrix(
    c(
      1, 1, 1,
      2, 2, NA,
      3, 3, 3,
      1, 1, 2,
      NA, 3, 2,
      3, 2, 3,
      1, 2, 1,
      2, 2, 3,
      3, 3, 3,
      1, 1, 1,
      2, 1, 2,
      3, 3, 2
    ),
    nrow = 12, byrow = TRUE
  )
  storage.mode(x) <- "integer"
  fit <- estimate_kappa_raw(x, method ="fiml", weight = "quadratic")
  psi <- fit$psi
  expect_equal(dim(psi), c(12L, 3L))
  expect_equal(unname(crossprod(psi) / 144), unname(vcov(fit)),
               tolerance = 1e-10)
})

test_that("influence() returns per-subject IFs for counts, continuous, and g-wise fits", {
  counts <- matrix(c(5, 5, 0, 8, 2, 0, 0, 3, 7, 4, 1, 5), nrow = 4, byrow = TRUE)
  storage.mode(counts) <- "integer"
  fit_counts <- kappa_counts(counts, estimator = "pairwise")
  psi_counts <- fit_counts$psi
  expect_equal(dim(psi_counts), c(4L, 2L))
  expect_equal(unname(crossprod(psi_counts) / 16), unname(vcov(fit_counts)),
               tolerance = 1e-10)

  xc <- matrix(stats::rnorm(40), nrow = 10)
  fit_cont <- kappa_continuous(xc, method = "available", weight = "quadratic")
  psi_cont <- fit_cont$psi
  expect_equal(dim(psi_cont), c(10L, 2L))
  expect_equal(unname(crossprod(psi_cont) / 100), unname(vcov(fit_cont)),
               tolerance = 1e-10)

  xg <- matrix(c(1, 1, 2, 1, 1, 2, 2, 3, 3), nrow = 3, byrow = TRUE)
  fit_gwise <- kappa_gwise(xg, distance = "nominal")
  psi_gwise <- fit_gwise$psi
  expect_equal(dim(psi_gwise), c(3L, 2L))
})

test_that("influence() returns per-subject IFs for FIML counts fits", {
  counts <- matrix(
    c(
      3, 1, 0,
      0, 2, 1,
      2, 0, 1,
      1, 2, 0,
      0, 0, 3,
      2, 1, 0,
      3, 0, 0,
      0, 1, 2
    ),
    nrow = 8, byrow = TRUE
  )
  storage.mode(counts) <- "integer"
  fit <- kappa_counts(counts, estimator = "cat_fiml", r_total = 4)
  psi <- fit$psi
  expect_equal(dim(psi), c(8L, 2L))
  expect_equal(unname(crossprod(psi) / 64), unname(vcov(fit)),
               tolerance = 1e-10)
})

test_that("influence() returns per-subject IFs for quadratic fits", {
  x <- matrix(
    c(0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1),
    nrow = 10, byrow = TRUE
  )
  storage.mode(x) <- "integer"

  xc <- matrix(stats::rnorm(40), nrow = 10)
  values <- c(1, 2, 3, 4, 5)[seq_len(min(ncol(xc), 5))]
  fit_quad <- kappa_quadratic(xc, values = values)
  psi_quad <- fit_quad$psi
  expect_equal(dim(psi_quad), c(10L, 3L))
  expect_equal(unname(crossprod(psi_quad) / 100), unname(vcov(fit_quad)),
               tolerance = 1e-10)
})

test_that("joint_vcov() assembles a block matrix and recovers per-fit vcov on the diagonal blocks", {
  set.seed(3)
  n <- 250L
  R <- 4L
  C <- 5L
  x <- matrix(sample.int(C, n * R, replace = TRUE), n, R)
  x[matrix(stats::runif(n * R) < 0.3, n, R)] <- NA
  storage.mode(x) <- "integer"

  ac  <- estimate_kappa_raw(x, method ="available")
  ipw <- estimate_kappa_raw(x, method ="ipw")
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

  fit1 <- estimate_kappa_raw(x1, method = "available")
  fit2 <- estimate_kappa_raw(x2, method = "available")
  expect_error(joint_vcov(fit1, fit2), "same number of subjects")

  expect_error(joint_vcov(fit1), "at least two")
})

test_that("joint_vcov() supports FIML fits", {
  x <- matrix(
    c(
      1, 1, 1,
      2, 2, NA,
      3, 3, 3,
      1, 1, 2,
      NA, 3, 2,
      3, 2, 3,
      1, 2, 1,
      2, 2, 3,
      3, 3, 3,
      1, 1, 1,
      2, 1, 2,
      3, 3, 2
    ),
    nrow = 12, byrow = TRUE
  )
  storage.mode(x) <- "integer"
  ac <- estimate_kappa_raw(x, method ="available", weight = "quadratic")
  fiml <- estimate_kappa_raw(x, method ="fiml", weight = "quadratic")
  V <- joint_vcov(ac = ac, fiml = fiml)

  expect_equal(dim(V), c(6L, 6L))
  expect_equal(unname(V[1:3, 1:3]), unname(vcov(ac)), tolerance = 1e-10)
  expect_equal(unname(V[4:6, 4:6]), unname(vcov(fiml)), tolerance = 1e-10)
})

test_that("Hausman contrast kappa_AC - kappa_IPW agrees between joint_vcov and direct computation", {
  set.seed(5)
  n <- 400L
  R <- 4L
  C <- 5L
  x <- matrix(sample.int(C, n * R, replace = TRUE), n, R)
  x[matrix(stats::runif(n * R) < 0.3, n, R)] <- NA
  storage.mode(x) <- "integer"

  ac  <- estimate_kappa_raw(x, method ="available")
  ipw <- estimate_kappa_raw(x, method ="ipw")
  V   <- joint_vcov(ac = ac, ipw = ipw)

  c_vec <- c(1, 0, 0, -1, 0, 0)
  v_contrast <- as.numeric(crossprod(c_vec, V %*% c_vec))

  # Directly: Var(kappa_AC.Conger - kappa_IPW.Conger) via the n x 6
  # stacked IF.
  psi <- cbind(ac$psi, ipw$psi)
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

  ac <- estimate_kappa_raw(x, method ="available")
  single <- wald_test(ac, contrast = "Conger", value = 0)
  expect_s3_class(single, "htest")
  expect_named(single$statistic, "X-squared")
  expect_equal(unname(single$estimate), unname(coef(ac)["Conger"]))
  expect_equal(unname(single$parameter), 1)
  expect_true(is.finite(single$p.value))

  ipw <- estimate_kappa_raw(x, method ="ipw")
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
  fit <- estimate_kappa_raw(x, method ="available")

  expect_error(wald_test(fit, contrast = "Nope"), "Unknown coefficient")
  expect_error(wald_test(fit, contrast = c(1, 0)), "length")

  expect_error(wald_test(fit, contrast = c("fit.Conger" = 1, "fit.Nope" = -1)),
               "Unknown coefficient")
})
