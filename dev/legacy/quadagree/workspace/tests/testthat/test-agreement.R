alleq <- \(x, y) isTRUE(all.equal(x, y))

test_that("population and x-based measures coincide", {
  x <- dat.zapf2016
  n <- nrow(x)
  values <- sort(unique(c(as.matrix(x))))
  sigma <- stats::cov(x) * (n - 1) / n
  mu <- colMeans(x)

  expect_equal(
    fleiss(x),
    fleiss_pop(mu, sigma)
  )
  expect_equal(
    conger(x),
    conger_pop(mu, sigma)
  )
  expect_equal(
    cohen(x),
    cohen_pop(mu, sigma)
  )
  expect_equal(
    bp(x, values, 1),
    bp_pop(mu, sigma, values, 1)
  )
})


test_that("Functions return errors for NA-saturated data", {
  expect_error(fleiss(dat.gwet2014))
  expect_error(cohen(dat.gwet2014))
  expect_error(conger(dat.gwet2014))
  expect_error(bp(dat.gwet2014))
})

test_that("BP agreement with irrCAC (long form)", {
  x <- dat.zapf2016
  est_irr <- (irrCAC::bp.coeff.raw(x, weights = "quadratic"))$est[4]
  est_quad <- bp(x, 1:5, 1)
  expect_equal(as.numeric(est_irr), est_quad)
})

test_that("BP agreement with irrCAC (wide form)", {
  x <- dat.fleiss1971
  est_irr <- irrCAC::bp.coeff.dist(x, weights = "quadratic")[2]
  est_quad <- bp_aggr(x)
  expect_equal(as.numeric(est_irr), est_quad)
})

test_that("BP agreement responds to kind argument", {
  x <- dat.fleiss1971
  est_1 <- bp_aggr(x, kind = 1)
  est_2 <- bp_aggr(x, kind = 2)
  expect_false(alleq(est_1, est_2))
})

test_that("BP agreement with irrCAC for se (wide form)", {
  x <- dat.fleiss1971
  n <- nrow(x)
  se_irr <- irrCAC::bp.coeff.dist(x, weights = "quadratic")[3]
  se_quad <- attr(bpci_aggr(x), "sd") / sqrt(n - 1)
  expect_equal(as.numeric(se_irr), se_quad)
})

test_that("Fleiss agreement with irrCAC for se (wide form)", {
  x <- dat.fleiss1971
  n <- nrow(x)
  se_irr <- irrCAC::fleiss.kappa.dist(x, weights = "quadratic")[3]
  se_quad <- attr(fleissci_aggr(x), "sd") / sqrt(n - 1)
  expect_equal(as.numeric(se_irr), se_quad)
})

test_that("Fleiss agreement with irrCAC for est (wide form)", {
  x <- dat.fleiss1971
  n <- nrow(x)
  est_irr <- c(irrCAC::fleiss.kappa.dist(x, weights = "quadratic")[2])
  est_quad <- fleiss_aggr(x)
  expect_equal(as.numeric(est_irr), est_quad)
})


test_that("Fleiss agreement with irrCAC (long form)", {
  x <- dat.zapf2016
  est_irr <- (irrCAC::fleiss.kappa.raw(x, weights = "quadratic"))$est[4]
  est_quad <- fleiss(x)
  expect_equal(as.numeric(est_irr), est_quad, tolerance = 1e-3)
})

test_that("Cohen agreement with irrCAC for se (long form)", {
  x <- dat.zapf2016
  n <- nrow(x)
  se_irr <- irrCAC::conger.kappa.raw(x, weights = "quadratic")$est[5]
  se_quad <- attr(cohenci(x), "sd") / sqrt(n - 3)
  expect_equal(as.numeric(se_irr), se_quad, tolerance = 1e-1)
})

test_that("Cohen agreement with irrCAC (long form)", {
  x <- dat.zapf2016
  est_irr <- irrCAC::conger.kappa.raw(x, weights = "quadratic")$est[4]
  est_quad <- attr(cohenci(x), "est")
  expect_equal(as.numeric(est_irr), est_quad, tolerance = 1e-7)
})
