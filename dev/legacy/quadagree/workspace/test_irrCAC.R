library(testthat)
set.seed(2024)

N_SUBJECTS <- 100
N_RATERS <- 5
CATEGORIES <- 1:5

x_raw <- as.matrix(
  replicate(N_RATERS, sample(CATEGORIES, N_SUBJECTS, replace = TRUE))
)

x_aggr <- t(apply(x_raw, 1, function(row) {
  table(factor(row, levels = CATEGORIES))
}))

mod_quadagree <- kappa_raw(x = x_raw)

mod_conger_irrcac <- irrCAC::conger.kappa.raw(x_raw, weights = "quadratic")
mod_fleiss_irrcac <- irrCAC::fleiss.kappa.raw(x_raw, weights = "quadratic")
mod_bp_irrcac <- irrCAC::bp.coeff.raw(x_raw, weights = "quadratic")

mod_quadagree_aggr <- kappa_aggr(x = x_aggr, values = CATEGORIES, r = N_RATERS)
mod_fleiss_irrcac_aggr <- irrCAC::fleiss.kappa.dist(x_aggr, weights = "quadratic")
mod_fleiss_bp_aggr <- irrCAC::bp.coeff.dist(x_aggr, weights = "quadratic")

test_that("kappa_raw estimates match irrCAC for complete data", {

  est_quadagree_conger <- mod_quadagree$conger$estimate
  est_irrcac_conger <- mod_conger_irrcac$est$coeff.val
  expect_equal(est_quadagree_conger, est_irrcac_conger, tolerance = 1e-4,
               label = "Conger's kappa estimate should match irrCAC")

  est_quadagree_fleiss <- mod_quadagree$fleiss$estimate
  est_irrcac_fleiss <- mod_fleiss_irrcac$est$coeff.val
  expect_equal(est_quadagree_fleiss, est_irrcac_fleiss, tolerance = 1e-4,
               label = "Fleiss' kappa estimate should match irrCAC")

  est_quadagree_bp <- mod_quadagree$bp$estimate
  est_irrcac_bp <- mod_bp_irrcac$est$coeff.val
  expect_equal(est_quadagree_bp, est_irrcac_bp, tolerance = 1e-4,
               label = "bp' kappa estimate should match irrCAC")
})


test_that("kappa_raw standard errors match irrCAC (with n-1 correction)", {

  n_eff <- mod_quadagree$fleiss$n_eff

  se_quadagree_conger_raw <- mod_quadagree$conger$std.err
  se_quadagree_conger_adj <- se_quadagree_conger_raw * sqrt(n_eff / (n_eff - 1))
  se_irrcac_conger <- mod_conger_irrcac$est$coeff.se
  expect_equal(se_quadagree_conger_adj, se_irrcac_conger, tolerance = 1e-3,
               label = "Conger's standard error should match irrCAC after n-1 adjustment")

  se_quadagree_fleiss_raw <- mod_quadagree$fleiss$std.err
  se_quadagree_fleiss_adj <- se_quadagree_fleiss_raw * sqrt(n_eff / (n_eff - 1))
  se_irrcac_fleiss <- mod_fleiss_irrcac$est$coeff.se
  expect_equal(se_quadagree_fleiss_adj, se_irrcac_fleiss, tolerance = 1e-4,
               label = "fleiss' standard error should match irrCAC after n-1 adjustment")

  se_quadagree_bp_raw <- mod_quadagree$bp$std.err
  se_quadagree_bp_adj <- se_quadagree_bp_raw * sqrt(n_eff / (n_eff - 1))
  se_irrcac_bp <- mod_bp_irrcac$est$coeff.se
  expect_equal(se_quadagree_bp_adj, se_irrcac_bp, tolerance = 1e-4,
               label = "bp' standard error should match irrCAC after n-1 adjustment")


})


test_that("kappa_aggr estimates match irrCAC for complete data", {

  est_quadagree_fleiss <- mod_quadagree_aggr$fleiss$estimate
  est_irrcac_fleiss <- mod_fleiss_irrcac_aggr$coeff
  expect_equal(est_quadagree_fleiss, est_irrcac_fleiss, tolerance = 1e-6,
               label = "Aggregated Fleiss' kappa estimate should match irrCAC")

  est_quadagree_bp_raw <- mod_quadagree_aggr$bp$estimate
  est_quadagree_bp_adj <- est_quadagree_bp_raw
  est_irrcac_bp <- mod_fleiss_bp_aggr$coeff
  expect_equal(est_quadagree_bp_adj, est_irrcac_bp, tolerance = 1e-6,
               label = "Aggregated BP estimate should match irrCAC")

})


test_that("kappa_aggr standard errors match irrCAC (with n-1 correction)", {

  # --- Renormalize your standard error ---
  n_eff <- mod_quadagree_aggr$fleiss$n_eff

  se_quadagree_fleiss_raw <- mod_quadagree_aggr$fleiss$std.err
  se_quadagree_fleiss_adj <- se_quadagree_fleiss_raw * sqrt(n_eff / (n_eff - 1))
  se_irrcac_fleiss <- mod_fleiss_irrcac_aggr$stderr
  expect_equal(se_quadagree_fleiss_adj, se_irrcac_fleiss, tolerance = 1e-6,
               label = "Aggregated Fleiss' standard error should match irrCAC after n-1 adjustment")


  se_quadagree_bp_raw <- mod_quadagree_aggr$bp$std.err
  se_quadagree_bp_adj <- se_quadagree_bp_raw * sqrt(n_eff / (n_eff - 1))
  se_irrcac_bp <- mod_fleiss_bp_aggr$stderr
  expect_equal(se_quadagree_bp_adj, se_irrcac_bp, tolerance = 1e-6,
               label = "Aggregated BP standard error should match irrCAC after n-1 adjustment")


})
