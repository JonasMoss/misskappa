

raw_data <- as.matrix(dat.zapf2016)
num_categories_raw <- 5
values_raw <- c(1, 2, 3, 4, 5)

counts_data <- as.matrix(dat.fleiss1971)
num_raters_counts <- 6
values_counts <- 1:5

em_opts <- list(
  tol = 1e-8,
  max_iter = 1000,
  prune_tol = 1e-9,
  start_alpha = 0.1
)

results_ml_raw <- kappa_ml_raw_cpp(
  x_r = raw_data,
  c = num_categories_raw,
  weight_type = "identity",
  values = values_raw,
  em_options = em_opts
)

n <- nrow(raw_data)
as.numeric(c(
  irrCAC::conger.kappa.raw(raw_data)$est["coeff.se"],
  irrCAC::fleiss.kappa.raw(raw_data)$est["coeff.se"],
  irrCAC::bp.coeff.raw(raw_data)$est["coeff.se"]
))^2 * (n - 1) / n

diag(results_ml_raw$vcov)



results_ml_counts <- kappa_ml_counts_cpp(
  x_r = counts_data,
  r = num_raters_counts,
  weight_type = "unweighted", # Using unweighted for this test
  values = values_counts,
  em_options = em_opts
)

n <- nrow(counts_data)
as.numeric(c(
  irrCAC::fleiss.kappa.dist(counts_data)["stderr"],
  irrCAC::bp.coeff.dist(counts_data)["stderr"]
))^2 * (n - 1) / n

diag(results_ml_counts$vcov)

