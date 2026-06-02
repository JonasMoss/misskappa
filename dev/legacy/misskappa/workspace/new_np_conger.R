# Load required packages
library(Rcpp)
library(RcppArmadillo)

# Compile the C++ code
sourceCpp("workspace/new_np_conger.cpp")

# --- Example Usage ---

# 1. Create Sample Data
set.seed(123)
n_subjects <- 100
n_raters <- 5
n_categories <- 4

# True marginal probabilities for each rater
true_p_jc <- t(apply(matrix(runif(n_raters * n_categories), n_raters), 1, function(x) x / sum(x)))

# Generate complete data
complete_data <- matrix(0, nrow = n_subjects, ncol = n_raters)
for (j in 1:n_raters) {
  complete_data[, j] <- sample(1:n_categories, size = n_subjects, replace = TRUE, prob = true_p_jc[j, ])
}

# Introduce MCAR missingness (approx 20%)
ratings_mat_missing <- complete_data
ratings_mat_missing[sample(1:(n_subjects * n_raters), size = 0.2 * n_subjects * n_raters)] <- NA

head(ratings_mat_missing)

# 2. Define a Loss Matrix (Quadratic Weights)
loss_mat <- matrix(0, n_categories, n_categories)
for (c1 in 1:n_categories) {
  for (c2 in 1:n_categories) {
    loss_mat[c1, c2] <- (c1 - c2)^2
  }
}
print(loss_mat)


# 3. Run the Function
kappa_results <- conger_kappa_discrete(ratings_mat_missing, loss_mat)

# 4. Print Results
print("--- Estimates ---")
print(kappa_results$estimates)

print("--- Asymptotic Covariance Matrix ---")
print(kappa_results$vcov)

# 5. Get Confidence Interval for Kappa
kappa_est <- kappa_results$estimates["kappaC"]
kappa_se <- sqrt(kappa_results$vcov["kappaC", "kappaC"])
ci_lower <- kappa_est - 1.96 * kappa_se
ci_upper <- kappa_est + 1.96 * kappa_se

cat(sprintf("\nConger's Kappa: %.4f (95%% CI: %.4f, %.4f)\n", kappa_est, ci_lower, ci_upper))
