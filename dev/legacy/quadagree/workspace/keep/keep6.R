kappa <- function(x) {
  x <- as.matrix(x)
  n <- nrow(x)
  r <- ncol(x)

  mu <- colMeans(x)
  sigma <- cov(x) * (n - 1) / n

  tr <- sum(diag(sigma))
  mu_bar <- sum(mu) / r
  d_mu <- sum((mu - mu_bar)^2)

  n_c <- sum(sigma) - tr
  n_f <- n_c - d_mu

  d_c <- (r - 1) * tr + r * d_mu
  d_f <- d_c - d_mu

  kappa_c <- n_c / d_c
  kappa_f <- n_f / d_f

  y <- sweep(x, 2, mu, "-")
  z <- cbind(y %*% (mu - mu_bar), rowSums(y)^2, rowSums(y^2))
  xi <- cov(z) * (n-1) / n

  cohen <- c(2 * r * kappa_c, -1, 1 + (r - 1) * kappa_c) / d_c
  fleiss <- c(2 * (1 + (r - 1) * kappa_f), -1, 1 + (r - 1) * kappa_f) / d_f
  v <- cbind(cohen, fleiss)
  acov <- t(v) %*% xi %*% v

  list(cohen = kappa_c,
       fleiss = kappa_f,
       cov = acov / (n-1))
}


get_xi_naive <- function(x, mu) {
  n <- nrow(x)
  r <- ncol(x)
  M <- !is.na(x)

  Y_hat <- sweep(x, 2, mu, "-")

  p1 <- colMeans(M)
  p2 <- crossprod(M) / n

  Z_star_matrix <- matrix(NA, nrow = n, ncol = 3)

  mu_bar <- mean(mu)
  mu_diff <- mu - mu_bar

  inv_p1 <- 1 / p1
  inv_p2 <- 1 / p2

  for (i in 1:n) {
    Mi <- M[i, ]
    if (sum(Mi) == 0) next

    Yi <- Y_hat[i, ]
    obs_raters <- which(Mi)

    z1 <- sum(mu_diff[Mi] * Yi[Mi] * inv_p1[Mi])
    z3 <- sum(Yi[Mi]^2 * inv_p1[Mi])

    z2_diag <- z3
    z2_offdiag <- 0
    if (length(obs_raters) > 1) {
      pairs <- combn(obs_raters, 2)
      for(k in 1:ncol(pairs)) {
        j1 <- pairs[1, k]; j2 <- pairs[2, k]
        z2_offdiag <- z2_offdiag + 2 * Yi[j1] * Yi[j2] * inv_p2[j1, j2]
      }
    }
    z2 <- z2_diag + z2_offdiag
    Z_star_matrix[i, ] <- c(z1, z2, z3)
  }

  complete_rows_idx <- which(!is.na(Z_star_matrix[,1]))
  Z_star_eff <- Z_star_matrix[complete_rows_idx, ]
  n_eff <- nrow(Z_star_eff)

  xi <- cov(Z_star_eff) * (n_eff - 1) / n_eff

  return(list(xi = xi, n_eff = n_eff))
}

kappa_naive_Z <- function(x) {
  x <- as.matrix(x)
  n <- nrow(x)
  r <- ncol(x)

  mu <- colMeans(x, na.rm = TRUE)
  sigma <- cov(x, use = "pairwise.complete.obs") * (n - 1) / n

  tr <- sum(diag(sigma))
  mu_bar <- sum(mu) / r
  d_mu <- sum((mu - mu_bar)^2)

  n_c <- sum(sigma) - tr
  n_f <- n_c - d_mu

  d_c <- (r - 1) * tr + r * d_mu
  d_f <- d_c - d_mu

  kappa_c <- n_c / d_c
  kappa_f <- n_f / d_f

  xi_res <- get_xi_naive(x, mu)
  xi <- xi_res$xi
  n_eff <- xi_res$n_eff

  cohen <- c(2 * r * kappa_c, -1, 1 + (r - 1) * kappa_c) / d_c
  fleiss <- c(2 * (1 + (r - 1) * kappa_f), -1, 1 + (r - 1) * kappa_f) / d_f
  v <- cbind(fleiss, cohen)

  scaled_acov <- t(v) %*% xi %*% v

  list(cohen = kappa_c,
       fleiss = kappa_f,
       scaled_acov = scaled_acov,
       n_eff = n_eff)
}

# --- Verification Run ---
set.seed(135)
n_sim <- 50000; r_sim <- 4; missing_prob <- 0.2
true_mu <- c(1, 1.5, 1.2, 2); true_sigma <- matrix(0.7, r_sim, r_sim); diag(true_sigma) <- 1
L <- chol(true_sigma)
complete_data <- matrix(rnorm(n_sim * r_sim), n_sim, r_sim) %*% L
complete_data <- sweep(complete_data, 2, true_mu, "+")
x_missing <- complete_data
x_missing[sample(1:(n_sim * r_sim), size = n_sim * r_sim * missing_prob)] <- NA

cat("\n--- Naive Z Method (minimal modification) ---\n")
print(kappa_naive_Z(x_missing)$scaled_acov)

cat("\n--- Correct Formal Method (for comparison) ---\n")
print(kappa_acov_formal_CORRECTED(x_missing)$scaled_acov)
