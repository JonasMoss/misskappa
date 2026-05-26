avar <- \(x, sigma, mu, type, fleiss, pi) {
  p <- colSums(!is.na(x)) / nrow(x)
  gamma <- gamma_est(x, sigma, type)
  mat <- cov_mat_kappa(p, mu, sigma, gamma, x, pi)
  avar_(mat, mu, sigma, fleiss)
}

avar_ <- function(mat, mu, sigma, fleiss) {
  r <- nrow(sigma)
  yy <- sum(diag(sigma))
  xx <- sum(sigma) - yy
  zz <- mean(mu^2) - mean(mu)^2

  vec <- if (fleiss) {
    mult <- 1 / ((r - 1) * (r * zz + yy)^2)
    c(r * zz + yy, r * zz - xx, -r * (xx + yy))
  } else {
    mult <- 1 / (r^2 * zz + (r - 1) * yy)^2
    c(r^2 * zz + (r - 1) * yy, -(r - 1) * xx, -r^2 * xx)
  }

  c(t(vec) %*% mat %*% vec) * mult^2
}

cov_mat_kappa <- function(p, mu, sigma, gamma, x, pi) {
  r <- length(p)

  # The covariances involving s.
  gamma_pi <- gamma * pi
  d <- get_diag_indices(r)
  ones <- rep(1, length(d))
  ones_minus_d <- ones - d
  cov_ss_ss <- t(ones_minus_d) %*% (gamma_pi) %*% ones_minus_d
  cov_ss_tr <- t(ones_minus_d) %*% (gamma_pi) %*% d
  cov_tr_tr <- t(d) %*% gamma_pi %*% d

  # The variance of mean(mu^2) - mean(mu)^2
  p_mat <- diag(1 / p)
  mu_middle <- (sigma + (p_mat - diag(r)) * diag(sigma))
  mu_vec <- (diag(r) - matrix(1 / r, r, r)) %*% mu
  cov_mu_mu <- c(4 / r^2 * t(mu_vec) %*% mu_middle %*% mu_vec)
  matrix(c(
    4 * cov_ss_ss, 2 * cov_ss_tr, 0,
    2 * cov_ss_tr, cov_tr_tr, 0,
    0, 0, cov_mu_mu
  ), nrow = 3, byrow = FALSE)
}

avar_bp <- \(x, type, c1, pi) {
  p <- colSums(!is.na(x)) / nrow(x)
  mu <- colMeans(x, na.rm = TRUE)
  sigma <- stats::cov(x, use = "pairwise.complete.obs")
  gamma <- gamma_est(x, sigma, type)
  var <- cov_mat_bp(p, mu, sigma, gamma, x, pi)
  avar_bp_(var, mu, sigma, c1)
}

cov_mat_bp <- function(p, mu, sigma, gamma, x, pi) {
  r <- length(p)

  # The variance of mean(mu^2) - mean(mu)^2
  p_mat <- diag(1 / p)
  mu_middle <- (sigma + (p_mat - diag(r)) * diag(sigma))
  mu_vec <- (diag(r) - matrix(1, r, r) / r) %*% mu
  var_a <- c(4 * t(mu_vec) %*% mu_middle %*% mu_vec)

  # The covariances involving Sigma.
  gamma_pi <- gamma * pi
  d <- get_diag_indices(r)
  ones <- rep(1, length(d))
  d_minus_ones <- d * (1 + 1 / r) - 2 / r * ones
  var_b <- c(t(d_minus_ones) %*% gamma_pi %*% d_minus_ones)
  var_a + var_b
}

avar_bp_ <- function(vars, mu, sigma, c1) {
  r <- nrow(sigma)
  4 / (r - 1)^2 / c1^2 * vars
}
