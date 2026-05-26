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
       cov = acov)
}

x <- irrCAC::cac.dist.g1g2[, 3:7]
kappa_dist(x) <- function(x, values = seq_len(ncol(x))) {
  n <- ncol(x)
  stopifnot(ncol(x) == length(values))

  x <- as.matrix(x)
  r <- sum(x[1, ])
  n = nrow(x)
  xt1 <- c(tcrossprod(values, x))
  xx = cbind(xt1, xt1^2, c(tcrossprod(values^2, x)))

  means <- colMeans(xx)
  theta <- tcrossprod(t(xx) - means) / n
  k <- r / (means[3] * r - means[1]^2)
  km <- k * (means[2] - means[1]^2)
  r_inv <- 1 / (r - 1)

  grad <- k * r_inv * c(2 * means[1] * (km / r - 1), 1, -km)
  est <- r_inv * (km - 1)
  var <- c(crossprod(grad, theta %*% grad))

  list(est = as.numeric(est), var = var / (n - 1))
}
