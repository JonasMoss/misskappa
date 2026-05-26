bp_aggr_prepare <- \(x, values, kind) {
  y <- as.matrix(x)
  r <- sum(y[1, ])
  c1 <- bp_get_c1(values, kind)
  xtx <- c(tcrossprod(values^2, y))
  xt12 <- c(tcrossprod(values, y))^2
  list(
    xx = cbind(xtx, xt12),
    n = nrow(y),
    r = r,
    c1 = c1
  )
}

bp_aggr_fun <- \(calc) {
  means <- colMeans(calc$xx)
  theta <- tcrossprod(t(calc$xx) - means) / calc$n

  calcr_inv_1 <- 1 / (calc$r - 1)
  calcr_inv <- 1 / calc$r
  calcc1_inv <- 1 / calc$c1
  k <- theta[1, 1] - 2 * theta[1, 2] * calcr_inv + theta[2, 2] * calcr_inv^2
  disagreement <- 2 * calcr_inv_1 * (means[1] - calcr_inv * means[2])

  est <- unname(1 - disagreement * calcc1_inv)
  var <- calcc1_inv^2 * 4 * calcr_inv_1^2 * k
  list(est = est, var = max(var, 0))
}

fleiss_aggr_prepare <- \(x, values) {
  y <- as.matrix(x)
  r <- sum(y[1, ])
  xtx <- c(tcrossprod(values^2, y))
  xt1 <- c(tcrossprod(values, y))
  xt12 <- xt1^2

  list(
    xx = cbind(xt1, xt12, xtx),
    n = nrow(y),
    r = r
  )
}

fleiss_aggr_fun <- \(calc) {
  means <- colMeans(calc$xx)
  theta <- tcrossprod(t(calc$xx) - means) / calc$n
  k <- calc$r / (means[3] * calc$r - means[1]^2)
  km <- k * (means[2] - means[1]^2)
  calc_r_inv <- 1 / (calc$r - 1)

  grad <- k * calc_r_inv * c(2 * means[1] * (km / calc$r - 1), 1, -km)
  est <- calc_r_inv * (km - 1)
  var <- c(crossprod(grad, theta %*% grad))

  list(est = est, var = var)
}
