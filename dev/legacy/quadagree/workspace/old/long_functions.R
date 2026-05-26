fleiss_prepare <- \(x, type) {
  x <- as.matrix(x)
  pi <- pi_mat_empirical(x)
  list(xx = as.matrix(x), type = type, n = nrow(x), pi = pi)
}

fleiss_fun <- \(calc) {
  n <- calc$n
  sigma <- stats::cov(calc$xx, use = "pairwise.complete.obs") * (n - 1) / n
  if (any(is.na(sigma))) stop("The data does not contain sufficient non-NAs.")
  mu <- colMeans(calc$xx, na.rm = TRUE)
  est <- fleiss_pop(mu, sigma)
  var <- avar(calc$xx, sigma, mu, calc$type, TRUE, calc$pi)
  list(est = est, var = var)
}

conger_prepare <- \(x, type) {
  x <- as.matrix(x)
  pi <- pi_mat_empirical(x)
  list(xx = x, type = type, n = nrow(x), pi = pi)
}

conger_fun <- \(calc) {
  n <- calc$n
  sigma <- stats::cov(calc$xx, use = "pairwise.complete.obs") * (n - 1) / n
  if (any(is.na(sigma))) stop("The data does not contain sufficient non-NAs.")
  mu <- colMeans(calc$xx, na.rm = TRUE)
  est <- conger_pop(mu, sigma)
  var <- avar(calc$xx, sigma, mu, calc$type, FALSE, calc$pi)
  list(est = est, var = var)
}

bp_prepare <- \(x, values, kind, type) {
  x <- as.matrix(x)
  pi <- pi_mat_empirical(x)
  if (is.null(values)) values <- stats::na.omit(unique(c(x)))
  c1 <- bp_get_c1(values, kind)
  list(xx = x, c1 = c1, type = type, n = nrow(x), pi = pi)
}

bp_fun <- \(calc) {
  est <- bp_calc(calc$xx, calc$c1)
  var <- avar_bp(calc$xx, calc$type, calc$c1, calc$pi)
  list(est = est, var = var)
}
