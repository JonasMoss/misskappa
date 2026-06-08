# Measurement models -> population item covariance Sigma. One-factor
# congeneric family (factor variance 1): Sigma = lambda lambda' + diag(psi).
# alpha is estimated qua alpha, so only Sigma matters for the estimand; the
# "wild" model deliberately leaves the congeneric comfort zone with a junk
# (near-zero loading) item, a dominant-variance item, and a reverse-keyed
# (negative loading) item.

make_model <- function(model, p) {
  switch(model,
    parallel = list(lambda = rep(0.70, p), psi = rep(0.51, p)),
    tau      = list(lambda = rep(0.70, p), psi = seq(0.30, 0.90, length.out = p)),
    congeneric_mild = list(lambda = seq(0.50, 0.85, length.out = p), psi = rep(0.50, p)),
    congeneric_wild = {
      lam <- seq(0.40, 0.90, length.out = p)
      lam[1] <- 0.05                       # junk item: near-zero loading
      if (p >= 3) lam[p] <- -0.60          # reverse-keyed item
      psi <- rep(0.50, p); psi[min(2L, p)] <- 4.0  # dominant-variance item
      list(lambda = lam, psi = psi)
    },
    stop("unknown model: ", model)
  )
}

sigma_from_model <- function(m) outer(m$lambda, m$lambda) + diag(m$psi)

# Scale residual variances by a constant to hit a target population alpha
# (alpha is monotone decreasing in the residual scale). NA target keeps the
# model's natural alpha.
scale_to_alpha <- function(m, target) {
  if (is.na(target)) return(m)
  f <- function(cscale) {
    mm <- m; mm$psi <- m$psi * cscale
    alpha_point_from_cov(sigma_from_model(mm)) - target
  }
  cscale <- tryCatch(stats::uniroot(f, c(1e-4, 1e4))$root, error = function(e) NA_real_)
  if (is.finite(cscale)) m$psi <- m$psi * cscale
  m
}

# Build the population Sigma for a (model, p, alpha_target) cell.
cell_sigma <- function(model, p, alpha_target = NA_real_) {
  m <- scale_to_alpha(make_model(model, p), alpha_target)
  sigma_from_model(m)
}
