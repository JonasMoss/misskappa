# Complete-data generators, all calibrated to the SAME population covariance
# Sigma so the population alpha is held fixed across distributions.
#
# Nonnormality uses magmaan's NATIVE (C++) independent generator with the
# Pearson marginal family -- this matches Foldnes & Olsson (2016) and, unlike
# Vale-Maurelli or a nonnormal-common-factor-with-normal-errors construction,
# does NOT enjoy asymptotic robustness, so the normal-theory FIML SE genuinely
# breaks. magmaan is an experiment-only data-generation dependency (Suggests);
# the misskappa package keeps no magmaan runtime dependency.

.mag <- function(name) get(name, envir = asNamespace("magmaan"))

rmvn <- function(n, Sigma, mu = rep(0, ncol(Sigma))) {
  p <- ncol(Sigma)
  matrix(stats::rnorm(n * p), n, p) %*% chol(Sigma) + rep(mu, each = n)
}

# Elliptical multivariate t, scaled so Cov = Sigma (df > 2). df = 3 has NO
# finite fourth moment -- the torture case for the pairwise SE's finite-
# fourth-moment requirement.
rmvt_scaled <- function(n, Sigma, df = 3, mu = rep(0, ncol(Sigma))) {
  scale <- Sigma * (df - 2) / df
  z <- rmvn(n, scale)
  w <- stats::rchisq(n, df) / df
  z / sqrt(w) + rep(mu, each = n)
}

# epsilon-contamination mixture; base rescaled so the mixture covariance = Sigma.
rmv_contam <- function(n, Sigma, eps = 0.10, k = 5, mu = rep(0, ncol(Sigma))) {
  p <- ncol(Sigma)
  Sb <- Sigma / (1 + eps * (k^2 - 1))
  out <- rmvn(n, Sb)
  hit <- stats::runif(n) < eps
  if (any(hit)) out[hit, ] <- rmvn(sum(hit), Sb * k^2)
  out + rep(mu, each = n)
}

# magmaan independent generator (Pearson family). Calibration depends only on
# (Sigma, skew, excess kurtosis); calibrate once per cell, redraw per rep.
ig_calibrate <- function(Sigma, skew, excess_kurt) {
  p <- ncol(Sigma)
  .mag("sim_ig_calibrate")(
    Sigma,
    target_skewness = rep(skew, p),
    target_excess_kurtosis = rep(excess_kurt, p),
    generator_family = "pearson"
  )
}
ig_draw <- function(calibration, n, seed, mu = NULL) {
  X <- .mag("sim_ig_draw")(calibration, n = as.integer(n), reps = 1L,
                           seed_base = as.integer(seed))$draws[[1L]]
  if (!is.null(mu)) X <- X + rep(mu, each = nrow(X))
  X
}

# Per-distribution nonnormality targets (per-item skew / excess kurtosis).
ig_targets <- function(dist) switch(dist,
  ig_skew  = list(skew = 1.0, ek = 2.0),
  ig_heavy = list(skew = 2.0, ek = 7.0),
  stop("no IG targets for dist ", dist))

# Threshold a latent continuous battery into C ordered categories (1..C) using
# population (theoretical) thresholds, so category probabilities are ~equal.
discretize <- function(X, n_cat, Sigma) {
  sds <- sqrt(diag(Sigma))
  br <- stats::qnorm(seq(0, 1, length.out = n_cat + 1))
  out <- X
  for (j in seq_len(ncol(X))) {
    z <- X[, j] / sds[j]
    out[, j] <- findInterval(z, br[-c(1, length(br))]) + 1L
  }
  out
}

# Generate one complete dataset for a cell. `cal` is a pre-built IG calibration
# (ignored for non-IG dists). `seed` drives reproducibility.
generate_complete <- function(cell, Sigma, seed, cal = NULL) {
  n <- cell$n; dist <- cell$dist
  set.seed(seed)
  X <- switch(dist,
    normal   = rmvn(n, Sigma),
    t3       = rmvt_scaled(n, Sigma, df = 3),
    contam   = rmv_contam(n, Sigma),
    ig_skew  = ig_draw(cal, n, seed),
    ig_heavy = ig_draw(cal, n, seed),
    discrete = rmvn(n, Sigma),
    stop("unknown dist: ", dist)
  )
  if (identical(dist, "discrete")) X <- discretize(X, cell$n_cat, Sigma)
  colnames(X) <- paste0("y", seq_len(ncol(X)))
  X
}

# Population alpha for a cell. Closed form for continuous dists; Monte Carlo
# reference for discrete (alpha of the observed integer-score covariance).
cell_truth_alpha <- function(cell, Sigma, ref_n = 2e5, seed = 999L) {
  if (!identical(cell$dist, "discrete")) return(alpha_point_from_cov(Sigma))
  set.seed(seed)
  Xc <- discretize(rmvn(ref_n, Sigma), cell$n_cat, Sigma)
  alpha_point_from_cov(stats::cov(Xc))
}
