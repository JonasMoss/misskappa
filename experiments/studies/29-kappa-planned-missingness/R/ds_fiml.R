# Dawid-Skene FIML prototype (experiment-local R implementation)
#
# Latent-class conditional-independence model: subject i carries a latent
# class k in 1..K (K = number of observed categories); rater j emits category
# c with probability theta_j(c | k), independently across raters given k.
# Conditional independence makes the observed-data likelihood factor over the
# observed raters, so MAR and planned designs need no C^R marginalisation --
# a missing rating simply drops out of the product. With ~K(K-1) + R*K*(C-1)
# parameters this sits between the saturated cat_fiml (C^R - 1, intrinsically
# biased at small n; study 29) and the exchangeable count FIML (rater
# exchangeability, Fleiss_exch only): it keeps rater-specific confusion
# matrices, so it targets the fixed-rater Conger / Fleiss estimands.
#
# Inference: delta method through the model-implied kappa map. Per-subject
# observed-data scores are analytic (Fisher's identity); the observed
# information comes from central differences of the analytic score; both are
# taken in redundant softmax coordinates, whose gauge directions (one per
# simplex) are exact null directions of the information and are removed by a
# rank-truncated spectral pseudo-inverse, mirroring the saturated cat_fiml
# info_rcond treatment. Scores and the kappa gradient are orthogonal to the
# gauge directions, so the pseudo-inverse delta method is exact in the
# quotient. Two flavours are kept: "ml" (bread only, correct under the model)
# and "robust" (bread-OPG-bread sandwich, correct for the pseudo-true value
# under misspecification). vcov() surfaces the robust flavour.
#
# Latent-class label switching is irrelevant for kappa because the implied
# rating distribution is invariant to class permutation.
#
# Sourced by the study-29 runner and by workbench probe 37-ds-fiml-se; keep
# this file free of dependencies on the surrounding script.

ds_em_fit <- function(Y, C, n_starts = 3L, max_iter = 2000L, tol = 1e-8) {
  n <- nrow(Y)
  R <- ncol(Y)
  K <- C
  cat_counts <- vapply(seq_len(C), function(c) rowSums(Y == c, na.rm = TRUE),
                       numeric(n))

  run_one <- function(W) {
    ll <- -Inf
    pi_k <- rep(1 / K, K)
    theta <- NULL
    for (iter in seq_len(max_iter)) {
      pi_k <- pmax(colMeans(W), 1e-12)
      pi_k <- pi_k / sum(pi_k)
      theta <- lapply(seq_len(R), function(j) {
        oj <- !is.na(Y[, j])
        Wo <- W[oj, , drop = FALSE]
        xo <- Y[oj, j]
        th <- matrix(0, K, C)
        for (c in seq_len(C)) {
          rows <- xo == c
          if (any(rows)) th[, c] <- colSums(Wo[rows, , drop = FALSE])
        }
        den <- colSums(Wo)
        th <- sweep(th, 1L, pmax(den, 1e-12), "/")
        th[den < 1e-12, ] <- 1 / C
        th <- pmax(th, 1e-12)
        th / rowSums(th)
      })
      logf <- matrix(log(pi_k), n, K, byrow = TRUE)
      for (j in seq_len(R)) {
        oj <- !is.na(Y[, j])
        logf[oj, ] <- logf[oj, ] + t(log(theta[[j]]))[Y[oj, j], , drop = FALSE]
      }
      mx <- apply(logf, 1L, max)
      lse <- mx + log(rowSums(exp(logf - mx)))
      ll_new <- sum(lse)
      W <- exp(logf - lse)
      if (is.finite(ll_new) && is.finite(ll) &&
          abs(ll_new - ll) < tol * (1 + abs(ll_new))) {
        ll <- ll_new
        break
      }
      ll <- ll_new
    }
    list(pi = pi_k, theta = theta, loglik = ll, iterations = iter)
  }

  # Anchored start: soft responsibilities from each subject's observed
  # category counts, which also aligns class k with category k. Extra starts
  # jitter the same anchor to probe local optima; best log-likelihood wins.
  starts <- lapply(seq_len(max(1L, n_starts)), function(s) {
    W0 <- cat_counts + 0.5
    if (s > 1L) W0 <- W0 * matrix(stats::rgamma(n * K, shape = 2, rate = 2), n, K)
    W0 / rowSums(W0)
  })
  fits <- lapply(starts, run_one)
  fits[[which.max(vapply(fits, `[[`, numeric(1), "loglik"))]]
}

ds_weight_matrix <- function(weight, values) {
  C <- length(values)
  if (weight == "nominal") return(diag(C))
  rng <- max(values) - min(values)
  if (weight == "quadratic") {
    if (rng^2 < 1e-9) return(diag(C))
    return(outer(values, values, function(a, b) 1 - (a - b)^2 / rng^2))
  }
  if (weight == "linear") {
    if (rng < 1e-9) return(diag(C))
    return(outer(values, values, function(a, b) 1 - abs(a - b) / rng))
  }
  stop("Unsupported ds_fiml weight: ", weight, call. = FALSE)
}

# Model-implied fixed-rater functionals, matching the saturated C++ map
# (estimate_fiml.cpp build_kappa_map): with disagreement L = 1 - W,
# pd averages pairwise disagreement over rater pairs, Conger's pe uses the
# product of per-rater marginals, Fleiss' pe the average marginal.
ds_implied_kappa <- function(fit, L) {
  R <- length(fit$theta)
  C <- ncol(fit$theta[[1L]])
  marg <- vapply(fit$theta, function(th) colSums(fit$pi * th), numeric(C))
  n_pairs <- R * (R - 1) / 2
  pd <- 0
  pec <- 0
  for (j in seq_len(R - 1L)) {
    for (l in (j + 1L):R) {
      P_jl <- t(fit$theta[[j]]) %*% (fit$pi * fit$theta[[l]])
      pd <- pd + sum(L * P_jl)
      pec <- pec + as.numeric(t(marg[, j]) %*% L %*% marg[, l])
    }
  }
  pd <- pd / n_pairs
  pec <- pec / n_pairs
  p_bar <- rowMeans(marg)
  pef <- as.numeric(t(p_bar) %*% L %*% p_bar)
  c(Conger = 1 - pd / pec, Fleiss = 1 - pd / pef)
}

# ---- inference ---------------------------------------------------------------
#
# Coordinate layout: eta = c(a, vec(B_1), ..., vec(B_R)) with a = log pi
# (length K) and B_j = log theta_j (K x C, column-major), so theta_j(k, c)
# sits at index K + (j-1)*K*C + (c-1)*K + k. pi = softmax(a); each row of
# theta_j is the softmax of the corresponding row of B_j.

ds_pack_eta <- function(pi_k, theta) {
  c(log(pi_k), unlist(lapply(theta, function(th) as.vector(log(th)))))
}

ds_unpack_eta <- function(eta, K, C, R) {
  softmax_rows <- function(B) {
    E <- exp(B - apply(B, 1L, max))
    E / rowSums(E)
  }
  a <- eta[seq_len(K)]
  pi_k <- exp(a - max(a))
  pi_k <- pi_k / sum(pi_k)
  theta <- vector("list", R)
  off <- K
  for (j in seq_len(R)) {
    theta[[j]] <- softmax_rows(matrix(eta[off + seq_len(K * C)], K, C))
    off <- off + K * C
  }
  list(pi = pi_k, theta = theta)
}

ds_posterior <- function(Y, pi_k, theta) {
  n <- nrow(Y)
  K <- length(pi_k)
  logf <- matrix(log(pi_k), n, K, byrow = TRUE)
  for (j in seq_len(ncol(Y))) {
    oj <- !is.na(Y[, j])
    logf[oj, ] <- logf[oj, ] + t(log(theta[[j]]))[Y[oj, j], , drop = FALSE]
  }
  mx <- apply(logf, 1L, max)
  lse <- mx + log(rowSums(exp(logf - mx)))
  list(W = exp(logf - lse), loglik_i = lse)
}

# Per-subject observed-data scores in the softmax coordinates, by Fisher's
# identity: with posterior responsibilities W,
#   s_i(a_k)    = W[i, k] - pi_k
#   s_i(b_jkc)  = W[i, k] * (1{Y_ij = c} - theta_j[k, c])   (0 if j missing).
# Rows of the b_jk blocks and the a block each sum to zero (the gauge).
ds_score_matrix <- function(Y, pi_k, theta, W = NULL) {
  n <- nrow(Y)
  R <- ncol(Y)
  K <- length(pi_k)
  C <- ncol(theta[[1L]])
  if (is.null(W)) W <- ds_posterior(Y, pi_k, theta)$W
  S <- matrix(0, n, K + R * K * C)
  S[, seq_len(K)] <- W - matrix(pi_k, n, K, byrow = TRUE)
  for (j in seq_len(R)) {
    oj <- which(!is.na(Y[, j]))
    if (!length(oj)) next
    Wo <- W[oj, , drop = FALSE]
    yo <- Y[oj, j]
    off <- K + (j - 1L) * K * C
    for (c in seq_len(C)) {
      cols <- off + (c - 1L) * K + seq_len(K)
      blk <- -Wo * matrix(theta[[j]][, c], length(oj), K, byrow = TRUE)
      hit <- yo == c
      if (any(hit)) blk[hit, ] <- blk[hit, ] + Wo[hit, , drop = FALSE]
      S[oj, cols] <- blk
    }
  }
  S
}

# Total score without materialising the n x p per-subject matrix; used inside
# the finite-difference information loop.
ds_score_total <- function(Y, pi_k, theta, W = NULL) {
  n <- nrow(Y)
  R <- ncol(Y)
  K <- length(pi_k)
  C <- ncol(theta[[1L]])
  if (is.null(W)) W <- ds_posterior(Y, pi_k, theta)$W
  tot <- numeric(K + R * K * C)
  tot[seq_len(K)] <- colSums(W) - n * pi_k
  for (j in seq_len(R)) {
    oj <- which(!is.na(Y[, j]))
    off <- K + (j - 1L) * K * C
    if (!length(oj)) next
    Wo <- W[oj, , drop = FALSE]
    yo <- Y[oj, j]
    Tk <- matrix(0, K, C)
    for (c in seq_len(C)) {
      hit <- yo == c
      if (any(hit)) Tk[, c] <- colSums(Wo[hit, , drop = FALSE])
    }
    tot[off + seq_len(K * C)] <- as.vector(Tk - colSums(Wo) * theta[[j]])
  }
  tot
}

ds_loglik_eta <- function(eta, Y, K, C, R) {
  par <- ds_unpack_eta(eta, K, C, R)
  sum(ds_posterior(Y, par$pi, par$theta)$loglik_i)
}

ds_score_total_eta <- function(eta, Y, K, C, R) {
  par <- ds_unpack_eta(eta, K, C, R)
  ds_score_total(Y, par$pi, par$theta)
}

# Observed information by central differences of the analytic score; FD error
# is O(h^2) and the score itself is exact. Symmetrised against FD noise.
ds_observed_info <- function(eta, Y, K, C, R, h = 1e-5) {
  p <- length(eta)
  Hc <- vapply(seq_len(p), function(r) {
    e <- numeric(p)
    e[r] <- h
    (ds_score_total_eta(eta + e, Y, K, C, R) -
       ds_score_total_eta(eta - e, Y, K, C, R)) / (2 * h)
  }, numeric(p))
  -(Hc + t(Hc)) / 2
}

# Rank-truncated spectral pseudo-inverse. The softmax gauge contributes
# 1 + R*K exact null directions that must be truncated, never inverted;
# anything truncated beyond the gauge is genuine weak identification and is
# surfaced through null_frac below.
ds_pinv <- function(M, rcond = 1e-6) {
  ei <- eigen(M, symmetric = TRUE)
  lam_max <- max(ei$values, 0)
  keep <- ei$values > rcond * lam_max
  V <- ei$vectors[, keep, drop = FALSE]
  list(pinv = V %*% (t(V) / ei$values[keep]),
       keep = keep, values = ei$values, vectors = ei$vectors)
}

ds_kappa_eta <- function(eta, L, K, C, R) {
  par <- ds_unpack_eta(eta, K, C, R)
  ds_implied_kappa(list(pi = par$pi, theta = par$theta), L)
}

# 2 x p Jacobian of the implied (Conger, Fleiss) map by central differences;
# the map is smooth and cheap, so FD is exact enough at h = 1e-6.
ds_kappa_jacobian <- function(eta, L, K, C, R, h = 1e-6) {
  p <- length(eta)
  vapply(seq_len(p), function(r) {
    e <- numeric(p)
    e[r] <- h
    (ds_kappa_eta(eta + e, L, K, C, R) -
       ds_kappa_eta(eta - e, L, K, C, R)) / (2 * h)
  }, numeric(2))
}

# Likelihood-side objects shared across weight matrices: per-subject scores,
# observed information, and its truncated pseudo-inverse (the bread).
ds_vcov_cache <- function(fit, Y, rcond = 1e-6) {
  K <- length(fit$pi)
  C <- ncol(fit$theta[[1L]])
  R <- length(fit$theta)
  eta <- ds_pack_eta(fit$pi, fit$theta)
  post <- ds_posterior(Y, fit$pi, fit$theta)
  S <- ds_score_matrix(Y, fit$pi, fit$theta, W = post$W)
  info <- ds_observed_info(eta, Y, K, C, R)
  pv <- ds_pinv(info, rcond)
  list(eta = eta, S = S, bread = pv$pinv, keep = pv$keep,
       values = pv$values, vectors = pv$vectors, K = K, C = C, R = R)
}

# Weight-specific delta method:
#   V_ml     = J bread J'                    (observed information only)
#   V_robust = J bread (sum_i s_i s_i') bread J'  = crossprod(psi)
# with per-subject influence contributions psi_i = J bread s_i. null_frac is
# the share of each kappa gradient lost to the truncated directions; the
# gauge lands there with zero overlap by construction, so values well above
# FD noise signal a weakly identified fit.
ds_kappa_vcov <- function(cache, L) {
  coefs <- c("Conger", "Fleiss")
  J <- ds_kappa_jacobian(cache$eta, L, cache$K, cache$C, cache$R)
  JB <- J %*% cache$bread
  V_ml <- JB %*% t(J)
  psi <- cache$S %*% t(JB)
  V_rb <- crossprod(psi)
  Vnull <- cache$vectors[, !cache$keep, drop = FALSE]
  proj <- J %*% Vnull
  null_frac <- rowSums(proj^2) / pmax(rowSums(J^2), .Machine$double.xmin)
  dimnames(V_ml) <- dimnames(V_rb) <- list(coefs, coefs)
  colnames(psi) <- coefs
  names(null_frac) <- coefs
  list(vcov_ml = V_ml, vcov_robust = V_rb, psi = psi, null_frac = null_frac)
}

fit_ds_fiml <- function(X, spec, weight, se = TRUE, rcond = 1e-6) {
  # mirrors the runner's score_values(); kept inline so this file is
  # self-contained for the probe
  vals <- sort(unique(as.integer(X[!is.na(X)]))) - 1L
  if (length(vals) < 2L) stop("Fewer than two observed categories.", call. = FALSE)
  cats <- vals + 1L
  Y <- matrix(match(X, cats), nrow(X), ncol(X))
  fit <- ds_em_fit(Y, C = length(cats))
  L <- 1 - ds_weight_matrix(weight, vals)
  out <- list(coefficients = ds_implied_kappa(fit, L), ds = fit)
  if (se) {
    cache <- ds_vcov_cache(fit, Y, rcond = rcond)
    vc <- ds_kappa_vcov(cache, L)
    out[c("vcov_ml", "vcov_robust", "psi", "null_frac")] <-
      vc[c("vcov_ml", "vcov_robust", "psi", "null_frac")]
  }
  structure(out, class = "ds_fiml_fit")
}

# Robust sandwich by default: correct for the pseudo-true value under
# misspecification and asymptotically equal to the ML flavour when the
# conditional-independence model holds.
vcov.ds_fiml_fit <- function(object, ...) {
  if (is.null(object$vcov_robust)) return(matrix(NA_real_, 0L, 0L))
  object$vcov_robust
}
