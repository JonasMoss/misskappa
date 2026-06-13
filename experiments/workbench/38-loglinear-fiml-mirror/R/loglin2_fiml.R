# All-two-way-interaction log-linear FIML: the categorical mirror of nt-fiml.
#
# Model: the joint law of the R rater ratings over {1..C}^R is the
# maximum-entropy (I-projection) distribution matching all univariate and
# bivariate margins -- equivalently the log-linear model with rater main
# effects and all pairwise interactions but no three-way-or-higher terms:
#
#   log p(y_1, ..., y_R) = const + sum_j u_j(y_j) + sum_{j<l} u_jl(y_j, y_l).
#
# Its sufficient statistics are exactly the 1-way and 2-way margins, which are
# exactly the statistics nominal Conger/Fleiss/BP kappa depend on (study 30
# margins-pinning; the ds_implied_kappa functional). So it is the discrete
# analogue of the saturated Gaussian (saturated in means + covariances,
# max-entropy beyond), and should inherit nt-fiml's duality: model-free
# consistency under MCAR, and consistency + efficiency under MAR when the
# truth has no genuine 3-way-or-higher rater interactions.
#
# Under item-level missingness the observed-data likelihood is fit by EM:
# the E-step distributes each observed pattern over its compatible full-table
# completions using the current fitted table; the M-step is an IPF projection
# of the expected full table onto the all-2-way model, which stats::loglin()
# performs directly. The full C^R table is materialised (as in saturated
# cat_fiml), so this is a moderate-C,R prototype.
#
# Point estimates only -- this probe tests the bias/efficiency duality, not
# SEs. If the mirror holds, delta-method SEs follow the validated DS-FIML
# recipe ([[project_ds_fiml_prototype]]).
#
# Sourced by the probe-38 runner; keep free of dependencies on the surrounding
# script.

# Disagreement matrix L = 1 - W from a weight name and integer score values,
# matching the study-29 ds_weight_matrix() convention.
loglin2_weight_matrix <- function(weight, values) {
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
  stop("Unsupported loglin2 weight: ", weight, call. = FALSE)
}

# Per-rater marginals (C x R) and pairwise joints from a full joint table p
# (array with dim rep(C, R)).
loglin2_margins <- function(p) {
  R <- length(dim(p))
  C <- dim(p)[1L]
  marg <- vapply(seq_len(R), function(j) as.numeric(apply(p, j, sum)), numeric(C))
  list(marg = marg, R = R, C = C)
}

# Model-implied fixed-rater (Conger, Fleiss) from a full joint table p and a
# disagreement matrix L. Identical functional to ds_implied_kappa, fed the
# log-linear-fitted joint instead of the latent-class-implied joint.
loglin2_implied_kappa <- function(p, L) {
  m <- loglin2_margins(p)
  R <- m$R
  marg <- m$marg
  n_pairs <- R * (R - 1) / 2
  pd <- 0
  pec <- 0
  for (j in seq_len(R - 1L)) {
    for (l in (j + 1L):R) {
      P_jl <- apply(p, c(j, l), sum)              # C x C bivariate joint
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

# Closed-form truth helper: implied kappa directly from per-rater marginals
# (C x R matrix `marg`) and a list of pairwise joints `P[[j]][[l]]` (C x C),
# for DGPs whose true margins are known analytically.
loglin2_kappa_from_margins <- function(marg, P, L) {
  R <- ncol(marg)
  n_pairs <- R * (R - 1) / 2
  pd <- 0
  pec <- 0
  for (j in seq_len(R - 1L)) {
    for (l in (j + 1L):R) {
      pd <- pd + sum(L * P[[j]][[l]])
      pec <- pec + as.numeric(t(marg[, j]) %*% L %*% marg[, l])
    }
  }
  pd <- pd / n_pairs
  pec <- pec / n_pairs
  p_bar <- rowMeans(marg)
  pef <- as.numeric(t(p_bar) %*% L %*% p_bar)
  c(Conger = 1 - pd / pec, Fleiss = 1 - pd / pef)
}

# Independence-model init: outer product of available-case 1-way margins.
loglin2_independence_table <- function(Y, C, R) {
  margs <- lapply(seq_len(R), function(j) {
    tab <- tabulate(Y[!is.na(Y[, j]), j], nbins = C) + 0.5  # Laplace floor
    tab / sum(tab)
  })
  p <- margs[[1L]]
  for (j in 2:R) p <- outer(p, margs[[j]])
  array(p, dim = rep(C, R))
}

# All pairwise margins of a full table p, stacked into one vector -- the
# log-linear sufficient statistics, and exactly the quantities kappa depends
# on, so EM convergence is monitored here rather than on the full-table
# log-likelihood (same linear rate, but this is the estimand-relevant target
# and avoids log(0) sensitivity on sparse tables).
loglin2_pair_margins <- function(p, all_pairs) {
  unlist(lapply(all_pairs, function(pr) as.numeric(apply(p, pr, sum))),
         use.names = FALSE)
}

# EM for the all-2-way log-linear model. Y is n x R with codes in 1..C, NA
# missing. Returns the fitted full table p (array dim rep(C,R)) and iteration
# count. Convergence is on the pairwise margins (tol). Complete data needs no
# EM: the expected table is the empirical table, so one IPF projection suffices.
loglin2_em_fit <- function(Y, C, max_iter = 300L, tol = 1e-7) {
  R <- ncol(Y)
  all_pairs <- utils::combn(R, 2L, simplify = FALSE)

  # Group subjects by observed pattern (NA-positions + observed values).
  keys <- apply(Y, 1L, function(r) paste(ifelse(is.na(r), "NA", r), collapse = "|"))
  grp <- split(seq_len(nrow(Y)), keys)
  patterns <- lapply(grp, function(rows) Y[rows[1L], ])
  pattern_n <- vapply(grp, length, integer(1))
  idx_list <- lapply(patterns, function(row) {
    lapply(seq_len(R), function(j) if (is.na(row[j])) seq_len(C) else as.integer(row[j]))
  })
  complete <- !any(is.na(Y))

  project <- function(E, start) {
    if (is.null(start)) start <- array(1, dim = rep(C, R))
    # suppressWarnings: loglin warns when its inner IPF hits the iteration cap
    # on sparse tables; the EM still converges on the margins (the convergence
    # criterion below), so the per-call warning is noise at simulation scale.
    fit <- suppressWarnings(stats::loglin(E, margin = all_pairs, start = start,
                                          fit = TRUE, print = FALSE,
                                          eps = 1e-9, iter = 400L))
    fit$fit / sum(fit$fit)
  }

  # Complete data: the expected full table is exactly the empirical count
  # table, so a single IPF projection is the MLE -- no EM loop.
  if (complete) {
    E <- array(0, dim = rep(C, R))
    for (g in seq_along(idx_list)) {
      idx <- idx_list[[g]]
      E <- do.call(`[<-`, c(list(E), idx, list(value = pattern_n[g])))
    }
    return(list(p = project(E, NULL), iterations = 1L))
  }

  dims <- rep(C, R)
  # One EM map p -> project(E-step(p)). Operates on a flat probability vector.
  em_map <- function(pv) {
    p <- array(pv, dim = dims)
    E <- array(0, dim = dims)
    for (g in seq_along(idx_list)) {
      idx <- idx_list[[g]]
      sub <- do.call(`[`, c(list(p), idx, list(drop = FALSE)))
      add <- sub * (pattern_n[g] / sum(sub))
      cur <- do.call(`[`, c(list(E), idx, list(drop = FALSE)))
      E <- do.call(`[<-`, c(list(E), idx, list(value = cur + add)))
    }
    as.numeric(project(E, pmax(p, 1e-12)))
  }
  clamp <- function(pv) { pv <- pmax(pv, 0); pv / sum(pv) }

  # SQUAREM (SqS3) acceleration of the EM map; the plain EM is linearly
  # convergent and slow on sparse tables, so a stabilised quadratic step cuts
  # the iteration count several-fold. Convergence is on the pairwise margins.
  p0 <- as.numeric(loglin2_independence_table(Y, C, R))
  m_prev <- loglin2_pair_margins(array(p0, dims), all_pairs)
  iter <- 0L
  for (sweep in seq_len(max_iter)) {
    p1 <- em_map(p0)
    p2 <- em_map(p1)
    iter <- iter + 2L
    r <- p1 - p0
    v <- (p2 - p1) - r
    sr <- sqrt(sum(r^2)); sv <- sqrt(sum(v^2))
    if (sv < 1e-14) {
      p0 <- p2
    } else {
      alpha <- min(-(sr / sv), -1)               # SqS3 step length, <= -1
      p_try <- clamp(p0 - 2 * alpha * r + alpha^2 * v)
      p_acc <- em_map(p_try)                       # stabilising EM step
      iter <- iter + 1L
      # Safeguard: a too-aggressive extrapolation can leave the simplex and
      # produce a non-finite map; fall back to the plain double-EM step p2.
      p0 <- if (all(is.finite(p_acc))) p_acc else p2
    }
    m_new <- loglin2_pair_margins(array(p0, dims), all_pairs)
    if (max(abs(m_new - m_prev)) < tol) { m_prev <- m_new; break }
    m_prev <- m_new
  }
  list(p = array(p0, dims), iterations = iter)
}

# Convenience wrapper: fit and map to (Conger, Fleiss). X is n x R with integer
# category codes (any labels) and NA; values are the integer score values used
# for weighting.
fit_loglin2_fiml <- function(X, weight = "nominal", values = NULL,
                             max_iter = 200L, tol = 1e-6) {
  vals <- if (is.null(values)) sort(unique(as.integer(X[!is.na(X)]))) - 1L else values
  C <- length(vals)
  cats <- vals + 1L
  Y <- matrix(match(X, cats), nrow(X), ncol(X))
  fit <- loglin2_em_fit(Y, C = C, max_iter = max_iter, tol = tol)
  L <- 1 - loglin2_weight_matrix(weight, vals)
  structure(list(coefficients = loglin2_implied_kappa(fit$p, L), fit = fit),
            class = "loglin2_fiml_fit")
}
