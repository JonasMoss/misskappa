# Pure-R reference implementation of the saturated normal FIML used as a test
# oracle for the C++ normal-FIML backend (alpha and quadratic kappa): the EM
# for (mu, Sigma), the analytic casewise scores, a finite-difference observed
# information, and the quadratic-kappa gradient. These used to live in the
# package namespace; they are test-only.
#
# Conventions match the C++ backend: the MLE covariance divides by n, and
# vech() is the column-major lower triangle including the diagonal, with
# off-diagonal scores carrying the factor 2 of the symmetric pair.

# ---- vech bookkeeping -------------------------------------------------------

# Column-major lower-triangular half-vectorisation, matching magmaan's
# `for (j) for (i in j:p)` order.
.amc_vech <- function(S) S[lower.tri(S, diag = TRUE)]

# p x p integer matrix giving, for i >= j, the position of (i, j) in vech();
# symmetric for convenient lookup.
.amc_vech_pos <- function(p) {
  pos <- matrix(0L, p, p)
  idx <- which(lower.tri(pos, diag = TRUE), arr.ind = TRUE) # column-major
  for (k in seq_len(nrow(idx))) {
    i <- idx[k, 1L]; j <- idx[k, 2L]
    pos[i, j] <- k
    pos[j, i] <- k
  }
  pos
}

# ---- saturated normal EM ----------------------------------------------------

# Missing-data patterns: list of {rows, o} where o is the observed column set.
.amc_patterns <- function(R) {
  key <- apply(R, 1L, function(r) paste0(as.integer(r), collapse = ""))
  lapply(split(seq_len(nrow(R)), key), function(rows) {
    list(rows = rows, o = which(R[rows[1L], ]))
  })
}

# EM for the saturated multivariate-normal mean and covariance (MLE, /n).
.amc_em <- function(X, patterns, tol = 1e-8, max_iter = 10000L) {
  n <- nrow(X); p <- ncol(X)
  mu <- colMeans(X, na.rm = TRUE)
  v <- apply(X, 2L, stats::var, na.rm = TRUE)
  v[!is.finite(v) | v <= 0] <- 1
  Sigma <- diag(v, p)
  converged <- FALSE
  iter <- 0L
  for (iter in seq_len(max_iter)) {
    T1 <- numeric(p)
    T2 <- matrix(0, p, p)
    for (g in patterns) {
      o <- g$o; rows <- g$rows; ng <- length(rows)
      m <- setdiff(seq_len(p), o)
      Xo <- X[rows, o, drop = FALSE]
      Xc <- matrix(0, ng, p)
      Xc[, o] <- Xo
      if (length(m) > 0L) {
        Soo_inv <- solve(Sigma[o, o, drop = FALSE])
        B <- Sigma[m, o, drop = FALSE] %*% Soo_inv          # regress m on o
        Xo_c <- sweep(Xo, 2L, mu[o], "-")
        Xc[, m] <- sweep(Xo_c %*% t(B), 2L, mu[m], "+")     # E[x_m | x_o]
        Cmm <- Sigma[m, m, drop = FALSE] -
          B %*% Sigma[o, m, drop = FALSE]                   # residual cov
        T2[m, m] <- T2[m, m] + ng * Cmm
      }
      T1 <- T1 + colSums(Xc)
      T2 <- T2 + crossprod(Xc)
    }
    mu_new <- T1 / n
    Sigma_new <- T2 / n - tcrossprod(mu_new)
    Sigma_new <- (Sigma_new + t(Sigma_new)) / 2
    delta <- max(abs(mu_new - mu), abs(Sigma_new - Sigma))
    mu <- mu_new
    Sigma <- Sigma_new
    if (delta < tol) { converged <- TRUE; break }
  }
  list(mu = mu, Sigma = Sigma, iterations = iter, converged = converged)
}

# ---- casewise scores --------------------------------------------------------

# n x q matrix of casewise scores of the observed-data Gaussian log-likelihood
# wrt theta = (mu, vech(Sigma)). Analytic; the off-diagonal vech score carries
# the factor 2 of the symmetric pair, matching .amc_alpha_grad_s().
.amc_score_matrix <- function(mu, Sigma, X, patterns, vech_pos, p, pstar) {
  n <- nrow(X)
  S <- matrix(0, n, p + pstar)
  for (g in patterns) {
    o <- g$o; rows <- g$rows
    no <- length(o)
    if (no == 0L) next
    Soinv <- solve(Sigma[o, o, drop = FALSE])
    E <- sweep(X[rows, o, drop = FALSE], 2L, mu[o], "-")
    U <- E %*% Soinv                                        # rows are (Soinv e)'
    S[rows, o] <- S[rows, o] + U                            # d/dmu
    half_const <- 0.5 * Soinv
    for (a in seq_len(no)) {
      for (b in seq_len(a)) {
        wt <- if (a == b) 1 else 2
        col <- p + vech_pos[o[a], o[b]]
        S[rows, col] <- S[rows, col] +
          wt * (0.5 * U[, a] * U[, b] - half_const[a, b])
      }
    }
  }
  S
}

.amc_unpack <- function(theta, p, pstar) {
  Sig <- matrix(0, p, p)
  Sig[lower.tri(Sig, diag = TRUE)] <- theta[(p + 1L):(p + pstar)]
  Sig <- Sig + t(Sig) - diag(diag(Sig))
  list(mu = theta[seq_len(p)], Sigma = Sig)
}

# Observed information per case, H = -(1/n) d(sum_i score_i)/dtheta, by central
# finite differences on the analytic total score. Sidesteps duplication-matrix
# bookkeeping; H and the sandwich are then exact functions of the analytic score.
.amc_information <- function(theta, X, patterns, vech_pos, p, pstar, h = 1e-5) {
  n <- nrow(X); q <- p + pstar
  g_of <- function(th) {
    pr <- .amc_unpack(th, p, pstar)
    colSums(.amc_score_matrix(pr$mu, pr$Sigma, X, patterns, vech_pos, p, pstar))
  }
  H <- matrix(0, q, q)
  for (k in seq_len(q)) {
    tp <- theta; tp[k] <- tp[k] + h
    tm <- theta; tm[k] <- tm[k] - h
    H[, k] <- -(g_of(tp) - g_of(tm)) / (2 * h) / n
  }
  (H + t(H)) / 2
}

# ---- coefficient gradient ---------------------------------------------------

# Conger, Fleiss, and Brennan-Prediger quadratic kappas as smooth maps of
# theta = (mu, vech Sigma), following kappa-missing-quadratic.tex (eqs.
# fleiss-t / conger-t). Returns the point estimates and the 3 x q gradient
# G = J_h A_t, in the same column-major vech convention (off-diagonal weight 2)
# as .amc_score_matrix(), so the sandwich contraction scores %*% (Hinv %*% t(G))
# is exact.
#
# Summaries: t1 = 1' Sigma 1, t2 = tr Sigma, t3 = mu' Q mu with Q = I - R^-1 11'.
#   kappa_F = (t1 - t2 - t3) / [(R-1)(t2 + t3)]
#   kappa_C = (t1 - t2)      / [(R-1) t2 + R t3]
.kqf_grad <- function(mu, Sigma, bp_chance) {
  R <- length(mu)
  pstar <- R * (R + 1L) / 2L
  q <- R + pstar

  t1 <- sum(Sigma)            # 1' Sigma 1
  t2 <- sum(diag(Sigma))      # tr Sigma
  mubar <- mean(mu)
  t3 <- sum((mu - mubar)^2)   # mu' Q mu

  dC <- (R - 1) * t2 + R * t3 # Conger denominator
  nC <- t1 - t2               # Conger numerator
  dF <- (R - 1) * (t2 + t3)   # Fleiss denominator
  nF <- t1 - t2 - t3          # Fleiss numerator
  d_obs <- (2 / (R - 1)) * (t2 + t3 - t1 / R)
  conger <- nC / dC
  fleiss <- nF / dF
  bp <- 1 - d_obs / bp_chance

  # A_t = d(t1, t2, t3) / d theta, a 3 x q matrix with columns
  # [mu (R) | vech(Sigma) (pstar)] in column-major lower-triangular order.
  idx <- which(lower.tri(Sigma, diag = TRUE), arr.ind = TRUE) # column-major
  diagonal <- idx[, 1L] == idx[, 2L]
  At <- matrix(0, 3L, q)
  At[1L, (R + 1L):q] <- ifelse(diagonal, 1, 2)  # t1: off-diagonals enter twice
  At[2L, (R + 1L):q] <- ifelse(diagonal, 1, 0)  # t2: diagonal only
  At[3L, seq_len(R)] <- 2 * (mu - mubar)        # t3: mean block, 2 Q mu

  # J_h = d(kappa_C, kappa_F, kappa_BP) / d(t1, t2, t3), a 3 x 3 matrix.
  Jh <- matrix(0, 3L, 3L)
  Jh[1L, ] <- c(1 / dC, -1 / dC - nC * (R - 1) / dC^2, -nC * R / dC^2)
  f_t23 <- -1 / dF - nF * (R - 1) / dF^2        # equal for t2 and t3
  Jh[2L, ] <- c(1 / dF, f_t23, f_t23)
  bp_factor <- -2 / (bp_chance * (R - 1))
  Jh[3L, ] <- c(-bp_factor / R, bp_factor, bp_factor)

  G <- Jh %*% At
  rownames(G) <- c("Conger", "Fleiss", "Brennan-Prediger")
  list(
    estimates = c(Conger = conger, Fleiss = fleiss,
                  `Brennan-Prediger` = bp),
    G = G,
    summaries = c(t1 = t1, t2 = t2, t3 = t3)
  )
}
