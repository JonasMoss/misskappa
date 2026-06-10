# Quadratic vector kappa via first and second moments.
#
# This is the vector analogue of kappa_quadratic_fiml(): for squared vector
# loss with a full feature metric W, Conger and Fleiss are smooth functions of
# the stacked rater-feature mean and covariance. The pairwise path is the
# pairwise-available covariance plug-in; the FIML path reuses the saturated
# normal EM / score / information helpers from alpha_continuous.R.

.kvq_flatten <- function(x) {
  dims <- dim(x)
  if (length(dims) != 3L) {
    stop("'x' must be a subjects-by-raters-by-features array.")
  }
  if (dims[1L] < 1L || dims[2L] < 2L || dims[3L] < 1L) {
    stop("'x' must have at least one subject, two raters, and one feature.")
  }
  if (!is.numeric(x)) stop("'x' must be numeric.")
  storage.mode(x) <- "double"

  n <- dims[1L]; R <- dims[2L]; p <- dims[3L]
  flat <- matrix(NA_real_, nrow = n, ncol = R * p)
  for (r in seq_len(R)) {
    for (l in seq_len(p)) {
      flat[, (r - 1L) * p + l] <- x[, r, l]
    }
  }
  list(X = flat, R = R, features = p)
}

.kvq_validate_W <- function(W, p) {
  if (is.null(W)) W <- diag(p)
  W <- as.matrix(W)
  if (!is.numeric(W) || any(dim(W) != c(p, p)) || any(!is.finite(W))) {
    stop("'W' must be a finite numeric features-by-features matrix.")
  }
  if (max(abs(W - t(W))) > 1e-10) stop("'W' must be symmetric.")
  W <- (W + t(W)) / 2
  ev <- eigen(W, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) < -1e-8) stop("'W' must be positive semidefinite.")
  W
}

.kvq_design <- function(R, p, W) {
  one <- matrix(1, R, R)
  P <- diag(R) - one / R
  list(
    B = kronecker(one, W),
    T = kronecker(diag(R), W),
    G = kronecker(P, W)
  )
}

.kvq_vech_grad <- function(A) {
  idx <- which(lower.tri(A, diag = TRUE), arr.ind = TRUE)
  ifelse(idx[, 1L] == idx[, 2L], A[idx], A[idx] + A[cbind(idx[, 2L], idx[, 1L])])
}

.kvq_grad <- function(mu, Sigma, R, p, W) {
  q <- R * p
  pstar <- q * (q + 1L) / 2L
  des <- .kvq_design(R, p, W)

  tB <- sum(des$B * Sigma)
  tT <- sum(des$T * Sigma)
  tG <- as.numeric(t(mu) %*% des$G %*% mu)

  dC <- (R - 1) * tT + R * tG
  nC <- tB - tT
  dF <- (R - 1) * (tT + tG)
  nF <- tB - tT - tG
  conger <- nC / dC
  fleiss <- nF / dF

  At <- matrix(0, 3L, q + pstar)
  At[1L, (q + 1L):(q + pstar)] <- .kvq_vech_grad(des$B)
  At[2L, (q + 1L):(q + pstar)] <- .kvq_vech_grad(des$T)
  At[3L, seq_len(q)] <- 2 * (des$G %*% mu)

  Jh <- matrix(0, 2L, 3L)
  Jh[1L, ] <- c(1 / dC, -1 / dC - nC * (R - 1) / dC^2, -nC * R / dC^2)
  f_t23 <- -1 / dF - nF * (R - 1) / dF^2
  Jh[2L, ] <- c(1 / dF, f_t23, f_t23)

  G <- Jh %*% At
  rownames(G) <- c("Conger", "Fleiss")
  list(
    estimates = c(Conger = conger, Fleiss = fleiss),
    G = G,
    Jh = Jh,
    design = des,
    summaries = c(B = tB, T = tT, G = tG)
  )
}

.kvq_pairwise <- function(X, R, p, W) {
  keep <- rowSums(is.finite(X)) > 0L
  if (!all(keep)) X <- X[keep, , drop = FALSE]
  n <- nrow(X)
  q <- ncol(X)
  M <- is.finite(X)
  if (any(colSums(M) == 0L)) {
    stop("every rater-feature cell must be observed for at least one subject.")
  }

  X0 <- X
  X0[!M] <- 0
  count <- colSums(M)
  mu <- colSums(X0) / count
  Y <- sweep(X, 2L, mu, "-")
  Y[!M] <- 0

  p1 <- count / n
  p2 <- crossprod(M) / n
  if (any(p2 <= 0)) stop("every rater-feature pair needs at least one overlap.")
  Sigma <- matrix(0, q, q)
  for (a in seq_len(q)) {
    for (b in a:q) {
      ok <- M[, a] & M[, b]
      Sigma[a, b] <- Sigma[b, a] <- mean(Y[ok, a] * Y[ok, b])
    }
  }

  grad <- .kvq_grad(mu, Sigma, R, p, W)
  des <- grad$design
  phi <- matrix(0, n, 3L)
  colnames(phi) <- c("B", "T", "G")
  g_mu <- 2 * (des$G %*% mu)
  for (i in seq_len(n)) {
    cov_if <- matrix(0, q, q)
    for (a in seq_len(q)) {
      if (!M[i, a]) next
      phi[i, 3L] <- phi[i, 3L] + g_mu[a] * Y[i, a] / p1[a]
      for (b in seq_len(q)) {
        if (!M[i, b]) next
        cov_if[a, b] <- (Y[i, a] * Y[i, b] - Sigma[a, b]) / p2[a, b]
      }
    }
    phi[i, 1L] <- sum(des$B * cov_if)
    phi[i, 2L] <- sum(des$T * cov_if)
  }

  psi <- phi %*% t(grad$Jh)
  colnames(psi) <- c("Conger", "Fleiss")
  vcov_mat <- crossprod(psi) / n^2
  dimnames(vcov_mat) <- list(colnames(psi), colnames(psi))
  list(mu = mu, Sigma = Sigma, estimates = grad$estimates,
       vcov = vcov_mat, psi = psi, summaries = grad$summaries)
}

.kvq_fiml <- function(X, R, p, W, em_options = list()) {
  X <- X[rowSums(is.finite(X)) > 0L, , drop = FALSE]
  Rmiss <- is.finite(X)
  n <- nrow(X)
  q <- ncol(X)
  if (any(colSums(Rmiss) == 0L)) {
    stop("every rater-feature cell must be observed for at least one subject.")
  }
  opt <- utils::modifyList(
    list(tol = 1e-8, max_iter = 10000L, fd_h = 1e-5), em_options
  )
  pstar <- q * (q + 1L) / 2L
  vech_pos <- .amc_vech_pos(q)
  patterns <- .amc_patterns(Rmiss)
  em <- .amc_em(X, patterns, tol = opt$tol, max_iter = opt$max_iter)
  if (!em$converged) {
    warning("saturated EM did not converge in ", opt$max_iter, " iterations.")
  }

  grad <- .kvq_grad(em$mu, em$Sigma, R, p, W)
  theta <- c(em$mu, .amc_vech(em$Sigma))
  scores <- .amc_score_matrix(em$mu, em$Sigma, X, patterns, vech_pos, q, pstar)
  H <- .amc_information(theta, X, patterns, vech_pos, q, pstar, h = opt$fd_h)
  Hinv <- solve(H)

  psi <- scores %*% (Hinv %*% t(grad$G))
  colnames(psi) <- c("Conger", "Fleiss")
  vcov_mat <- crossprod(psi) / n^2
  dimnames(vcov_mat) <- list(colnames(psi), colnames(psi))
  list(mu = em$mu, Sigma = em$Sigma, estimates = grad$estimates,
       vcov = vcov_mat, psi = psi, summaries = grad$summaries,
       iterations = em$iterations, converged = em$converged)
}

#' Quadratic vector agreement from covariance moments
#'
#' @description
#' Internal backend for vector-valued quadratic agreement. For a
#' subjects-by-raters-by-features array, estimates Conger and Fleiss kappas
#' generated by the squared vector loss `(x - y)' W (x - y)`, where `W` is a
#' full symmetric positive-semidefinite feature-weight matrix. The `"pairwise"`
#' method uses pairwise-available covariance moments and is MCAR-oriented; the
#' `"nt_fiml"` method fits the saturated normal mean/covariance by EM and is
#' the vector analogue of [kappa_quadratic_fiml()].
#'
#' Both methods need each rater-feature cell to be observed at least once. The
#' `"pairwise"` method additionally requires every rater-feature pair to
#' overlap at least once, because each covariance entry is estimated from
#' directly co-observed rows. Direct `"nt_fiml"` callers should enforce the
#' same complete pairwise co-observation condition when the saturated
#' covariance functional is the target.
#'
#' @param x Numeric array with dimensions subjects, raters, features.
#' @param method `"pairwise"` or `"nt_fiml"`.
#' @param W Optional features-by-features symmetric positive-semidefinite
#'   weight matrix. Defaults to the identity matrix.
#' @param em_options Used only by `"nt_fiml"`; named list with `tol`,
#'   `max_iter`, and `fd_h`.
#'
#' @return A `misskappa_estimate` object with `Conger` and `Fleiss`
#'   coefficients, covariance matrix, per-subject influence functions, and
#'   fitted mean/covariance moments.
#'
#' @keywords internal
kappa_vector_quadratic <- function(x,
                                   method = c("pairwise", "nt_fiml"),
                                   W = NULL,
                                   em_options = list()) {
  method <- match.arg(method)
  prep <- .kvq_flatten(x)
  W <- .kvq_validate_W(W, prep$features)
  fit <- if (method == "pairwise") {
    .kvq_pairwise(prep$X, prep$R, prep$features, W)
  } else {
    .kvq_fiml(prep$X, prep$R, prep$features, W, em_options = em_options)
  }

  structure(
    list(
      estimates = fit$estimates,
      vcov = fit$vcov,
      psi = fit$psi,
      method = paste0("vector-quadratic-", method),
      weight = "quadratic",
      W = W,
      moments = list(mu = fit$mu, Sigma = fit$Sigma,
                     summaries = fit$summaries,
                     iterations = fit$iterations,
                     converged = fit$converged)
    ),
    class = "misskappa_estimate"
  )
}
