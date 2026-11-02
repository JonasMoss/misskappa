# Quadratic vector kappa via first and second moments.
#
# This is the vector analogue of .kappa_quadratic_fiml(): for squared vector
# loss with a full feature metric W, Conger and Fleiss are smooth functions of
# the stacked rater-feature mean and covariance. The pairwise path is the
# pairwise-available covariance plug-in; the FIML path reuses the saturated
# normal EM / score / information machinery of the C++ normal-FIML backend.

.kvq_flatten <- function(x) {
  dims <- dim(x)
  if (length(dims) != 3L) {
    stop("'x' must be a subjects-by-raters-by-features array.", call. = FALSE)
  }
  if (dims[1L] < 1L || dims[2L] < 2L || dims[3L] < 1L) {
    stop("'x' must have at least one subject, two raters, and one feature.",
         call. = FALSE)
  }
  if (!is.numeric(x)) stop("'x' must be numeric.", call. = FALSE)
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
    stop("'W' must be a finite numeric features-by-features matrix.",
         call. = FALSE)
  }
  if (max(abs(W - t(W))) > 1e-10) stop("'W' must be symmetric.", call. = FALSE)
  W <- (W + t(W)) / 2
  ev <- eigen(W, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) < -1e-8) stop("'W' must be positive semidefinite.", call. = FALSE)
  W
}

.kvq_validate_ranges <- function(ranges, p) {
  if (is.null(ranges)) return(numeric(0))
  ranges <- as.matrix(ranges)
  if (!is.numeric(ranges) || nrow(ranges) != p || ncol(ranges) != 2L ||
      any(!is.finite(ranges))) {
    stop("'ranges' must be a finite features-by-2 matrix of (lo, hi) per feature.",
         call. = FALSE)
  }
  widths <- ranges[, 2L] - ranges[, 1L]
  if (any(widths <= 0)) stop("'ranges' must have hi > lo for every feature.",
                             call. = FALSE)
  as.numeric(widths)
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
    stop("every rater-feature cell must be observed for at least one subject.",
         call. = FALSE)
  }

  X0 <- X
  X0[!M] <- 0
  count <- colSums(M)
  mu <- colSums(X0) / count
  Y <- sweep(X, 2L, mu, "-")
  Y[!M] <- 0

  p1 <- count / n
  p2 <- crossprod(M) / n
  if (any(p2 <= 0)) stop("every rater-feature pair needs at least one overlap.",
                         call. = FALSE)
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

.kvq_fiml <- function(X, R, p, W, em_options = list(), widths = numeric(0)) {
  X <- X[rowSums(is.finite(X)) > 0L, , drop = FALSE]
  Rmiss <- is.finite(X)
  if (any(colSums(Rmiss) == 0L)) {
    stop("every rater-feature cell must be observed for at least one subject.",
         call. = FALSE)
  }
  opt <- utils::modifyList(
    list(tol = 1e-8, max_iter = 100000L), em_options
  )
  covariance <- .normal_covariance_model(opt)
  mean_model <- .normal_mean_model(opt)
  # Vector-valued FIML supports the saturated covariance and the block
  # compound-symmetry covariance Sigma = J_R (x) Psi + I_R (x) E (selected by the
  # feature block size in C++). The equal-mean model is not yet generalized.
  if (mean_model != "free") {
    stop("equal-mean normal FIML is not yet implemented for ",
         "vector-valued normal-FIML kappa.", call. = FALSE)
  }
  opt$covariance <- covariance
  opt$mean <- mean_model
  raw <- rcpp_kappa_vector_quadratic_fiml(X, p, W, opt, widths)
  nm <- c("Conger", "Fleiss", "Brennan-Prediger")
  estimates <- as.numeric(raw$estimates)
  names(estimates) <- nm
  grad <- .kvq_grad(as.numeric(raw$mu), raw$Sigma, R, p, W)
  psi <- raw$psi
  colnames(psi) <- nm
  vcov_mat <- raw$vcov
  dimnames(vcov_mat) <- list(nm, nm)
  list(mu = as.numeric(raw$mu), Sigma = raw$Sigma, estimates = estimates,
       vcov = vcov_mat, psi = psi, summaries = grad$summaries,
       iterations = as.integer(raw$iterations), converged = isTRUE(raw$converged))
}

#' Quadratic vector agreement from covariance moments
#'
#' @description
#' Internal backend for vector-valued quadratic agreement. For a
#' subjects-by-raters-by-features array, estimates Conger, Fleiss, and
#' Brennan-Prediger kappas generated by the squared vector loss
#' `(x - y)' W (x - y)`, where `W` is a full symmetric positive-semidefinite
#' feature-weight matrix. The `"pairwise"` method uses pairwise-available
#' covariance moments and is MCAR-oriented (Conger and Fleiss only). The
#' `"nt_fiml"` method fits the saturated normal mean/covariance by EM and is
#' the vector analogue of [.kappa_quadratic_fiml()]. The `"nt_fiml_cs"` method
#' fits the same normal FIML under the block compound-symmetry covariance
#' `Sigma = J_R \%x\% Psi + I_R \%x\% E`. With a quadratic weight this
#' is the covariance/variance-components route used for Lin's concordance
#' correlation coefficient under unbalanced repeated measures. The full matrix
#' `W` couples features in the loss, separate from the missing-data machinery.
#'
#' Both methods need each rater-feature cell to be observed at least once. The
#' `"pairwise"` method additionally requires every rater-feature pair to
#' overlap at least once, because each covariance entry is estimated from
#' directly co-observed rows. Direct `"nt_fiml"` callers should enforce the
#' same complete pairwise co-observation condition when the saturated
#' covariance functional is the target.
#'
#' The full-`W` route is internal because the current normal-FIML backend is
#' numerically strict: it inverts per-pattern covariance blocks and can return a
#' numerical error when the fitted covariance is singular or nearly singular
#' (for example, a deterministic or non-varying feature within a missingness
#' pattern). There is no ridge or rank-aware fallback yet. The block
#' compound-symmetry route (`"nt_fiml_cs"`) is the structured counterpart and
#' needs only a connected co-observation graph.
#'
#' @param x Numeric array with dimensions subjects, raters, features.
#' @param method `"pairwise"`, `"nt_fiml"`, or `"nt_fiml_cs"`.
#' @param W Optional features-by-features symmetric positive-semidefinite
#'   weight matrix. Defaults to the identity matrix.
#' @param em_options Used by `"nt_fiml"` / `"nt_fiml_cs"`. A named list with `tol`
#'   and `max_iter`. The covariance structure is selected by `method` (the
#'   estimator name), so `mean` must be `"free"`. The exchangeable equal-mean
#'   model has not been generalized to vector-valued ratings.
#' @param ranges Optional features-by-2 matrix of `(lo, hi)` score bounds, one
#'   row per feature, that fixes the Brennan-Prediger chance baseline
#'   `sum_a W_aa (hi_a - lo_a)^2 / 6` instead of estimating each width from the
#'   observed data range. Fixing the range makes the chance denominator
#'   deterministic, so the Brennan-Prediger standard error is exact and the
#'   estimand no longer drifts with the sample range. It is the vector analogue
#'   of the scalar `values` argument. `NULL` (default) keeps the observed-range
#'   behaviour. Used by the `"nt_fiml"` / `"nt_fiml_cs"` methods only (the
#'   `"pairwise"` method does not report Brennan-Prediger).
#'
#' @return A `misskappa_estimate` object with `Conger`, `Fleiss`, and
#'   `Brennan-Prediger` coefficients (`"pairwise"` reports `Conger` and `Fleiss`
#'   only), covariance matrix, per-subject influence functions, and fitted
#'   mean/covariance moments.
#'
#' @references
#' Carrasco, J. L., & Jover, L. (2003). Estimating the generalized concordance
#' correlation coefficient through variance components. \emph{Biometrics},
#' 59(4), 849--858.
#'
#' King, T. S., & Chinchilli, V. M. (2001). A generalized concordance
#' correlation coefficient for continuous and categorical data.
#' \emph{Statistics in Medicine}, 20(14), 2131--2147.
#'
#' King, T. S., Chinchilli, V. M., & Carrasco, J. L. (2007). A repeated
#' measures concordance correlation coefficient. \emph{Statistics in Medicine},
#' 26(16), 3095--3113.
#'
#' Lin, L. I.-K. (1989). A concordance correlation coefficient to evaluate
#' reproducibility. \emph{Biometrics}, 45(1), 255--268.
#'
#' @noRd
.kappa_vector_quadratic <- function(x,
                                   method = c("pairwise", "nt_fiml", "nt_fiml_cs"),
                                   W = NULL,
                                   em_options = list(),
                                   ranges = NULL) {
  method <- match.arg(method)
  prep <- .kvq_flatten(x)
  W <- .kvq_validate_W(W, prep$features)
  widths <- .kvq_validate_ranges(ranges, prep$features)
  fit <- if (method == "pairwise") {
    if (length(widths)) {
      warning("'ranges' is ignored for method = \"pairwise\" ",
              "(it only fixes the Brennan-Prediger baseline, which ",
              "\"pairwise\" does not report).", call. = FALSE)
    }
    .kvq_pairwise(prep$X, prep$R, prep$features, W)
  } else {
    # the estimator name selects the covariance structure, mirroring the scalar
    # surface: nt_fiml is saturated, nt_fiml_cs the block compound symmetry.
    em_options$covariance <- if (method == "nt_fiml_cs") {
      "variance_components"
    } else {
      "saturated"
    }
    .kvq_fiml(prep$X, prep$R, prep$features, W, em_options = em_options,
              widths = widths)
  }

  .new_estimate(fit, names(fit$estimates),
                method = paste0("vector-quadratic-", method),
                weight = "quadratic", W = W,
                moments = list(mu = fit$mu, Sigma = fit$Sigma,
                               summaries = fit$summaries,
                               iterations = fit$iterations,
                               converged = fit$converged))
}
