# Quadratic (Conger / Fleiss / Brennan-Prediger) kappa for continuous ratings
# under ignorable missingness, via a normal FIML / EM covariance and a
# delta-method standard error. This is the quadratic-kappa counterpart of
# alpha_continuous(): the two share the normal EM, the analytic
# casewise scores, the observed information, and the sandwich contraction. The
# only coefficient-specific piece is the gradient .kqf_grad().
#
# The EM, score, and information calculations live in the C++ normal-FIML
# backend shared with alpha_continuous().
#
# Unlike alpha, the quadratic kappas depend on the rater means through
# t3 = sum_j (mu_j - mubar)^2, so the gradient carries a non-zero mean block.

.normal_covariance_model <- function(em_options) {
  cov <- em_options$covariance
  if (is.null(cov)) return("saturated")
  if (!is.character(cov) || length(cov) != 1L || is.na(cov)) {
    stop("em_options$covariance must be a single string.", call. = FALSE)
  }
  cov <- gsub("-", "_", tolower(trimws(cov)), fixed = TRUE)
  if (cov %in% c("saturated")) return("saturated")
  if (cov %in% c("variance_components")) return("variance_components")
  stop("unsupported em_options$covariance: ", em_options$covariance,
       call. = FALSE)
}

.quadratic_bp_chance <- function(x, values = NULL) {
  if (is.null(values)) {
    observed <- x[is.finite(x)]
    if (!length(observed)) stop("at least one rating must be observed.", call. = FALSE)
    rng <- range(observed)
    chance <- diff(rng)^2 / 6
  } else {
    if (!is.numeric(values)) stop("'values' must be numeric.", call. = FALSE)
    values <- as.numeric(values)
    if (length(values) < 2L || any(!is.finite(values))) {
      stop("'values' must contain at least two finite scores.", call. = FALSE)
    }
    C <- length(values)
    chance <- (2 / C^2) * (C * sum(values^2) - sum(values)^2)
  }
  if (!is.finite(chance) || chance <= 0) {
    stop("Brennan-Prediger's quadratic chance disagreement is zero.",
         call. = FALSE)
  }
  chance
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

# ---- public estimator -------------------------------------------------------

#' Quadratic kappa under missing data (normal FIML)
#'
#' @description
#' Backend for `kappa(estimator = "nt_fiml")`. Estimates the quadratically
#' weighted Conger, Fleiss, and Brennan-Prediger agreement coefficients for
#' continuous (or numerically scored) ratings with missing entries, via a
#' multivariate-normal covariance fitted by full-information maximum likelihood
#' (the EM algorithm), valid under ignorable (MCAR or MAR) missingness. The
#' default covariance model is saturated; `em_options = list(covariance =
#' "variance_components")` fits the compound-symmetry (fixed-rater) covariance
#' model \eqn{\Sigma = \psi J + e I}. It treats raters as fixed; it is not ten
#' Hove's random-rater MLE-RE / ICC. The coefficients are smooth functions of
#' the fitted mean and covariance, and the standard error is a delta-method
#' contraction of their gradient with the asymptotic covariance of the fitted
#' moments. This is the quadratic-kappa counterpart of [alpha_continuous()] and
#' shares its EM, casewise-score, and sandwich machinery.
#'
#' This backend is normally reached through [kappa()], which checks that every
#' rater is observed and, for the default saturated covariance model, every
#' rater pair is jointly observed before dispatch. Direct callers should
#' enforce the same fixed-rater condition; for `covariance =
#' "variance_components"` this is a connected rater co-observation graph.
#'
#' @param x A subjects-by-raters numeric matrix or data frame; `NA` marks
#'   missing entries. Rows that are entirely missing are dropped.
#' @param values Optional finite numeric score vector defining the uniform
#'   Brennan-Prediger chance baseline. If `NULL`, the continuous uniform
#'   baseline on the observed score range is used.
#' @param em_options Named list tuning the EM fit: `tol` (default `1e-8`),
#'   `max_iter` (default `10000`), and `covariance`. `covariance = "saturated"`
#'   fits the unrestricted normal FIML model; `"variance_components"` fits the
#'   compound-symmetry (fixed-rater) covariance model
#'   \eqn{\Sigma = \psi J + e I} with \eqn{\psi \ge 0}. `fd_h` is accepted for
#'   backward compatibility and ignored because the observed information is
#'   analytic. Pass any subset.
#'
#' @return An object of class `misskappa_estimate` carrying the `Conger`,
#'   `Fleiss`, and `Brennan-Prediger` coefficients and their 3x3 asymptotic
#'   covariance. Additional fields: `moments` (the fitted `mu`, `Sigma`, EM
#'   `iterations`, and `converged` flag) and `psi` (per-subject
#'   influence-function rows). Methods: `print`, `coef`, `vcov`, `confint`,
#'   `as.data.frame`, and `stats::influence`.
#'
#' @keywords internal
kappa_quadratic_fiml <- function(x, values = NULL, em_options = list()) {
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a matrix or data frame.")
  }
  X <- as.matrix(x)
  if (!is.numeric(X)) stop("'x' must be numeric.")
  storage.mode(X) <- "double"
  p <- ncol(X)
  if (p < 2L) stop("quadratic kappa requires at least two raters.")

  R <- is.finite(X)
  keep <- rowSums(R) > 0L
  if (!all(keep)) {
    X <- X[keep, , drop = FALSE]
    R <- R[keep, , drop = FALSE]
  }
  if (any(colSums(R) == 0L)) {
    stop("every rater must be observed for at least one subject.")
  }

  opt <- utils::modifyList(
    list(tol = 1e-8, max_iter = 10000L), em_options
  )
  covariance <- .normal_covariance_model(opt)
  bp_chance <- .quadratic_bp_chance(X, values)

  fit <- .normal_fiml_from_cpp(
    rcpp_kappa_quadratic_fiml(X, values, opt),
    estimate_names = c("Conger", "Fleiss", "Brennan-Prediger"),
    method = if (covariance == "variance_components") {
      "quadratic-fiml-variance-components"
    } else {
      "quadratic-fiml"
    },
    weight = "quadratic"
  )
  fit$covariance <- covariance
  fit$bp_chance <- bp_chance
  fit
}
