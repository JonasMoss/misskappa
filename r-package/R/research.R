# Research-only entry points.
#
# Nothing exported reaches the functions in this file. They are comparators
# (available-case, Gwet), complete-data special cases, alternative estimators
# under study (compound-symmetry moments, closed-form quadratic), and
# utilities (multi-weight FIML, joint influence covariance) kept for the
# companion papers' simulations and for parity checks. They are tested but
# unsupported: unexported, named `research_*()`, reached with
# `misskappa:::research_*()`, and free to change or disappear in any release.
# The C++ counterparts live in `misskappa::research` (research.hpp).
#
# Plumbing that the public functions call (dot-prefixed, e.g.
# `.kappa_raw_fit()`) lives next to the public code, not here.

#' Categorical raw-rating kappa with any reweighting (research)
#'
#' @description
#' The categorical estimator behind [agreement()] with the full set of C++
#' reweightings: `"ipw"` and `"fiml"` (as in [agreement()]) plus the
#' available-case (`"available"`) and Gwet (`"gwet"`) comparators. `x` holds
#' integer category codes; see [agreement()] for `categories` and `values`.
#'
#' @return A `misskappa_estimate` with Conger, Fleiss, and Brennan-Prediger.
#' @noRd
research_kappa_raw <- function(x,
                               method = c("available", "ipw", "fiml", "gwet"),
                               weight = "identity",
                               values = NULL,
                               em_options = list(),
                               categories = NULL) {
  method <- match.arg(method)
  .kappa_raw_fit(x, method = method, weight = weight, values = values,
                 em_options = em_options, categories = categories)
}

#' G-wise kappa including the complete-data estimator (research)
#'
#' @description
#' [agreement()] with `g > 2` uses the IPW and categorical FIML g-wise
#' estimators. This also exposes the complete-data estimator
#' (`method = "complete"`), which IPW reduces to when nothing is missing, and
#' the raw distance names (`"nominal"`, `"absolute"`, `"quadratic"`,
#' `"hubert"`).
#'
#' @return A `misskappa_estimate` with Conger and Fleiss coefficients.
#' @noRd
research_kappa_gwise <- function(x,
                                 distance = c("nominal", "absolute",
                                              "quadratic", "hubert"),
                                 method = c("complete", "ipw", "fiml"),
                                 g = NULL,
                                 em_options = list(),
                                 max_chance_tuples = 5000000L) {
  .kappa_gwise(x, distance = match.arg(distance), method = match.arg(method),
               g = g, em_options = em_options,
               max_chance_tuples = max_chance_tuples)
}

#' Multi-weight saturated categorical FIML kappa (research)
#'
#' @description
#' Fits the saturated multinomial FIML model once, then maps the fitted
#' distribution to several agreement-weight functionals. This is intended for
#' simulations and validation code that need nominal, linear, and quadratic
#' Cat-FIML coefficients for the same data matrix. The public [agreement()] API
#' deliberately keeps one `weight` per call.
#'
#' This fixed-rater functional requires every rater to be observed and every
#' rater pair to be jointly observed by at least one subject. Otherwise the
#' requested saturated pairwise agreement functional is not identified.
#'
#' @param x A subjects-by-raters matrix of integer category codes. `NA`s
#'   indicate missing entries.
#' @param weights Character vector containing any of `"identity"`,
#'   `"nominal"`, `"linear"`, or `"quadratic"`. `"nominal"` is an alias for
#'   `"identity"`.
#' @param values Optional numeric category scores for metric weights.
#' @param em_options Named list of categorical EM options.
#'
#' @return A named list of `misskappa_estimate` objects.
#'
#' @noRd
research_kappa_fiml_multi <- function(x,
                                      weights = c("identity", "linear", "quadratic"),
                                      values = NULL,
                                      em_options = list()) {
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a matrix or data frame.", call. = FALSE)
  }
  x_mat <- as.matrix(x)
  if (!is.numeric(x_mat)) stop("'x' must be numeric.", call. = FALSE)
  storage.mode(x_mat) <- "integer"

  choices <- c("identity", "nominal", "linear", "quadratic")
  weights <- match.arg(weights, choices, several.ok = TRUE)
  weights <- unique(ifelse(weights == "nominal", "identity", weights))

  raw <- rcpp_kappa_fiml_multi(
    x = x_mat,
    weight_types = weights,
    values = values,
    em_options = em_options
  )

  out <- lapply(seq_along(raw), function(i) {
    .new_estimate(raw[[i]], c("Conger", "Fleiss", "Brennan-Prediger"),
                  method = "fiml", weight = weights[[i]])
  })
  names(out) <- weights
  out
}


#' Closed-form quadratic-loss kappa estimator (research)
#'
#' @description
#' Closed-form moment-based estimator for the quadratically-weighted
#' Conger / Fleiss / Brennan-Prediger family. Treats categorical ratings
#' as numeric scores and uses per-rater means and covariances.
#'
#' Use this when the quadratic weighting is the right loss and a parametric
#' moment-based estimate is preferred over the U-statistic version
#' (`agreement(x, estimator = "pairwise")`). The two
#' agree under the usual asymptotics. The closed form can be faster on
#' large `n`.
#'
#' With missing ratings, every rater must be observed at least once and every
#' rater pair must be jointly observed by at least one subject.
#'
#' @param x A subjects-by-raters numeric matrix of category scores. `NA`
#'   marks missing entries.
#' @param values Length-C numeric vector of category scores. The quadratic
#'   loss is `(values[i] - values[j])^2 / (max - min)^2`.
#'   The covariance is the empirical covariance of the reduced quadratic
#'   moment summaries from their row-wise estimating equations.
#'
#' @return A `misskappa_estimate` object with `Conger`, `Fleiss`,
#'   `Brennan-Prediger` coefficients, the 3x3 vcov, and per-subject
#'   influence-function rows.
#'
#' @noRd
research_kappa_quadratic <- function(x, values) {
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a matrix or data frame.", call. = FALSE)
  }
  x_mat <- as.matrix(x)
  storage.mode(x_mat) <- "double"
  if (!is.numeric(values)) stop("'values' must be numeric.", call. = FALSE)

  out <- rcpp_kappa_quadratic(
    x = x_mat,
    values = as.numeric(values)
  )
  .new_estimate(out, c("Conger", "Fleiss", "Brennan-Prediger"),
                method = "quadratic", weight = "quadratic")
}


#' Continuous-rating agreement coefficients, all methods (research)
#'
#' @description
#' Estimates Conger / Fleiss weighted agreement coefficients for raw
#' real-valued ratings under missing data. The MCAR moment routes
#' (`"available"`, `"ipw"`, `"gwet"`) cover any continuous loss kernel. The
#' normal-FIML route (`"fiml"`, quadratic kernel only) fits a multivariate-
#' normal covariance by EM and is valid under ignorable (MCAR or MAR)
#' missingness. The default covariance model is saturated. `em_options =
#' list(covariance = "variance_components")` fits the compound-symmetry
#' (fixed-rater) covariance model. Missing entries are encoded as `NA` (or any
#' non-finite value).
#'
#' For the saturated missing-data routes here, every rater must be observed at
#' least once and every rater pair must be jointly observed by at least one
#' subject. The compound-symmetry normal-FIML covariance model only requires
#' a connected rater co-observation graph.
#'
#' @param x A subjects-by-raters numeric matrix.
#' @param method One of `"available"`, `"ipw"`, `"gwet"`, or `"fiml"`.
#'   `"fiml"` requires `weight = "quadratic"` and dispatches to
#'   [.kappa_quadratic_fiml()].
#' @param weight Continuous loss kernel: `"identity"` (binary 0/1),
#'   `"linear"`, `"quadratic"`, `"radical"`, or `"ratio"`. The data range
#'   `[min(x), max(x)]` is used to parameterise the kernel. Only
#'   `"quadratic"` is available with `method = "fiml"`.
#' @param values Optional finite numeric score vector used only for
#'   `method = "fiml"` to define Brennan-Prediger's uniform chance baseline.
#'   If `NULL`, the continuous uniform baseline on the observed score range is
#'   used.
#' @param em_options Used only when `method = "fiml"`: named list with `tol`,
#'   `max_iter`, `covariance`, and `mean`. `covariance = "saturated"` fits the
#'   unrestricted covariance (default). `"variance_components"` fits
#'   the compound-symmetry (fixed-rater) covariance model
#'   \eqn{\Sigma = \psi J + e I}. `mean = "free"` estimates one mean per
#'   rater/item (default). `mean = "equal"` estimates a common mean and
#'   requires `covariance = "variance_components"`.
#'
#' @return A `misskappa_estimate` object carrying named coefficients
#'   (`Conger`, `Fleiss`, `Brennan-Prediger`) and the 3x3 asymptotic
#'   covariance matrix for `method = "fiml"`. Other continuous methods carry
#'   `Conger` and `Fleiss`.
#'
#' @noRd
research_kappa_continuous <- function(x,
                                      method = c("available", "ipw", "gwet", "fiml"),
                                      weight = c("quadratic", "linear",
                                                 "identity", "unweighted",
                                                 "radical", "ratio"),
                                      values = NULL,
                                      em_options = list()) {
  method <- match.arg(method)
  weight <- match.arg(weight)

  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a matrix or data frame.", call. = FALSE)
  }
  x_mat <- as.matrix(x)
  if (!is.numeric(x_mat)) stop("'x' must be numeric.", call. = FALSE)
  storage.mode(x_mat) <- "double"

  if (method == "fiml") {
    if (weight != "quadratic") {
      stop("method = \"fiml\" is only defined for the quadratic kernel; the ",
           "normal-FIML coefficients are functions of the mean and covariance ",
           "only under quadratic weighting.", call. = FALSE)
    }
    return(.kappa_quadratic_fiml(x_mat, values = values, em_options = em_options))
  }

  .kappa_continuous_fit(x_mat, method, weight)
}

# Free-mean compound-symmetry (CS) pairwise moment estimator for the quadratic
# kappas. This is the distribution-free MCAR moment companion to the
# compound-symmetry normal FIML (`estimator = "nt_fiml_cs"`): the
# "pairwise + compound-symmetry" cell of the moment-vs-FIML grid.
#
# Model Y_ij = mu_j + b_i + eps_ij, Var(b_i) = psi, Var(eps_ij) = e, so
# Sigma = psi J + e I (compound symmetry) with FREE rater means mu_j. Recipe:
# rater-center by the available column means, then the same pooled
# within-subject pairwise moments as the counts/exchangeable case, plus a
# mean-spread term t3 with a finite-sample bias correction.
#
#   mu_j = available column mean of rater j
#   r_ij = Y_ij - mu_j
#   e    = sum_i sum_{j<k obs} (r_ij - r_ik)^2 / (2 sum_i P_i),  P_i = C(m_i, 2)
#   psi  = sum_i sum_{j<k obs} r_ij r_ik / sum_i P_i
#   t3   = sum_j (mu_j - mubar)^2 - tr(A V),  A = I - 11'/R, V = Cov(mu_hat)
#
# Coefficients (CS-implied t1 = R^2 psi + R e, t2 = R(psi + e)):
#   Conger = psi / (psi + e + t3/(R-1))
#   Fleiss = (R(R-1) psi - t3) / ((R-1)(R(psi+e) + t3))
#   BP     = 1 - (2e + 2 t3/(R-1)) / bp_chance
#
# Inference is the analytic i.i.d.-subjects sandwich. The estimator is a stacked
# M-estimator in theta = (mu_1..mu_R, psi, e) whose expected cross-Jacobian
# d g_{psi,e} / d mu vanishes at the truth (E[r] = 0), so A is block diagonal and
# the per-subject influence is closed form:
#   IF_mu_j = obs_ij (Y_ij - mu_j) / n_j
#   IF_psi  = (sum_{j<k} r_ij r_ik         - psi P_i) / Ptot
#   IF_e    = (sum_{j<k} (r_ij - r_ik)^2/2 - e   P_i) / Ptot
# t3 enters the coefficients with d t3 / d mu_j = 2(mu_j - mubar); the O(1/n)
# bias correction does not affect the first-order influence. bp_chance is a fixed
# function of the rating support (treated as constant, as in the other quadratic
# paths). Validated against a subject-level bootstrap.
#
# Internal since 2.0.0: no longer selectable through agreement(); reachable for
# simulations as misskappa:::research_kappa_pairwise_cs().

# t3 sampling bias = tr(A V), A = I - 11'/R, V = Cov(mu_hat):
#   V_jj = (psi + e)/n_j,   V_jk = psi * n_jk / (n_j n_k)   (n_jk = co-observations)
.cs_pairwise_t3_bias <- function(psi, e, obs) {
  R <- ncol(obs)
  n_j <- colSums(obs)
  n_jk <- crossprod(obs)                        # R x R co-observation counts
  V <- psi * (n_jk / outer(n_j, n_j))
  diag(V) <- (psi + e) / n_j
  A <- diag(R) - matrix(1 / R, R, R)
  sum(diag(A %*% V))
}

# Coefficient gradients d(Conger, Fleiss, Brennan-Prediger) / d(psi, e, t3).
.cs_pairwise_kappa_jac <- function(psi, e, t3, R, bp_chance) {
  Rm1 <- R - 1
  D <- psi + e + t3 / Rm1
  NF <- R * Rm1 * psi - t3
  DF <- Rm1 * (R * (psi + e) + t3)
  rbind(
    Conger = c(psi = (e + t3 / Rm1) / D^2, e = -psi / D^2, t3 = -psi / (D^2 * Rm1)),
    Fleiss = c(psi = R * Rm1 * (DF - NF) / DF^2, e = -R * Rm1 * NF / DF^2,
               t3 = (-DF - NF * Rm1) / DF^2),
    `Brennan-Prediger` = c(psi = 0, e = -2 / bp_chance, t3 = -2 / (Rm1 * bp_chance))
  )
}

# ---- entry point -------------------------------------------------------------

research_kappa_pairwise_cs <- function(x_mat, values = NULL) {
  X <- as.matrix(x_mat); storage.mode(X) <- "double"
  obs <- is.finite(X)
  keep <- rowSums(obs) > 0
  X <- X[keep, , drop = FALSE]; obs <- obs[keep, , drop = FALSE]
  n <- nrow(X); R <- ncol(X)
  if (R < 2L) stop("quadratic kappa requires at least two raters.",
                   call. = FALSE)
  n_j <- colSums(obs)
  if (any(n_j == 0)) {
    stop("every rater must be observed by at least one subject.", call. = FALSE)
  }

  # mean vector and rater-centered residuals
  mu <- colSums(ifelse(obs, X, 0)) / n_j
  Xc <- X - matrix(mu, n, R, byrow = TRUE)

  # pooled within-subject pairwise moments, plus per-subject scores for the IF
  IF_mu <- matrix(0, n, R)
  for (j in seq_len(R)) IF_mu[, j] <- ifelse(obs[, j], Xc[, j], 0) / n_j[j]
  g_psi <- numeric(n); g_e <- numeric(n); P <- numeric(n)
  num_psi <- 0; num_e <- 0
  for (i in seq_len(n)) {
    o <- which(obs[i, ]); m <- length(o)
    if (m < 2L) next
    ri <- Xc[i, o]; s1 <- sum(ri); s2 <- sum(ri^2); Pi <- m * (m - 1) / 2
    P[i] <- Pi
    g_psi[i] <- (s1^2 - s2) / 2          # sum_{j<k} r_ij r_ik (psi P_i subtracted below)
    g_e[i]   <- (m * s2 - s1^2) / 2      # sum_{j<k} (r_ij - r_ik)^2 / 2
    num_psi <- num_psi + g_psi[i]
    num_e   <- num_e + g_e[i]
  }
  Ptot <- sum(P)
  if (Ptot <= 0) {
    stop("need at least one subject with two or more observed ratings.",
         call. = FALSE)
  }
  psi <- num_psi / Ptot
  e   <- num_e / Ptot

  # mean-spread term with finite-sample bias correction
  mubar <- mean(mu)
  t3 <- max(sum((mu - mubar)^2) - .cs_pairwise_t3_bias(psi, e, obs), 0)

  bp_chance <- .quadratic_bp_chance(X, values)
  estimates <- c(
    Conger = psi / (psi + e + t3 / (R - 1)),
    Fleiss = (R * (R - 1) * psi - t3) / ((R - 1) * (R * (psi + e) + t3)),
    `Brennan-Prediger` = 1 - (2 * e + 2 * t3 / (R - 1)) / bp_chance
  )

  # per-subject influence rows (theta_hat - theta ~ sum_i IF_i)
  IF_psi <- (g_psi - psi * P) / Ptot
  IF_e   <- (g_e - e * P) / Ptot
  IF_t3  <- as.numeric(IF_mu %*% (2 * (mu - mubar)))
  J <- .cs_pairwise_kappa_jac(psi, e, t3, R, bp_chance)
  IF <- cbind(IF_psi, IF_e, IF_t3) %*% t(J)         # n x 3, columns Conger/Fleiss/BP
  colnames(IF) <- names(estimates)

  psi_mat <- n * IF                                  # package convention: vcov = crossprod(psi)/n^2
  .new_estimate(
    list(estimates = estimates, vcov = crossprod(psi_mat) / n^2, psi = psi_mat),
    names(estimates), method = "pairwise_cs", weight = "quadratic",
    bp_chance = bp_chance,
    moments = list(mu = mu, psi = psi, e = e, t3 = t3, n = n, R = R))
}

#' Joint asymptotic covariance across kappa fits on the same data
#'
#' @description
#' Given two or more `misskappa_estimate` objects fit on the **same
#' subjects in the same order**, returns the joint asymptotic covariance
#' of their stacked coefficient vector. The estimators can differ
#' (available-case, IPW, Gwet), use different weight schemes, or be
#' computed on different rater subsets of the same subjects.
#'
#' The joint matrix is assembled from per-subject influence functions:
#' `V = (1 / n^2) * crossprod(cbind(psi_1, ..., psi_K))`. This is the
#' standard nonparametric IF-based joint asymptotic covariance.
#'
#' All `misskappa_estimate` inputs expose per-subject influence functions, so
#' any combination of estimators, weight schemes, or rater subsets fit on the
#' same subjects can be stacked.
#'
#' Row alignment is the caller's responsibility: the helper only checks
#' that the number of subjects matches across inputs and stacks `psi`
#' column-wise.
#'
#' @param ... Two or more `misskappa_estimate` objects. They may optionally be
#'   named, in which case the names are used as prefixes for the output
#'   `dimnames`.
#'
#' @return A square numeric matrix of dimension `K_total x K_total`, where
#'   `K_total` is the total number of coefficients across all inputs.
#'   Row / column names are `<prefix>.<coefficient>` for each fit.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 500
#' R <- 4
#' x <- matrix(sample.int(5, n * R, replace = TRUE), n, R)
#' x[matrix(runif(n * R) < 0.3, n, R)] <- NA
#' ac  <- research_kappa_raw(x, method = "available")
#' ipw <- agreement(x, estimator = "ipw")
#' V   <- research_joint_vcov(ac = ac, ipw = ipw)
#'
#' # Hausman-style Wald test of kappa_AC = kappa_IPW (Conger).
#' delta <- coef(ac)[["Conger"]] - coef(ipw)[["Conger"]]
#' se <- sqrt(V["ac.Conger", "ac.Conger"]
#'          + V["ipw.Conger", "ipw.Conger"]
#'          - 2 * V["ac.Conger", "ipw.Conger"])
#' z <- delta / se
#' }
#'
#' @noRd
research_joint_vcov <- function(...) {
  fits <- list(...)
  if (length(fits) < 2L)
    stop("research_joint_vcov() requires at least two fits.", call. = FALSE)
  if (!all(vapply(fits, inherits, logical(1), "misskappa_estimate")))
    stop("All inputs must be 'misskappa_estimate' objects.", call. = FALSE)

  psis <- lapply(fits, function(f) {
    m <- f$psi
    if (is.null(m) || prod(dim(m)) == 0L) NULL else m
  })
  null_idx <- vapply(psis, is.null, logical(1))
  if (any(null_idx))
    stop("Inputs ", paste(which(null_idx), collapse = ", "),
         " do not carry per-subject influence functions (fit$psi).",
         call. = FALSE)

  ns <- vapply(psis, nrow, integer(1))
  if (length(unique(ns)) != 1L)
    stop("All fits must have the same number of subjects; got ",
         paste(ns, collapse = ", "), ".", call. = FALSE)
  n <- ns[[1L]]

  psi_stack <- do.call(cbind, psis)
  V <- crossprod(psi_stack) / (as.numeric(n) * as.numeric(n))

  prefixes <- names(fits)
  if (is.null(prefixes) || any(!nzchar(prefixes)))
    prefixes <- paste0("fit", seq_along(fits))
  block_names <- unname(unlist(Map(function(p, mat) paste(p, colnames(mat),
                                                          sep = "."),
                                   prefixes, psis)))
  dimnames(V) <- list(block_names, block_names)
  V
}
