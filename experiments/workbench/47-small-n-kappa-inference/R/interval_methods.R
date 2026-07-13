# Experiment-local interval methods for complete-data quadratic Fleiss kappa.
#
# The subject moments are z_i = (mean_j y_ij, mean_j y_ij^2,
# mean_{j<k}(y_ij-y_ik)^2). If s = E(z), then
#
#   kappa_F = 1 - s_disagreement / {2(s_second - s_mean^2)}.
#
# This gives a three-moment robust GMM profile. The same coefficient is the
# intraclass correlation in a common-mean compound-symmetry normal working
# model, whose nuisance mean and scale can be profiled analytically.

.fleiss_subject_moments <- function(x) {
  x <- as.matrix(x)
  n <- nrow(x)
  R <- ncol(x)
  if (n < 4L || R < 2L || any(!is.finite(x))) {
    stop("The focused interval methods need a complete n-by-R matrix, n >= 4 and R >= 2.",
         call. = FALSE)
  }

  pair_index <- utils::combn(R, 2L)
  disagreement <- rowMeans(vapply(
    seq_len(ncol(pair_index)),
    function(k) (x[, pair_index[1L, k]] - x[, pair_index[2L, k]])^2,
    numeric(n)
  ))
  cbind(mean = rowMeans(x), second = rowMeans(x^2),
        disagreement = disagreement)
}

.fleiss_from_moments <- function(moment) {
  chance <- 2 * (moment[["second"]] - moment[["mean"]]^2)
  if (!is.finite(chance) || chance <= 0) return(NA_real_)
  1 - moment[["disagreement"]] / chance
}

.fleiss_moment_fit <- function(x) {
  z <- .fleiss_subject_moments(x)
  n <- nrow(z)
  moment <- colMeans(z)
  estimate <- .fleiss_from_moments(moment)
  chance <- 2 * (moment[["second"]] - moment[["mean"]]^2)
  a <- moment[["disagreement"]]
  gradient <- c(
    mean = -4 * a * moment[["mean"]] / chance^2,
    second = 2 * a / chance^2,
    disagreement = -1 / chance
  )
  centered <- sweep(z, 2L, moment, "-")
  gamma <- crossprod(centered) / n
  influence <- as.numeric(centered %*% gradient)
  se <- sqrt(sum(influence^2)) / n
  list(z = z, moment = moment, estimate = estimate, gradient = gradient,
       gamma = gamma, influence = influence, se = se)
}

.sym_inverse <- function(x, rcond = 1e-10) {
  x <- (x + t(x)) / 2
  eig <- eigen(x, symmetric = TRUE)
  cutoff <- rcond * max(abs(eig$values), 1)
  keep <- eig$values > cutoff
  if (!any(keep)) stop("Moment covariance has no positive eigenvalues.", call. = FALSE)
  inv <- eig$vectors[, keep, drop = FALSE] %*%
    (diag(1 / eig$values[keep], nrow = sum(keep)) %*%
       t(eig$vectors[, keep, drop = FALSE]))
  list(inverse = inv, rank = sum(keep), values = eig$values)
}

.gmm_profile_factory <- function(moment_fit) {
  n <- nrow(moment_fit$z)
  inv <- .sym_inverse(moment_fit$gamma)
  W <- inv$inverse
  observed <- moment_fit$moment
  start <- c(observed[["mean"]],
             log(max(observed[["second"]] - observed[["mean"]]^2, 1e-8)))

  statistic <- function(kappa0) {
    if (!is.finite(kappa0) || kappa0 >= 1) return(Inf)
    objective <- function(par) {
      mu <- par[[1L]]
      variance <- exp(par[[2L]])
      candidate <- c(
        mean = mu,
        second = mu^2 + variance,
        disagreement = 2 * (1 - kappa0) * variance
      )
      delta <- observed - candidate
      value <- as.numeric(n * crossprod(delta, W %*% delta))
      if (is.finite(value)) value else .Machine$double.xmax^0.25
    }
    opt <- stats::optim(start, objective, method = "BFGS",
                        control = list(reltol = 1e-11, maxit = 200L))
    if (opt$convergence != 0L || !is.finite(opt$value)) return(NA_real_)
    max(0, opt$value)
  }

  list(statistic = statistic, rank = inv$rank,
       min_eigenvalue = min(inv$values), max_eigenvalue = max(inv$values))
}

# Profile state for the common-mean, compound-symmetry normal working model.
# Parameterization: Sigma = tau2 * {(1-kappa) I + kappa J}.
.cs_profile_state <- function(x, kappa0) {
  x <- as.matrix(x)
  n <- nrow(x)
  R <- ncol(x)
  lower <- -1 / (R - 1)
  if (!is.finite(kappa0) || kappa0 <= lower || kappa0 >= 1) return(NULL)

  mu <- mean(x)
  residual <- sweep(x, 2L, mu, "-")
  row_sum <- rowSums(residual)
  ss <- rowSums(residual^2)
  between <- row_sum^2 / R
  within <- pmax(ss - between, 0)
  a_within <- 1 - kappa0
  a_between <- 1 + (R - 1) * kappa0
  quadratic <- within / a_within + between / a_between
  tau2 <- mean(quadratic) / R
  if (!is.finite(tau2) || tau2 <= 0) return(NULL)

  logdet_A <- (R - 1) * log(a_within) + log(a_between)
  loglik_i <- -0.5 * (R * log(2 * pi) + R * log(tau2) +
                        logdet_A + quadratic / tau2)

  score_scale <- -0.5 * (R - quadratic / tau2)
  d_logdet <- -(R - 1) / a_within + (R - 1) / a_between
  d_quadratic <- within / a_within^2 -
    (R - 1) * between / a_between^2
  score_kappa <- -0.5 * (d_logdet + d_quadratic / tau2)

  A_scale_scale <- R / 2
  A_scale_kappa <- d_logdet / 2
  A_kappa_kappa <- 0.5 * ((R - 1) / a_within^2 +
                           (R - 1)^2 / a_between^2)
  efficient_score <- score_kappa -
    (A_scale_kappa / A_scale_scale) * score_scale
  A_eff <- A_kappa_kappa - A_scale_kappa^2 / A_scale_scale
  B_eff <- mean((efficient_score - mean(efficient_score))^2)

  list(loglik = sum(loglik_i), loglik_i = loglik_i, mu = mu, tau2 = tau2,
       efficient_score = efficient_score, A_eff = A_eff, B_eff = B_eff,
       scale = B_eff / A_eff)
}

.cs_test_factory <- function(x, center) {
  n <- nrow(x)
  center_state <- .cs_profile_state(x, center)
  if (is.null(center_state)) stop("Normal-profile center is outside the CS parameter space.",
                                  call. = FALSE)

  state_at <- function(kappa0) .cs_profile_state(x, kappa0)

  profile_lr <- function(kappa0) {
    state <- state_at(kappa0)
    if (is.null(state)) return(Inf)
    max(0, 2 * (center_state$loglik - state$loglik))
  }
  profile_lr_robust <- function(kappa0) {
    state <- state_at(kappa0)
    if (is.null(state) || !is.finite(state$scale) || state$scale <= 0) return(NA_real_)
    profile_lr(kappa0) / state$scale
  }
  score_normal <- function(kappa0) {
    state <- state_at(kappa0)
    if (is.null(state) || !is.finite(state$A_eff) || state$A_eff <= 0) return(NA_real_)
    n * mean(state$efficient_score)^2 / state$A_eff
  }
  score_robust <- function(kappa0) {
    state <- state_at(kappa0)
    if (is.null(state) || !is.finite(state$B_eff) || state$B_eff <= 0) return(NA_real_)
    n * mean(state$efficient_score)^2 / state$B_eff
  }
  score_robust_t <- function(kappa0) {
    state <- state_at(kappa0)
    if (is.null(state)) return(NA_real_)
    variance <- stats::var(state$efficient_score)
    if (!is.finite(variance) || variance <= 0) return(NA_real_)
    n * mean(state$efficient_score)^2 / variance
  }

  list(profile_lr = profile_lr, profile_lr_robust = profile_lr_robust,
       score_normal = score_normal, score_robust = score_robust,
       score_robust_t = score_robust_t,
       center_scale = center_state$scale)
}

.invert_scalar_test <- function(center, statistic, critical, lower, upper,
                                scan_points = 18L) {
  if (!is.finite(center) || center <= lower || center >= upper) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  objective <- function(x) {
    value <- statistic(x)
    if (!is.finite(value)) return(.Machine$double.xmax^0.25)
    value - critical
  }

  find_limit <- function(boundary, side) {
    # Concentrate points near the estimate, then extend to the parameter bound.
    fractions <- c(seq(0, 0.25, length.out = scan_points / 2L + 1L)[-1L],
                   seq(0.30, 1, length.out = scan_points / 2L))
    grid <- center + fractions * (boundary - center)
    previous <- center
    f_previous <- objective(previous)
    for (point in grid) {
      f_point <- objective(point)
      if (is.finite(f_point) && f_point >= 0 && f_previous <= 0) {
        interval <- sort(c(previous, point))
        root <- tryCatch(
          stats::uniroot(objective, interval = interval, tol = 1e-7)$root,
          error = function(e) NA_real_
        )
        return(root)
      }
      previous <- point
      f_previous <- f_point
    }
    boundary
  }

  c(lower = find_limit(lower, "lower"), upper = find_limit(upper, "upper"))
}

.focused_intervals <- function(x, package_estimate, package_se, level = 0.95) {
  n <- nrow(x)
  R <- ncol(x)
  alpha <- 1 - level
  zcrit <- stats::qnorm(1 - alpha / 2)
  chi_critical <- stats::qchisq(level, df = 1)
  t_critical <- stats::qt(1 - alpha / 2, df = n - 3L)^2
  lower_bound <- -1 / (R - 1) + 1e-6
  upper_bound <- 1 - 1e-6

  wald <- c(lower = package_estimate - zcrit * package_se,
            upper = package_estimate + zcrit * package_se)
  if (abs(package_estimate) < 1 && is.finite(package_se)) {
    fisher_se <- package_se / (1 - package_estimate^2)
    fisher <- tanh(atanh(package_estimate) + c(lower = -1, upper = 1) *
                     zcrit * fisher_se)
  } else {
    fisher <- c(lower = NA_real_, upper = NA_real_)
  }

  moment_fit <- .fleiss_moment_fit(x)
  center <- moment_fit$estimate
  empty <- c(lower = NA_real_, upper = NA_real_)

  # Exact observed agreement puts the normal working MLE at singular Sigma and
  # removes the disagreement direction from the empirical moment covariance.
  # Keep the valid Wald result, but mark profile/score methods unavailable.
  if (!is.finite(center) || center <= lower_bound || center >= upper_bound) {
    intervals <- rbind(
      wald_z = wald, wald_fisher = fisher,
      profile_lr = empty, profile_lr_robust = empty,
      score_normal = empty, score_robust = empty, score_robust_t = empty,
      moment_profile = empty, moment_profile_f = empty
    )
    return(list(
      intervals = intervals, moment_estimate = center, moment_se = moment_fit$se,
      package_delta_estimate = package_estimate - center,
      package_delta_se = package_se - moment_fit$se,
      cs_center_scale = NA_real_, gmm_rank = NA_integer_,
      gmm_min_eigenvalue = NA_real_, gmm_max_eigenvalue = NA_real_
    ))
  }

  cs <- .cs_test_factory(x, center)
  gmm <- .gmm_profile_factory(moment_fit)

  tests <- list(
    profile_lr = list(fn = cs$profile_lr, critical = chi_critical),
    profile_lr_robust = list(fn = cs$profile_lr_robust, critical = chi_critical),
    score_normal = list(fn = cs$score_normal, critical = chi_critical),
    score_robust = list(fn = cs$score_robust, critical = chi_critical),
    score_robust_t = list(fn = cs$score_robust_t, critical = t_critical),
    moment_profile = list(fn = gmm$statistic, critical = chi_critical),
    moment_profile_f = list(fn = gmm$statistic, critical = t_critical)
  )
  inverted <- lapply(tests, function(spec) {
    .invert_scalar_test(center, spec$fn, spec$critical,
                        lower_bound, upper_bound)
  })

  intervals <- rbind(wald_z = wald, wald_fisher = fisher,
                     do.call(rbind, inverted))
  list(
    intervals = intervals,
    moment_estimate = moment_fit$estimate,
    moment_se = moment_fit$se,
    package_delta_estimate = package_estimate - moment_fit$estimate,
    package_delta_se = package_se - moment_fit$se,
    cs_center_scale = cs$center_scale,
    gmm_rank = gmm$rank,
    gmm_min_eigenvalue = gmm$min_eigenvalue,
    gmm_max_eigenvalue = gmm$max_eigenvalue
  )
}
