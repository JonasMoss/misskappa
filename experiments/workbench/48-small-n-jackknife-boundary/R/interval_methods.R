# Experiment-local small-sample intervals for complete-data quadratic Fleiss
# kappa. The implementation works on one row of subject moments per item so
# the nested bootstrap/jackknife remains cheap and transparent.

.fleiss_subject_moments <- function(x) {
  x <- as.matrix(x)
  n <- nrow(x)
  R <- ncol(x)
  if (n < 4L || R < 2L || any(!is.finite(x))) {
    stop("Expected a complete n-by-R matrix with n >= 4 and R >= 2.",
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

.fleiss_moment_fit <- function(z) {
  z <- as.matrix(z)
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
  influence <- as.numeric(centered %*% gradient)
  se <- sqrt(sum(influence^2)) / n
  list(moment = moment, estimate = estimate, influence = influence, se = se)
}

.leave_one_out_functional <- function(z, fn = .fleiss_from_moments) {
  z <- as.matrix(z)
  n <- nrow(z)
  total <- colSums(z)
  vapply(seq_len(n), function(i) fn((total - z[i, ]) / (n - 1)), numeric(1))
}

.jackknife_se <- function(leave_one_out) {
  n <- length(leave_one_out)
  if (n < 2L || any(!is.finite(leave_one_out))) return(NA_real_)
  sqrt((n - 1) / n * sum((leave_one_out - mean(leave_one_out))^2))
}

.clip_interval <- function(interval, lower_bound, upper_bound = 1) {
  interval <- sort(as.numeric(interval))
  if (length(interval) != 2L || any(!is.finite(interval))) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  c(lower = max(lower_bound, interval[[1L]]),
    upper = min(upper_bound, interval[[2L]]))
}

.fisher_interval <- function(estimate, se, critical, lower_bound) {
  if (!is.finite(estimate) || !is.finite(se) || se < 0 || abs(estimate) >= 1) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  eta <- atanh(estimate)
  eta_se <- se / (1 - estimate^2)
  .clip_interval(tanh(eta + c(-1, 1) * critical * eta_se), lower_bound)
}

.direct_fisher_jackknife <- function(estimate, leave_one_out, critical,
                                     lower_bound, boundary_tolerance) {
  if (!is.finite(estimate) || abs(estimate) >= 1 - boundary_tolerance ||
      any(!is.finite(leave_one_out)) ||
      any(abs(leave_one_out) >= 1 - boundary_tolerance)) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  eta_se <- .jackknife_se(atanh(leave_one_out))
  .clip_interval(
    tanh(atanh(estimate) + c(-1, 1) * critical * eta_se), lower_bound
  )
}

.bca_interval <- function(estimate, boot_estimates, leave_one_out, level,
                          lower_bound) {
  boot_estimates <- boot_estimates[is.finite(boot_estimates)]
  if (length(boot_estimates) < 20L || any(!is.finite(leave_one_out))) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  alpha <- 1 - level
  less <- sum(boot_estimates < estimate)
  ties <- sum(boot_estimates == estimate)
  probability <- (less + 0.5 * ties) / length(boot_estimates)
  probability <- min(max(probability, 0.5 / length(boot_estimates)),
                     1 - 0.5 / length(boot_estimates))
  z0 <- stats::qnorm(probability)

  centered <- mean(leave_one_out) - leave_one_out
  denominator <- 6 * sum(centered^2)^(3 / 2)
  acceleration <- if (is.finite(denominator) && denominator > 0) {
    sum(centered^3) / denominator
  } else {
    0
  }
  z_alpha <- stats::qnorm(c(alpha / 2, 1 - alpha / 2))
  adjusted <- stats::pnorm(
    z0 + (z0 + z_alpha) / (1 - acceleration * (z0 + z_alpha))
  )
  adjusted <- pmin(pmax(adjusted, 1 / (length(boot_estimates) + 1)),
                   length(boot_estimates) / (length(boot_estimates) + 1))
  .clip_interval(
    stats::quantile(boot_estimates, adjusted, names = FALSE, type = 6),
    lower_bound
  )
}

.bootstrap_t_interval <- function(estimate, se, bootstrap_t, level,
                                  lower_bound, min_valid_fraction) {
  valid <- is.finite(bootstrap_t)
  if (!is.finite(se) || se <= 0 ||
      mean(valid) < min_valid_fraction || sum(valid) < 20L) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  alpha <- 1 - level
  quantiles <- stats::quantile(
    bootstrap_t[valid], c(1 - alpha / 2, alpha / 2),
    names = FALSE, type = 6
  )
  .clip_interval(estimate - quantiles * se, lower_bound)
}

.chance_from_moments <- function(moment, score_range) {
  2 * (moment[["second"]] - moment[["mean"]]^2) / score_range^2
}

.boundary_interval <- function(z, score_range, level, lower_bound) {
  n <- nrow(z)
  alpha <- 1 - level
  normalized_disagreement <- z[, "disagreement"] / score_range^2
  event_count <- sum(normalized_disagreement > 0)

  # A = E(D) <= Pr(D > 0) because normalized D is in [0, 1]. The
  # one-sided Clopper-Pearson limit therefore upper-bounds A as well.
  p_upper <- if (event_count >= n) {
    1
  } else {
    stats::qbeta(1 - alpha / 2, event_count + 1, n - event_count)
  }

  chance_fn <- function(moment) .chance_from_moments(moment, score_range)
  chance <- chance_fn(colMeans(z))
  chance_loo <- .leave_one_out_functional(z, chance_fn)
  chance_se <- .jackknife_se(chance_loo)
  chance_lower <- chance - stats::qt(1 - alpha / 2, n - 1) * chance_se

  lower <- if (!is.finite(chance_lower) || chance_lower <= 0) {
    lower_bound
  } else {
    max(lower_bound, 1 - p_upper / chance_lower)
  }
  c(lower = lower, upper = 1)
}

.focused_intervals <- function(x, package_estimate, package_se,
                               bootstrap_reps = 399L, level = 0.95,
                               min_valid_fraction = 0.70,
                               bootstrap_seed = NULL,
                               score_range = max(x) - min(x)) {
  x <- as.matrix(x)
  n <- nrow(x)
  R <- ncol(x)
  lower_bound <- -1 / (R - 1)
  boundary_tolerance <- 1e-12
  z_critical <- stats::qnorm(1 - (1 - level) / 2)
  t_critical <- stats::qt(1 - (1 - level) / 2, n - 1)

  z <- .fleiss_subject_moments(x)
  fit <- .fleiss_moment_fit(z)
  estimate <- fit$estimate
  leave_one_out <- .leave_one_out_functional(z)
  jackknife_se <- .jackknife_se(leave_one_out)
  if (!is.finite(score_range) || score_range <= 0) {
    stop("score_range must be positive.", call. = FALSE)
  }
  fallback <- .boundary_interval(z, score_range, level, lower_bound)

  wald <- .clip_interval(
    package_estimate + c(-1, 1) * z_critical * package_se, lower_bound
  )
  fisher_if <- .fisher_interval(package_estimate, package_se, z_critical,
                                lower_bound)
  jack_t <- .clip_interval(
    estimate + c(-1, 1) * t_critical * jackknife_se, lower_bound
  )
  jack_fisher_delta <- .fisher_interval(
    estimate, jackknife_se, t_critical, lower_bound
  )
  jack_fisher_direct <- .direct_fisher_jackknife(
    estimate, leave_one_out, t_critical, lower_bound, boundary_tolerance
  )
  jack_fisher_direct_hybrid <- jack_fisher_direct
  direct_fallback <- any(!is.finite(jack_fisher_direct_hybrid))
  if (direct_fallback) jack_fisher_direct_hybrid <- fallback

  if (!is.null(bootstrap_seed)) set.seed(bootstrap_seed)
  boot_estimates <- rep(NA_real_, bootstrap_reps)
  bootstrap_t <- rep(NA_real_, bootstrap_reps)
  boot_boundary <- rep(FALSE, bootstrap_reps)
  for (b in seq_len(bootstrap_reps)) {
    z_boot <- z[sample.int(n, n, replace = TRUE), , drop = FALSE]
    boot_fit <- .fleiss_moment_fit(z_boot)
    boot_estimates[[b]] <- boot_fit$estimate
    boot_loo <- .leave_one_out_functional(z_boot)
    boot_se <- .jackknife_se(boot_loo)
    boot_boundary[[b]] <- is.finite(boot_fit$estimate) &&
      boot_fit$estimate >= 1 - boundary_tolerance
    if (is.finite(boot_fit$estimate) && is.finite(boot_se) && boot_se > 0) {
      bootstrap_t[[b]] <- (boot_fit$estimate - estimate) / boot_se
    }
  }

  bca <- .bca_interval(estimate, boot_estimates, leave_one_out, level,
                       lower_bound)
  bca_boundary <- is.finite(estimate) &&
    estimate >= 1 - boundary_tolerance
  bca_hybrid <- if (bca_boundary) fallback else bca

  boot_t <- .bootstrap_t_interval(
    estimate, jackknife_se, bootstrap_t, level, lower_bound,
    min_valid_fraction
  )
  boot_t_fallback <- any(!is.finite(boot_t))
  boot_t_hybrid <- if (boot_t_fallback) fallback else boot_t

  intervals <- rbind(
    wald_if_z = wald,
    wald_if_fisher = fisher_if,
    jack_t = jack_t,
    jack_fisher_delta_t = jack_fisher_delta,
    jack_fisher_direct_t = jack_fisher_direct,
    jack_fisher_direct_boundary = jack_fisher_direct_hybrid,
    bca = bca,
    bca_boundary = bca_hybrid,
    bootstrap_t_jack = boot_t,
    bootstrap_t_jack_boundary = boot_t_hybrid
  )
  fallback_used <- c(
    wald_if_z = FALSE,
    wald_if_fisher = FALSE,
    jack_t = FALSE,
    jack_fisher_delta_t = FALSE,
    jack_fisher_direct_t = FALSE,
    jack_fisher_direct_boundary = direct_fallback,
    bca = FALSE,
    bca_boundary = bca_boundary,
    bootstrap_t_jack = FALSE,
    bootstrap_t_jack_boundary = boot_t_fallback
  )

  list(
    intervals = intervals,
    fallback_used = fallback_used[rownames(intervals)],
    moment_estimate = estimate,
    moment_se = fit$se,
    jackknife_se = jackknife_se,
    disagreement_events = sum(z[, "disagreement"] > 0),
    exact_agreement = all(z[, "disagreement"] == 0),
    leave_one_out_boundary = any(leave_one_out >= 1 - boundary_tolerance),
    valid_bootstrap_t_fraction = mean(is.finite(bootstrap_t)),
    bootstrap_boundary_fraction = mean(boot_boundary),
    package_delta_estimate = package_estimate - estimate,
    package_delta_se = package_se - fit$se
  )
}
