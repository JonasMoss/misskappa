#!/usr/bin/env Rscript
#
# Experiment 08: literature bootstrap intervals for quadratic kappa.
#
# Complete categorical Perreault-Leigh simulations for quadratically weighted
# kappa. The interval menu is restricted to named literature-style methods:
# Moss (2024) t-Wald intervals on identity/arcsine/Fisher scales, and standard
# nonparametric bootstrap percentile/basic/BCa/bootstrap-t intervals.

usage <- function(status = 0L) {
  cat(paste0(
    "Usage: Rscript run_experiment.R [options]\n\n",
    "Options:\n",
    "  --help              Show this help and exit.\n",
    "  --smoke             Cheap mechanical check.\n",
    "  --n-grid CSV        Sample sizes. Default: 10,20,40.\n",
    "  --r-grid CSV        Rater counts. Default: 2,5.\n",
    "  --reps N            Monte Carlo replicates per design cell. Default: 200.\n",
    "  --boot-reps B       Bootstrap resamples per replicate. Default: 199.\n",
    "  --target-grid CSV   True quadratic kappas. Default: 0.8,0.9.\n",
    "  --seed-base N       Base seed for deterministic runs.\n",
    "  --out-dir PATH      Output directory. Default: script-local results/.\n",
    "  --progress          Print one line per design cell.\n\n",
    "Examples:\n",
    "  Rscript run_experiment.R --smoke\n",
    "  Rscript run_experiment.R --reps 100 --boot-reps 199 --target-grid 0.8,0.9 --n-grid 10,20 --r-grid 2,5 --progress\n"
  ))
  quit(save = "no", status = status)
}

parse_int_csv <- function(x, min_value, arg_name) {
  out <- as.integer(strsplit(x, ",", fixed = TRUE)[[1L]])
  if (any(is.na(out)) || any(out < min_value)) {
    stop(arg_name, " must contain integers >= ", min_value, ".", call. = FALSE)
  }
  unique(out)
}

parse_num_csv <- function(x, arg_name) {
  out <- as.numeric(strsplit(x, ",", fixed = TRUE)[[1L]])
  if (any(is.na(out)) || any(out <= 0) || any(out >= 1)) {
    stop(arg_name, " must contain numbers in (0, 1).", call. = FALSE)
  }
  unique(out)
}

parse_args <- function(argv) {
  opts <- list(
    smoke = FALSE,
    n_grid = c(10L, 20L, 40L),
    r_grid = c(2L, 5L),
    reps = 200L,
    boot_reps = 199L,
    target_grid = c(0.8, 0.9),
    seed_base = 808080L,
    out_dir = NULL,
    progress = FALSE
  )
  i <- 1L
  while (i <= length(argv)) {
    arg <- argv[[i]]
    if (arg == "--help" || arg == "-h") usage(0L)
    if (arg == "--smoke") {
      opts$smoke <- TRUE
      i <- i + 1L
      next
    }
    if (arg == "--progress") {
      opts$progress <- TRUE
      i <- i + 1L
      next
    }
    needs_value <- c(
      "--n-grid", "--r-grid", "--reps", "--boot-reps", "--target-grid",
      "--seed-base", "--out-dir"
    )
    if (arg %in% needs_value) {
      if (i == length(argv)) stop(arg, " needs a value.", call. = FALSE)
      val <- argv[[i + 1L]]
      if (arg == "--n-grid") opts$n_grid <- parse_int_csv(val, 2L, "--n-grid")
      if (arg == "--r-grid") opts$r_grid <- parse_int_csv(val, 2L, "--r-grid")
      if (arg == "--reps") opts$reps <- as.integer(val)
      if (arg == "--boot-reps") opts$boot_reps <- as.integer(val)
      if (arg == "--target-grid") opts$target_grid <- parse_num_csv(val, "--target-grid")
      if (arg == "--seed-base") opts$seed_base <- as.integer(val)
      if (arg == "--out-dir") opts$out_dir <- val
      i <- i + 2L
      next
    }
    stop("Unknown argument: ", arg, call. = FALSE)
  }
  if (opts$smoke) {
    opts$n_grid <- c(10L, 20L)
    opts$r_grid <- 2L
    opts$reps <- 12L
    opts$boot_reps <- 39L
    opts$target_grid <- c(0.8, 0.9)
  }
  if (is.na(opts$reps) || opts$reps < 1L) stop("--reps must be >= 1.", call. = FALSE)
  if (is.na(opts$boot_reps) || opts$boot_reps < 19L) {
    stop("--boot-reps must be >= 19.", call. = FALSE)
  }
  if (is.na(opts$seed_base)) stop("--seed-base must be an integer.", call. = FALSE)
  opts
}

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE))
} else {
  getwd()
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))
if (is.null(opts$out_dir)) opts$out_dir <- file.path(script_dir, "results")

category_count <- 5L
values <- -2:2
alpha <- 0.05
coef_names <- c("Conger", "Fleiss")
analytic_methods <- c("moss_basic", "moss_arcsine", "moss_fisher")
bootstrap_methods <- c(
  "boot_percentile", "boot_basic", "boot_bca", "boot_t",
  "boot_t_arcsine", "boot_t_fisher"
)
methods <- c(analytic_methods, bootstrap_methods)

clip_unit <- function(x) pmin(pmax(x, -1 + 1e-10), 1 - 1e-10)

simulate_perreault_leigh <- function(n, R, knowledge_prob) {
  truth <- sample.int(category_count, n, replace = TRUE)
  x <- matrix(NA_integer_, nrow = n, ncol = R)
  for (j in seq_len(R)) {
    knows <- stats::runif(n) < knowledge_prob
    guesses <- sample.int(category_count, n, replace = TRUE)
    x[, j] <- ifelse(knows, truth, guesses)
  }
  matrix(values[x], nrow = n, ncol = R)
}

quadratic_loss <- function(values) {
  range_sq <- (max(values) - min(values))^2
  outer(values, values, function(a, b) (a - b)^2 / range_sq)
}

quadratic_if_complete <- function(x, values) {
  n <- nrow(x)
  R <- ncol(x)
  C <- length(values)
  L <- quadratic_loss(values)
  idx <- matrix(match(as.integer(x), values), nrow = n, ncol = R)
  if (any(is.na(idx))) stop("Rating outside value support.", call. = FALSE)

  subject_counts <- matrix(0, nrow = n, ncol = C)
  rater_counts <- matrix(0, nrow = R, ncol = C)
  h_dN <- numeric(n)
  for (i in seq_len(n)) {
    for (j in seq_len(R)) {
      a <- idx[i, j]
      subject_counts[i, a] <- subject_counts[i, a] + 1
      rater_counts[j, a] <- rater_counts[j, a] + 1
    }
    if (R > 1L) {
      for (j in seq_len(R - 1L)) {
        for (k in (j + 1L):R) {
          h_dN[i] <- h_dN[i] + L[idx[i, j], idx[i, k]]
        }
      }
    }
  }

  pair_count <- choose(R, 2)
  psi_dN <- mean(h_dN)
  d_hat <- psi_dN / pair_count
  phi_d <- (h_dN - psi_dN) / pair_count

  category_totals <- colSums(subject_counts)
  subject_L <- subject_counts %*% L
  row_FN <- as.numeric(subject_L %*% category_totals)
  col_FN <- as.numeric((subject_counts %*% t(L)) %*% category_totals)
  psi_FN <- sum(row_FN) / (n * n)
  d_fleiss <- psi_FN / (R * R)
  phi_f_chance <- ((row_FN / n - psi_FN) + (col_FN / n - psi_FN)) / (R * R)

  suffix_counts <- matrix(0, nrow = R, ncol = C)
  prefix_counts <- matrix(0, nrow = R, ncol = C)
  running <- numeric(C)
  for (j in R:1L) {
    suffix_counts[j, ] <- running
    running <- running + rater_counts[j, ]
  }
  running <- numeric(C)
  for (j in seq_len(R)) {
    prefix_counts[j, ] <- running
    running <- running + rater_counts[j, ]
  }

  row_CN <- numeric(n)
  col_CN <- numeric(n)
  for (i in seq_len(n)) {
    for (j in seq_len(R)) {
      a <- idx[i, j]
      row_CN[i] <- row_CN[i] + sum(L[a, ] * suffix_counts[j, ])
      col_CN[i] <- col_CN[i] + sum(prefix_counts[j, ] * L[, a])
    }
  }
  psi_CN <- sum(row_CN) / (n * n)
  d_conger <- psi_CN / pair_count
  phi_c_chance <- ((row_CN / n - psi_CN) + (col_CN / n - psi_CN)) / pair_count

  estimates <- c(
    Conger = 1 - d_hat / d_conger,
    Fleiss = 1 - d_hat / d_fleiss
  )
  phi <- cbind(
    Conger = -phi_d / d_conger + d_hat * phi_c_chance / (d_conger^2),
    Fleiss = -phi_d / d_fleiss + d_hat * phi_f_chance / (d_fleiss^2)
  )
  phi <- sweep(phi, 2L, colMeans(phi), "-")
  se <- sqrt(colMeans(phi^2) / n)
  list(estimates = estimates, se = se)
}

analytic_interval <- function(theta, se, n, method) {
  c_t <- stats::qt(1 - alpha / 2, df = n - 1)
  if (!is.finite(theta) || !is.finite(se) || se < 0) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  if (method == "moss_basic") {
    return(c(lower = theta - c_t * se, upper = theta + c_t * se))
  }

  theta_c <- clip_unit(theta)
  if (method == "moss_arcsine") {
    se_t <- se / sqrt(1 - theta_c^2)
    return(sort(sin(asin(theta_c) + c(-1, 1) * c_t * se_t)))
  }
  if (method == "moss_fisher") {
    se_t <- se / (1 - theta_c^2)
    return(sort(tanh(atanh(theta_c) + c(-1, 1) * c_t * se_t)))
  }
  stop("Unknown analytic method: ", method, call. = FALSE)
}

quantile_interval <- function(x, probs) {
  x <- x[is.finite(x)]
  if (length(x) < 10L) return(c(lower = NA_real_, upper = NA_real_))
  as.numeric(stats::quantile(x, probs = probs, names = FALSE, type = 7))
}

bca_interval <- function(theta, boot_theta, jack_theta) {
  boot_theta <- boot_theta[is.finite(boot_theta)]
  jack_theta <- jack_theta[is.finite(jack_theta)]
  B <- length(boot_theta)
  if (!is.finite(theta) || B < 20L || length(jack_theta) < 3L) {
    return(c(lower = NA_real_, upper = NA_real_))
  }

  prop_less <- mean(boot_theta < theta)
  prop_less <- pmin(pmax(prop_less, 1 / (2 * B)), 1 - 1 / (2 * B))
  z0 <- stats::qnorm(prop_less)

  jack_mean <- mean(jack_theta)
  diffs <- jack_mean - jack_theta
  denom <- 6 * sum(diffs^2)^(3 / 2)
  accel <- if (is.finite(denom) && denom > 0) sum(diffs^3) / denom else 0

  z_alpha <- stats::qnorm(c(alpha / 2, 1 - alpha / 2))
  adj <- stats::pnorm(z0 + (z0 + z_alpha) / (1 - accel * (z0 + z_alpha)))
  adj <- pmin(pmax(adj, 1 / (2 * B)), 1 - 1 / (2 * B))
  quantile_interval(boot_theta, adj)
}

transform_parts <- function(kind) {
  if (kind == "identity") {
    return(list(
      g = function(x) x,
      inv = function(x) x,
      deriv = function(x) rep(1, length(x))
    ))
  }
  if (kind == "arcsine") {
    return(list(
      g = function(x) asin(clip_unit(x)),
      inv = function(x) sin(x),
      deriv = function(x) 1 / sqrt(1 - clip_unit(x)^2)
    ))
  }
  if (kind == "fisher") {
    return(list(
      g = function(x) atanh(clip_unit(x)),
      inv = function(x) tanh(x),
      deriv = function(x) 1 / (1 - clip_unit(x)^2)
    ))
  }
  stop("Unknown transform: ", kind, call. = FALSE)
}

bootstrap_t_interval <- function(theta, se, boot_theta, boot_se, kind) {
  if (!is.finite(theta) || !is.finite(se) || se <= 0) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  parts <- transform_parts(kind)
  ok <- is.finite(boot_theta) & is.finite(boot_se) & boot_se > 0
  if (sum(ok) < 20L) return(c(lower = NA_real_, upper = NA_real_))

  t_boot <- (parts$g(boot_theta[ok]) - parts$g(theta)) /
    (parts$deriv(boot_theta[ok]) * boot_se[ok])
  t_boot <- t_boot[is.finite(t_boot)]
  if (length(t_boot) < 20L) return(c(lower = NA_real_, upper = NA_real_))

  q <- stats::quantile(t_boot, probs = c(alpha / 2, 1 - alpha / 2), names = FALSE, type = 7)
  se_t <- parts$deriv(theta) * se
  limits_t <- parts$g(theta) - c(q[[2L]], q[[1L]]) * se_t
  sort(parts$inv(limits_t))
}

make_bootstrap_intervals <- function(theta, se, boot_theta, boot_se, jack_theta) {
  percentile <- quantile_interval(boot_theta, c(alpha / 2, 1 - alpha / 2))
  basic_q <- quantile_interval(boot_theta, c(alpha / 2, 1 - alpha / 2))
  basic <- if (all(is.finite(basic_q))) {
    c(lower = 2 * theta - basic_q[[2L]], upper = 2 * theta - basic_q[[1L]])
  } else {
    c(lower = NA_real_, upper = NA_real_)
  }

  list(
    boot_percentile = percentile,
    boot_basic = basic,
    boot_bca = bca_interval(theta, boot_theta, jack_theta),
    boot_t = bootstrap_t_interval(theta, se, boot_theta, boot_se, "identity"),
    boot_t_arcsine = bootstrap_t_interval(theta, se, boot_theta, boot_se, "arcsine"),
    boot_t_fisher = bootstrap_t_interval(theta, se, boot_theta, boot_se, "fisher")
  )
}

bootstrap_cache <- function(x, B) {
  n <- nrow(x)
  p <- length(coef_names)
  boot_theta <- matrix(NA_real_, nrow = B, ncol = p, dimnames = list(NULL, coef_names))
  boot_se <- matrix(NA_real_, nrow = B, ncol = p, dimnames = list(NULL, coef_names))
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    fit_b <- quadratic_if_complete(x[idx, , drop = FALSE], values)
    boot_theta[b, ] <- fit_b$estimates[coef_names]
    boot_se[b, ] <- fit_b$se[coef_names]
  }

  jack_theta <- matrix(NA_real_, nrow = n, ncol = p, dimnames = list(NULL, coef_names))
  for (i in seq_len(n)) {
    fit_j <- quadratic_if_complete(x[-i, , drop = FALSE], values)
    jack_theta[i, ] <- fit_j$estimates[coef_names]
  }
  list(theta = boot_theta, se = boot_se, jack = jack_theta)
}

summarise_intervals <- function(intervals, estimates) {
  intervals$covered <- intervals$lower <= intervals$truth & intervals$truth <= intervals$upper
  intervals$length <- intervals$upper - intervals$lower
  intervals$miss_below <- intervals$upper < intervals$truth
  intervals$miss_above <- intervals$lower > intervals$truth

  split_key <- interaction(
    intervals$target, intervals$R, intervals$n, intervals$coefficient, intervals$method,
    drop = TRUE
  )
  pieces <- split(intervals, split_key)
  out <- lapply(pieces, function(d) {
    est <- estimates[estimates$target == d$target[1] &
                       estimates$R == d$R[1] & estimates$n == d$n[1] &
                       estimates$coefficient == d$coefficient[1], ]
    valid <- is.finite(d$lower) & is.finite(d$upper)
    data.frame(
      target = d$target[1],
      R = d$R[1],
      n = d$n[1],
      coefficient = d$coefficient[1],
      method = d$method[1],
      source = d$source[1],
      reps = nrow(d),
      valid_reps = sum(valid),
      coverage = mean(d$covered[valid]),
      miss_below = mean(d$miss_below[valid]),
      miss_above = mean(d$miss_above[valid]),
      mean_length = mean(d$length[valid]),
      median_length = stats::median(d$length[valid]),
      mean_estimate = mean(est$estimate, na.rm = TRUE),
      bias = mean(est$estimate, na.rm = TRUE) - d$target[1],
      mean_se = mean(est$se, na.rm = TRUE),
      mean_boot_reps_used = mean(d$boot_reps_used[valid], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

method_source <- function(method) {
  if (startsWith(method, "moss_")) return("Moss 2024 transform")
  if (method == "boot_bca") return("Efron 1987 BCa")
  if (method %in% c("boot_t", "boot_t_arcsine", "boot_t_fisher")) {
    return("Efron 1987 bootstrap-t")
  }
  "nonparametric bootstrap"
}

estimate_rows <- list()
interval_rows <- list()
estimate_pos <- 1L
interval_pos <- 1L

for (target_idx in seq_along(opts$target_grid)) {
  target <- opts$target_grid[[target_idx]]
  knowledge_prob <- sqrt(target)
  for (R in opts$r_grid) {
    for (n in opts$n_grid) {
      if (opts$progress) {
        message(sprintf(
          "Running target=%.2f R=%d n=%d reps=%d B=%d",
          target, R, n, opts$reps, opts$boot_reps
        ))
      }
      for (rep in seq_len(opts$reps)) {
        seed <- opts$seed_base + 10000000L * target_idx + 100000L * R + 1000L * n + rep
        set.seed(seed)
        x <- simulate_perreault_leigh(n, R, knowledge_prob)
        fit <- quadratic_if_complete(x, values)
        boot <- bootstrap_cache(x, opts$boot_reps)

        for (coefficient in coef_names) {
          theta <- unname(fit$estimates[coefficient])
          se <- unname(fit$se[coefficient])
          estimate_rows[[estimate_pos]] <- data.frame(
            target = target,
            knowledge_prob = knowledge_prob,
            R = R, n = n, rep = rep, seed = seed,
            coefficient = coefficient,
            estimate = theta,
            se = se,
            truth = target,
            stringsAsFactors = FALSE
          )
          estimate_pos <- estimate_pos + 1L

          boot_theta <- boot$theta[, coefficient]
          boot_se <- boot$se[, coefficient]
          jack_theta <- boot$jack[, coefficient]
          boot_reps_used <- sum(is.finite(boot_theta))

          ci_by_method <- list()
          for (method in analytic_methods) {
            ci_by_method[[method]] <- analytic_interval(theta, se, n, method)
          }
          ci_by_method <- c(
            ci_by_method,
            make_bootstrap_intervals(theta, se, boot_theta, boot_se, jack_theta)
          )

          for (method in methods) {
            ci <- ci_by_method[[method]]
            interval_rows[[interval_pos]] <- data.frame(
              target = target,
              knowledge_prob = knowledge_prob,
              R = R, n = n, rep = rep, seed = seed,
              coefficient = coefficient,
              method = method,
              source = method_source(method),
              lower = ci[[1L]],
              upper = ci[[2L]],
              truth = target,
              boot_reps_used = if (startsWith(method, "moss_")) NA_integer_ else boot_reps_used,
              stringsAsFactors = FALSE
            )
            interval_pos <- interval_pos + 1L
          }
        }
      }
    }
  }
}

estimates <- do.call(rbind, estimate_rows)
intervals <- do.call(rbind, interval_rows)
summary <- summarise_intervals(intervals, estimates)
summary <- summary[order(summary$target, summary$coefficient, summary$R, summary$n, summary$method), ]

dir.create(opts$out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(estimates, file.path(opts$out_dir, "estimates.csv"), row.names = FALSE)
write.csv(intervals, file.path(opts$out_dir, "intervals.csv"), row.names = FALSE)
write.csv(summary, file.path(opts$out_dir, "summary.csv"), row.names = FALSE)

metadata <- data.frame(
  key = c(
    "experiment", "n_grid", "r_grid", "reps", "boot_reps", "target_grid",
    "knowledge_prob_grid", "seed_base", "category_count", "methods", "source"
  ),
  value = c(
    "08-quadratic-bootstrap-literature",
    paste(opts$n_grid, collapse = ","),
    paste(opts$r_grid, collapse = ","),
    as.character(opts$reps),
    as.character(opts$boot_reps),
    paste(opts$target_grid, collapse = ","),
    paste(sqrt(opts$target_grid), collapse = ","),
    as.character(opts$seed_base),
    as.character(category_count),
    paste(methods, collapse = ","),
    "Moss 2024 transforms; Efron 1987 bootstrap percentile/basic/BCa/bootstrap-t"
  ),
  stringsAsFactors = FALSE
)
write.csv(metadata, file.path(opts$out_dir, "metadata.csv"), row.names = FALSE)

cat("Wrote:\n")
cat(" ", file.path(opts$out_dir, "estimates.csv"), "\n")
cat(" ", file.path(opts$out_dir, "intervals.csv"), "\n")
cat(" ", file.path(opts$out_dir, "summary.csv"), "\n")
cat(" ", file.path(opts$out_dir, "metadata.csv"), "\n")
