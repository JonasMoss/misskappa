#!/usr/bin/env Rscript
#
# Experiment 12: rare-disagreement boundary intervals for quadratic kappa.
#
# The smooth delta corrections in experiment 11 did not repair high-agreement
# coverage. This experiment keeps bootstrap out of the loop and treats observed
# quadratic disagreement as a lattice count. It compares Moss/Fisher baselines
# against Garwood, Jeffreys, and Anscombe count intervals for the observed
# disagreement count, optionally widened by a delta interval for the chance
# disagreement denominator.

usage <- function(status = 0L) {
  cat(paste0(
    "Usage: Rscript run_experiment.R [options]\n\n",
    "Options:\n",
    "  --help              Show this help and exit.\n",
    "  --smoke             Cheap mechanical check.\n",
    "  --n-grid CSV        Sample sizes. Default: 10,20,40,100.\n",
    "  --r-grid CSV        Rater counts. Default: 2,5.\n",
    "  --reps N            Monte Carlo replicates per design cell. Default: 1000.\n",
    "  --target-grid CSV   True quadratic kappas. Default: 0.8,0.9,0.95.\n",
    "  --seed-base N       Base seed for deterministic runs.\n",
    "  --out-dir PATH      Output directory. Default: script-local results/.\n",
    "  --progress          Print one line per design cell.\n\n",
    "Examples:\n",
    "  Rscript run_experiment.R --smoke\n",
    "  Rscript run_experiment.R --reps 500 --target-grid 0.8,0.9,0.95 --n-grid 10,20,40 --r-grid 2,5 --progress\n"
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
    n_grid = c(10L, 20L, 40L, 100L),
    r_grid = c(2L, 5L),
    reps = 1000L,
    target_grid = c(0.8, 0.9, 0.95),
    seed_base = 121212L,
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
    needs_value <- c("--n-grid", "--r-grid", "--reps", "--target-grid", "--seed-base", "--out-dir")
    if (arg %in% needs_value) {
      if (i == length(argv)) stop(arg, " needs a value.", call. = FALSE)
      val <- argv[[i + 1L]]
      if (arg == "--n-grid") opts$n_grid <- parse_int_csv(val, 2L, "--n-grid")
      if (arg == "--r-grid") opts$r_grid <- parse_int_csv(val, 2L, "--r-grid")
      if (arg == "--reps") opts$reps <- as.integer(val)
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
    opts$reps <- 40L
    opts$target_grid <- c(0.9, 0.95)
  }
  if (is.na(opts$reps) || opts$reps < 1L) stop("--reps must be >= 1.", call. = FALSE)
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
range_sq <- (max(values) - min(values))^2
alpha <- 0.05
coef_names <- c("Conger", "Fleiss")
methods <- c(
  "moss_fisher", "delta2_fisher",
  "garwood_plug", "garwood_chance",
  "jeffreys_plug", "jeffreys_chance",
  "anscombe_plug", "anscombe_chance"
)

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
  outer(values, values, function(a, b) (a - b)^2 / range_sq)
}

fit_quadratic_boundary <- function(x, values) {
  n <- nrow(x)
  R <- ncol(x)
  C <- length(values)
  L <- quadratic_loss(values)
  idx <- matrix(match(as.integer(x), values), nrow = n, ncol = R)
  if (any(is.na(idx))) stop("Rating outside value support.", call. = FALSE)

  subject_counts <- matrix(0, nrow = n, ncol = C)
  rater_counts <- matrix(0, nrow = R, ncol = C)
  h_dN <- numeric(n)
  loss_units_by_item <- numeric(n)
  for (i in seq_len(n)) {
    for (j in seq_len(R)) {
      a <- idx[i, j]
      subject_counts[i, a] <- subject_counts[i, a] + 1
      rater_counts[j, a] <- rater_counts[j, a] + 1
    }
    if (R > 1L) {
      for (j in seq_len(R - 1L)) {
        for (k in (j + 1L):R) {
          loss_units <- (idx[i, j] - idx[i, k])^2
          loss_units_by_item[i] <- loss_units_by_item[i] + loss_units
          h_dN[i] <- h_dN[i] + loss_units / range_sq
        }
      }
    }
  }

  pair_count <- choose(R, 2)
  count_denom <- range_sq * pair_count * n
  loss_units_total <- sum(loss_units_by_item)
  psi_dN <- mean(h_dN)
  d_hat <- psi_dN / pair_count
  phi_a <- (h_dN - psi_dN) / pair_count

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

  out <- list()
  for (coefficient in coef_names) {
    B <- if (coefficient == "Conger") d_conger else d_fleiss
    phi_b <- if (coefficient == "Conger") phi_c_chance else phi_f_chance
    theta <- 1 - d_hat / B

    phi_ab <- cbind(A = phi_a, B = phi_b)
    phi_ab <- sweep(phi_ab, 2L, colMeans(phi_ab), "-")
    sigma_n <- crossprod(phi_ab) / n
    sigma_unbiased <- crossprod(phi_ab) / (n - 1)
    grad <- c(-1 / B, d_hat / (B^2))
    H <- matrix(c(0, 1 / (B^2), 1 / (B^2), -2 * d_hat / (B^3)), nrow = 2L)

    v_if <- as.numeric(t(grad) %*% sigma_n %*% grad) / n
    tr_delta <- sum((H %*% sigma_unbiased) * t(H %*% sigma_unbiased))
    v_delta2 <- as.numeric(t(grad) %*% sigma_unbiased %*% grad) / n +
      max(0, 0.5 * tr_delta / (n^2))
    bias2 <- 0.5 * sum(H * sigma_unbiased) / n

    phi_b <- phi_b - mean(phi_b)
    se_b <- sqrt(mean(phi_b^2) / n)

    out[[coefficient]] <- data.frame(
      coefficient = coefficient,
      estimate = theta,
      observed_disagreement = d_hat,
      chance_disagreement = B,
      loss_units_total = loss_units_total,
      count_denom = count_denom,
      zero_disagreement = as.integer(loss_units_total == 0),
      mean_loss_units_item = mean(loss_units_by_item),
      se_if_n = sqrt(max(0, v_if)),
      se_delta2_var = sqrt(max(0, v_delta2)),
      bias_delta2 = bias2,
      chance_se = se_b,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

fisher_interval <- function(theta, se, n, center_shift = 0) {
  if (!is.finite(theta) || !is.finite(se) || se < 0) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  center <- theta - center_shift
  center_c <- clip_unit(center)
  se_t <- se / (1 - center_c^2)
  c_t <- stats::qt(1 - alpha / 2, df = n - 1)
  sort(tanh(atanh(center_c) + c(-1, 1) * c_t * se_t))
}

garwood_count_interval <- function(m) {
  lower <- if (m <= 0) 0 else 0.5 * stats::qchisq(alpha / 2, df = 2 * m)
  upper <- 0.5 * stats::qchisq(1 - alpha / 2, df = 2 * (m + 1))
  c(lower = lower, upper = upper)
}

jeffreys_count_interval <- function(m) {
  stats::qgamma(c(alpha / 2, 1 - alpha / 2), shape = m + 0.5, rate = 1)
}

anscombe_count_interval <- function(m) {
  z <- stats::qnorm(1 - alpha / 2)
  center <- sqrt(m + 3 / 8)
  lower <- max(0, center - z / 2)^2
  upper <- (center + z / 2)^2
  c(lower = lower, upper = upper)
}

count_to_kappa_interval <- function(count_ci, count_denom, B_lower, B_upper) {
  eps <- 1e-10
  B_lower <- max(eps, B_lower)
  B_upper <- max(eps, B_upper)
  A_lower <- count_ci[[1L]] / count_denom
  A_upper <- count_ci[[2L]] / count_denom
  lower <- 1 - A_upper / B_lower
  upper <- 1 - A_lower / B_upper
  c(lower = max(-1, min(1, lower)), upper = max(-1, min(1, upper)))
}

chance_bounds <- function(B, se_b, n, use_chance) {
  if (!use_chance || !is.finite(se_b) || se_b <= 0) return(c(lower = B, upper = B))
  c_t <- stats::qt(1 - alpha / 2, df = n - 1)
  c(lower = max(1e-10, B - c_t * se_b), upper = max(1e-10, B + c_t * se_b))
}

make_interval <- function(fit_row, n, method) {
  theta <- fit_row$estimate
  if (method == "moss_fisher") {
    return(fisher_interval(theta, fit_row$se_if_n, n))
  }
  if (method == "delta2_fisher") {
    return(fisher_interval(theta, fit_row$se_delta2_var, n, center_shift = fit_row$bias_delta2))
  }

  count_method <- sub("_(plug|chance)$", "", method)
  use_chance <- grepl("_chance$", method)
  count_ci <- switch(
    count_method,
    garwood = garwood_count_interval(fit_row$loss_units_total),
    jeffreys = jeffreys_count_interval(fit_row$loss_units_total),
    anscombe = anscombe_count_interval(fit_row$loss_units_total),
    stop("Unknown count method: ", method, call. = FALSE)
  )
  B_ci <- chance_bounds(fit_row$chance_disagreement, fit_row$chance_se, n, use_chance)
  count_to_kappa_interval(count_ci, fit_row$count_denom, B_ci[[1L]], B_ci[[2L]])
}

method_source <- function(method) {
  if (method == "moss_fisher") return("Moss 2024 Fisher transform")
  if (method == "delta2_fisher") return("Moss 2024 Fisher transform; second-order delta variance")
  if (startsWith(method, "garwood")) return("Garwood exact Poisson count interval")
  if (startsWith(method, "jeffreys")) return("Jeffreys Poisson count interval")
  if (startsWith(method, "anscombe")) return("Anscombe Poisson variance-stabilizing interval")
  "unknown"
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
      zero_disagreement_rate = mean(est$zero_disagreement, na.rm = TRUE),
      mean_loss_units = mean(est$loss_units_total, na.rm = TRUE),
      mean_chance_disagreement = mean(est$chance_disagreement, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
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
        message(sprintf("Running target=%.2f R=%d n=%d reps=%d", target, R, n, opts$reps))
      }
      for (rep in seq_len(opts$reps)) {
        seed <- opts$seed_base + 10000000L * target_idx + 100000L * R + 1000L * n + rep
        set.seed(seed)
        x <- simulate_perreault_leigh(n, R, knowledge_prob)
        fit <- fit_quadratic_boundary(x, values)

        for (row_idx in seq_len(nrow(fit))) {
          coefficient <- fit$coefficient[[row_idx]]
          fit_row <- fit[row_idx, ]
          estimate_rows[[estimate_pos]] <- data.frame(
            target = target,
            knowledge_prob = knowledge_prob,
            R = R, n = n, rep = rep, seed = seed,
            coefficient = coefficient,
            estimate = fit_row$estimate,
            observed_disagreement = fit_row$observed_disagreement,
            chance_disagreement = fit_row$chance_disagreement,
            loss_units_total = fit_row$loss_units_total,
            count_denom = fit_row$count_denom,
            zero_disagreement = fit_row$zero_disagreement,
            mean_loss_units_item = fit_row$mean_loss_units_item,
            se_if_n = fit_row$se_if_n,
            se_delta2_var = fit_row$se_delta2_var,
            bias_delta2 = fit_row$bias_delta2,
            chance_se = fit_row$chance_se,
            truth = target,
            stringsAsFactors = FALSE
          )
          estimate_pos <- estimate_pos + 1L

          for (method in methods) {
            ci <- make_interval(fit_row, n, method)
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
    "experiment", "n_grid", "r_grid", "reps", "target_grid", "knowledge_prob_grid",
    "seed_base", "category_count", "methods", "source"
  ),
  value = c(
    "12-quadratic-rare-disagreement",
    paste(opts$n_grid, collapse = ","),
    paste(opts$r_grid, collapse = ","),
    as.character(opts$reps),
    paste(opts$target_grid, collapse = ","),
    paste(sqrt(opts$target_grid), collapse = ","),
    as.character(opts$seed_base),
    as.character(category_count),
    paste(methods, collapse = ","),
    "Moss 2024 Fisher transform; rare-disagreement Poisson count intervals"
  ),
  stringsAsFactors = FALSE
)
write.csv(metadata, file.path(opts$out_dir, "metadata.csv"), row.names = FALSE)

cat("Wrote:\n")
cat(" ", file.path(opts$out_dir, "estimates.csv"), "\n")
cat(" ", file.path(opts$out_dir, "intervals.csv"), "\n")
cat(" ", file.path(opts$out_dir, "summary.csv"), "\n")
cat(" ", file.path(opts$out_dir, "metadata.csv"), "\n")
