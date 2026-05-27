#!/usr/bin/env Rscript
#
# Experiment 07: quadratic kappa coverage and Edgeworth-like intervals.
#
# Complete categorical Perreault-Leigh simulations for quadratically weighted
# kappa. Reproduces the Moss (2024) basic/arcsine/Fisher comparison and adds
# simple Cornish-Fisher corrections driven by estimated influence skewness.

usage <- function(status = 0L) {
  cat(paste0(
    "Usage: Rscript run_experiment.R [options]\n\n",
    "Options:\n",
    "  --help              Show this help and exit.\n",
    "  --smoke             Cheap mechanical check.\n",
    "  --n-grid CSV        Sample sizes. Default: 10,20,40,100.\n",
    "  --r-grid CSV        Rater counts. Default: 2,5.\n",
    "  --reps N            Replicates per design cell. Default: 1000.\n",
    "  --target K          True quadratic kappa. Default: 0.8.\n",
    "  --seed-base N       Base seed for deterministic runs.\n",
    "  --out-dir PATH      Output directory. Default: script-local results/.\n",
    "  --progress          Print one line per design cell.\n\n",
    "Examples:\n",
    "  Rscript run_experiment.R --smoke\n",
    "  Rscript run_experiment.R --reps 400 --n-grid 10,20,40 --r-grid 2,5 --progress\n"
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

parse_args <- function(argv) {
  opts <- list(
    smoke = FALSE,
    n_grid = c(10L, 20L, 40L, 100L),
    r_grid = c(2L, 5L),
    reps = 1000L,
    target = 0.8,
    seed_base = 707070L,
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
    needs_value <- c("--n-grid", "--r-grid", "--reps", "--target", "--seed-base", "--out-dir")
    if (arg %in% needs_value) {
      if (i == length(argv)) stop(arg, " needs a value.", call. = FALSE)
      val <- argv[[i + 1L]]
      if (arg == "--n-grid") opts$n_grid <- parse_int_csv(val, 2L, "--n-grid")
      if (arg == "--r-grid") opts$r_grid <- parse_int_csv(val, 2L, "--r-grid")
      if (arg == "--reps") opts$reps <- as.integer(val)
      if (arg == "--target") opts$target <- as.numeric(val)
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
  }
  if (is.na(opts$reps) || opts$reps < 1L) stop("--reps must be >= 1.", call. = FALSE)
  if (is.na(opts$target) || opts$target <= 0 || opts$target >= 1) {
    stop("--target must be in (0, 1).", call. = FALSE)
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

suppressPackageStartupMessages({
  library(misskappa)
})

category_count <- 5L
values <- seq_len(category_count)
alpha <- 0.05
coef_names <- c("Conger", "Fleiss", "Brennan-Prediger")
methods <- c("basic", "arcsine", "fisher", "cf_if", "cf_student_plus", "cf_student_minus")

clip_unit <- function(x) pmin(pmax(x, -1 + 1e-10), 1 - 1e-10)

simulate_perreault_leigh <- function(n, R, knowledge_prob) {
  truth <- sample.int(category_count, n, replace = TRUE)
  x <- matrix(NA_integer_, nrow = n, ncol = R)
  for (j in seq_len(R)) {
    knows <- stats::runif(n) < knowledge_prob
    guesses <- sample.int(category_count, n, replace = TRUE)
    x[, j] <- ifelse(knows, truth, guesses)
  }
  x
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

  d_bp <- mean(L)
  estimates <- c(
    Conger = 1 - d_hat / d_conger,
    Fleiss = 1 - d_hat / d_fleiss,
    `Brennan-Prediger` = 1 - d_hat / d_bp
  )
  phi <- cbind(
    Conger = -phi_d / d_conger + d_hat * phi_c_chance / (d_conger^2),
    Fleiss = -phi_d / d_fleiss + d_hat * phi_f_chance / (d_fleiss^2),
    `Brennan-Prediger` = -phi_d / d_bp
  )
  phi <- sweep(phi, 2L, colMeans(phi), "-")
  se <- sqrt(colMeans(phi^2) / n)
  skew <- colMeans(phi^3) / (colMeans(phi^2)^(3 / 2))
  list(estimates = estimates, se = se, skew = skew)
}

make_interval <- function(theta, se, skew, n, method) {
  c_t <- stats::qt(1 - alpha / 2, df = n - 1)
  if (!is.finite(se) || se < 0) return(c(lower = NA_real_, upper = NA_real_))
  if (method == "basic") return(c(lower = theta - c_t * se, upper = theta + c_t * se))

  theta_c <- clip_unit(theta)
  if (method == "arcsine") {
    se_t <- se / sqrt(1 - theta_c^2)
    return(sort(sin(asin(theta_c) + c(-1, 1) * c_t * se_t)))
  }
  if (method == "fisher") {
    se_t <- se / (1 - theta_c^2)
    return(sort(tanh(atanh(theta_c) + c(-1, 1) * c_t * se_t)))
  }

  z <- stats::qnorm(c(alpha / 2, 1 - alpha / 2))
  gamma_n <- skew / sqrt(n)
  if (method == "cf_if") {
    q <- z + (gamma_n / 6) * (z^2 - 1)
  } else if (method == "cf_student_plus") {
    q <- z + (gamma_n / 6) * (2 * z^2 + 1)
  } else if (method == "cf_student_minus") {
    q <- z - (gamma_n / 6) * (2 * z^2 + 1)
  } else {
    stop("Unknown interval method: ", method, call. = FALSE)
  }
  c(lower = theta - q[2] * se, upper = theta - q[1] * se)
}

summarise_intervals <- function(intervals, estimates) {
  intervals$covered <- intervals$lower <= intervals$truth & intervals$truth <= intervals$upper
  intervals$length <- intervals$upper - intervals$lower
  split_key <- interaction(intervals$R, intervals$n, intervals$coefficient, intervals$method, drop = TRUE)
  pieces <- split(intervals, split_key)
  out <- lapply(pieces, function(d) {
    est <- estimates[estimates$R == d$R[1] & estimates$n == d$n[1] &
                       estimates$coefficient == d$coefficient[1], ]
    valid <- is.finite(d$lower) & is.finite(d$upper)
    data.frame(
      R = d$R[1],
      n = d$n[1],
      coefficient = d$coefficient[1],
      method = d$method[1],
      reps = nrow(d),
      valid_reps = sum(valid),
      coverage = mean(d$covered[valid]),
      mean_length = mean(d$length[valid]),
      median_length = stats::median(d$length[valid]),
      mean_estimate = mean(est$estimate, na.rm = TRUE),
      sd_estimate = stats::sd(est$estimate, na.rm = TRUE),
      mean_se = mean(est$se, na.rm = TRUE),
      mean_skew_if = mean(est$skew_if, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

knowledge_prob <- sqrt(opts$target)
estimate_rows <- list()
interval_rows <- list()
estimate_pos <- 1L
interval_pos <- 1L

for (R in opts$r_grid) {
  for (n in opts$n_grid) {
    if (opts$progress) {
      message(sprintf("Running R=%d n=%d reps=%d", R, n, opts$reps))
    }
    for (rep in seq_len(opts$reps)) {
      seed <- opts$seed_base + 100000L * R + 1000L * n + rep
      set.seed(seed)
      x <- simulate_perreault_leigh(n, R, knowledge_prob)

      infl <- quadratic_if_complete(x, values)

      for (coefficient in coef_names) {
        theta <- unname(infl$estimates[coefficient])
        se <- unname(infl$se[coefficient])
        skew <- unname(infl$skew[coefficient])
        estimate_rows[[estimate_pos]] <- data.frame(
          R = R, n = n, rep = rep, seed = seed,
          coefficient = coefficient,
          estimate = theta,
          se = se,
          if_estimate = theta,
          if_se = se,
          skew_if = skew,
          truth = opts$target,
          stringsAsFactors = FALSE
        )
        estimate_pos <- estimate_pos + 1L

        for (method in methods) {
          ci <- make_interval(theta, se, skew, n, method)
          interval_rows[[interval_pos]] <- data.frame(
            R = R, n = n, rep = rep, seed = seed,
            coefficient = coefficient,
            method = method,
            lower = ci[[1L]],
            upper = ci[[2L]],
            truth = opts$target,
            stringsAsFactors = FALSE
          )
          interval_pos <- interval_pos + 1L
        }
      }
    }
  }
}

estimates <- do.call(rbind, estimate_rows)
intervals <- do.call(rbind, interval_rows)
summary <- summarise_intervals(intervals, estimates)
summary <- summary[order(summary$coefficient, summary$R, summary$n, summary$method), ]

dir.create(opts$out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(estimates, file.path(opts$out_dir, "estimates.csv"), row.names = FALSE)
write.csv(intervals, file.path(opts$out_dir, "intervals.csv"), row.names = FALSE)
write.csv(summary, file.path(opts$out_dir, "summary.csv"), row.names = FALSE)

metadata <- data.frame(
  key = c(
    "experiment", "n_grid", "r_grid", "reps", "target", "knowledge_prob",
    "seed_base", "category_count", "methods", "source"
  ),
  value = c(
    "07-quadratic-edgeworth-coverage",
    paste(opts$n_grid, collapse = ","),
    paste(opts$r_grid, collapse = ","),
    as.character(opts$reps),
    as.character(opts$target),
    as.character(knowledge_prob),
    as.character(opts$seed_base),
    as.character(category_count),
    paste(methods, collapse = ","),
    "Moss 2024 Perreault-Leigh complete categorical quadratic setup"
  ),
  stringsAsFactors = FALSE
)
write.csv(metadata, file.path(opts$out_dir, "metadata.csv"), row.names = FALSE)

cat("Wrote:\n")
cat("  ", file.path(opts$out_dir, "estimates.csv"), "\n", sep = "")
cat("  ", file.path(opts$out_dir, "intervals.csv"), "\n", sep = "")
cat("  ", file.path(opts$out_dir, "summary.csv"), "\n", sep = "")
cat("  ", file.path(opts$out_dir, "metadata.csv"), "\n", sep = "")
