#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(misskappa)
})

script_dir <- local({
  cmd <- commandArgs(FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) == 0L) getwd()
  else dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE))
})

usage <- function() {
  cat(
    "Usage: Rscript run_experiment.R [options]\n",
    "\n",
    "Options:\n",
    "  --help, -h           Show this help text.\n",
    "  --smoke              Cheap verification run. Defaults to n=150, reps=3.\n",
    "  --n N                Subjects per Monte Carlo replicate. Default: 500.\n",
    "  --reps B             Monte Carlo replicates per mechanism. Default: 30.\n",
    "  --seed-base S        Integer seed base. Default: 41004.\n",
    "  --out-dir DIR        Output directory. Default: this experiment's results/.\n",
    "  --weight NAME        misskappa weight name. Default: identity.\n",
    "\n",
    "Writes metadata.csv, replicates.csv, summary.csv, and category_shift.csv.\n",
    sep = ""
  )
}

parse_args <- function(argv) {
  args <- list(
    smoke = FALSE,
    n = 500L,
    reps = 30L,
    seed_base = 41004L,
    out_dir = file.path(script_dir, "results"),
    weight = "identity"
  )
  explicit <- list(n = FALSE, reps = FALSE)

  set_value <- function(key, value) {
    key <- gsub("-", "_", key, fixed = TRUE)
    if (key == "n") {
      args$n <<- as.integer(value)
      explicit$n <<- TRUE
    } else if (key == "reps") {
      args$reps <<- as.integer(value)
      explicit$reps <<- TRUE
    } else if (key == "seed_base") {
      args$seed_base <<- as.integer(value)
    } else if (key == "out_dir") {
      args$out_dir <<- value
    } else if (key == "weight") {
      args$weight <<- value
    } else {
      stop("Unknown option: --", key)
    }
  }

  i <- 1L
  while (i <= length(argv)) {
    arg <- argv[[i]]
    if (arg %in% c("--help", "-h")) {
      usage()
      quit(save = "no", status = 0L)
    } else if (arg == "--smoke") {
      args$smoke <- TRUE
    } else if (grepl("^--[^=]+=", arg)) {
      key <- sub("^--([^=]+)=.*$", "\\1", arg)
      value <- sub("^--[^=]+=", "", arg)
      set_value(key, value)
    } else if (grepl("^--", arg)) {
      key <- sub("^--", "", arg)
      if (i == length(argv)) stop("Missing value for option: ", arg)
      i <- i + 1L
      set_value(key, argv[[i]])
    } else {
      stop("Unexpected argument: ", arg)
    }
    i <- i + 1L
  }

  if (args$smoke) {
    if (!explicit$n) args$n <- 150L
    if (!explicit$reps) args$reps <- 3L
  }
  if (is.na(args$n) || args$n < 10L) stop("'n' must be an integer >= 10.")
  if (is.na(args$reps) || args$reps < 1L) stop("'reps' must be an integer >= 1.")
  if (is.na(args$seed_base)) stop("'seed-base' must be an integer.")
  args
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

C <- 4L
R <- 5L
p_truth <- c(0.20, 0.30, 0.30, 0.20)
coef_names <- c("Fleiss", "Brennan-Prediger")
em_options <- list(max_iter = 20000L, tol = 1e-7)

make_dgp <- function(mechanism) {
  if (mechanism == "m1_rater_mcar") {
    list(
      mechanism = mechanism,
      label = "M1: heterogeneous raters, rater MCAR",
      seed_offset = 1000L,
      p_truth = p_truth,
      skill = c(0.96, 0.91, 0.84, 0.73, 0.62),
      guess = rbind(
        c(0.06, 0.16, 0.31, 0.47),
        c(0.11, 0.23, 0.33, 0.33),
        c(0.20, 0.30, 0.30, 0.20),
        c(0.32, 0.32, 0.22, 0.14),
        c(0.48, 0.29, 0.16, 0.07)
      ),
      observe_prob_rater = c(0.95, 0.83, 0.64, 0.39, 0.24),
      observe_prob_value = rep(NA_real_, C)
    )
  } else if (mechanism == "m2_value_dropout") {
    list(
      mechanism = mechanism,
      label = "M2: value-dependent dropout",
      seed_offset = 2000L,
      p_truth = p_truth,
      skill = rep(0.86, R),
      guess = matrix(rep(p_truth, each = R), nrow = R, byrow = TRUE),
      observe_prob_rater = rep(NA_real_, R),
      observe_prob_value = c(0.35, 0.55, 0.75, 0.92)
    )
  } else {
    stop("Unknown mechanism: ", mechanism)
  }
}

dgps <- list(
  m1_rater_mcar = make_dgp("m1_rater_mcar"),
  m2_value_dropout = make_dgp("m2_value_dropout")
)

simulate_complete <- function(n, dgp) {
  truth <- sample.int(C, n, replace = TRUE, prob = dgp$p_truth)
  x_star <- matrix(NA_integer_, nrow = n, ncol = R)
  for (j in seq_len(R)) {
    correct <- stats::runif(n) < dgp$skill[[j]]
    guessed <- sample.int(C, n, replace = TRUE, prob = dgp$guess[j, ])
    x_star[, j] <- ifelse(correct, truth, guessed)
  }
  list(x_star = x_star, truth = truth)
}

apply_missing <- function(x_star, dgp) {
  keep <- matrix(FALSE, nrow = nrow(x_star), ncol = ncol(x_star))
  if (dgp$mechanism == "m1_rater_mcar") {
    for (j in seq_len(R)) {
      keep[, j] <- stats::runif(nrow(x_star)) < dgp$observe_prob_rater[[j]]
    }
  } else if (dgp$mechanism == "m2_value_dropout") {
    prob <- matrix(dgp$observe_prob_value[as.integer(x_star)], nrow = nrow(x_star))
    keep <- matrix(stats::runif(length(x_star)), nrow = nrow(x_star)) < prob
  } else {
    stop("Unknown mechanism: ", dgp$mechanism)
  }
  x <- x_star
  x[!keep] <- NA_integer_
  x
}

ratings_to_counts <- function(x, C) {
  counts <- matrix(0L, nrow = nrow(x), ncol = C)
  for (k in seq_len(C)) counts[, k] <- rowSums(x == k, na.rm = TRUE)
  colnames(counts) <- paste0("cat", seq_len(C))
  counts
}

loss_matrix <- function(weight, C) {
  values <- seq_len(C)
  if (weight %in% c("identity", "unweighted")) {
    L <- matrix(1, nrow = C, ncol = C)
    diag(L) <- 0
    L
  } else if (weight == "linear") {
    abs(outer(values, values, "-")) / (C - 1)
  } else if (weight == "quadratic") {
    (outer(values, values, "-") / (C - 1))^2
  } else {
    stop("Analytic truth in this runner supports identity, linear, and quadratic weights.")
  }
}

estimate_truth <- function(dgp) {
  L <- loss_matrix(args$weight, C)
  q_by_truth <- array(0, dim = c(C, R, C))
  for (t in seq_len(C)) {
    for (j in seq_len(R)) {
      q_by_truth[t, j, ] <- (1 - dgp$skill[[j]]) * dgp$guess[j, ]
      q_by_truth[t, j, t] <- q_by_truth[t, j, t] + dgp$skill[[j]]
    }
  }

  d_pair <- 0
  n_pairs <- 0L
  for (j in seq_len(R - 1L)) {
    for (k in (j + 1L):R) {
      d_jk <- 0
      for (t in seq_len(C)) {
        d_jk <- d_jk + dgp$p_truth[[t]] *
          sum(outer(q_by_truth[t, j, ], q_by_truth[t, k, ]) * L)
      }
      d_pair <- d_pair + d_jk
      n_pairs <- n_pairs + 1L
    }
  }
  d_hat <- d_pair / n_pairs

  marginal_by_rater <- matrix(0, nrow = R, ncol = C)
  for (j in seq_len(R)) {
    for (t in seq_len(C)) {
      marginal_by_rater[j, ] <- marginal_by_rater[j, ] +
        dgp$p_truth[[t]] * q_by_truth[t, j, ]
    }
  }
  pooled <- colMeans(marginal_by_rater)
  d_fleiss <- sum(outer(pooled, pooled) * L)
  d_bp <- sum(L) / (C * C)

  estimates <- c(
    Fleiss = 1 - d_hat / d_fleiss,
    `Brennan-Prediger` = 1 - d_hat / d_bp
  )
  data.frame(
    mechanism = dgp$mechanism,
    coefficient = coef_names,
    truth = as.numeric(estimates[coef_names]),
    stringsAsFactors = FALSE
  )
}

fit_estimator <- function(estimator, x, counts) {
  tryCatch({
    estimates <- switch(
      estimator,
      raw_ipw = misskappa::kappa(x, method = "ipw", weight = args$weight)$estimates[coef_names],
      raw_fiml = misskappa::kappa(
        x, method = "fiml", weight = args$weight, em_options = em_options
      )$estimates[coef_names],
      counts_fiml = misskappa::kappa_counts(
        counts, method = "fiml", weight = args$weight, r_total = R,
        em_options = em_options
      )$estimates[coef_names],
      stop("Unknown estimator: ", estimator)
    )
    data.frame(
      estimator = estimator,
      coefficient = names(estimates),
      estimate = as.numeric(estimates),
      error = "",
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      estimator = estimator,
      coefficient = coef_names,
      estimate = NA_real_,
      error = conditionMessage(e),
      stringsAsFactors = FALSE
    )
  })
}

category_diagnostics <- function(x_star, x, dgp, rep_id) {
  full_tab <- tabulate(as.integer(x_star), nbins = C) / length(x_star)
  observed_values <- as.integer(x[!is.na(x)])
  obs_tab <- if (length(observed_values) == 0L) {
    rep(NA_real_, C)
  } else {
    tabulate(observed_values, nbins = C) / length(observed_values)
  }
  data.frame(
    mechanism = dgp$mechanism,
    replicate = rep_id,
    category = seq_len(C),
    full_share = full_tab,
    observed_share = obs_tab,
    shift = obs_tab - full_tab,
    stringsAsFactors = FALSE
  )
}

run_replicate <- function(dgp, rep_id, truth) {
  set.seed(args$seed_base + dgp$seed_offset + rep_id)
  complete <- simulate_complete(args$n, dgp)
  x <- apply_missing(complete$x_star, dgp)
  counts <- ratings_to_counts(x, C)

  fits <- do.call(rbind, lapply(c("counts_fiml", "raw_fiml", "raw_ipw"), fit_estimator,
                                x = x, counts = counts))
  truth_match <- truth$truth[match(fits$coefficient, truth$coefficient)]
  fits$mechanism <- dgp$mechanism
  fits$replicate <- rep_id
  fits$n <- args$n
  fits$truth <- truth_match
  fits$bias <- fits$estimate - fits$truth
  fits$observed_fraction <- mean(!is.na(x))
  fits$mean_observed_raters <- mean(rowSums(!is.na(x)))
  fits$zero_observed_fraction <- mean(rowSums(!is.na(x)) == 0L)
  fits <- fits[, c(
    "mechanism", "replicate", "n", "estimator", "coefficient", "estimate",
    "truth", "bias", "observed_fraction", "mean_observed_raters",
    "zero_observed_fraction", "error"
  )]

  list(
    fits = fits,
    category = category_diagnostics(complete$x_star, x, dgp, rep_id)
  )
}

summarise_replicates <- function(df) {
  groups <- split(df, interaction(df$mechanism, df$estimator, df$coefficient, drop = TRUE))
  rows <- lapply(groups, function(d) {
    ok <- !is.na(d$estimate)
    est <- d$estimate[ok]
    bias <- d$bias[ok]
    data.frame(
      mechanism = d$mechanism[[1L]],
      estimator = d$estimator[[1L]],
      coefficient = d$coefficient[[1L]],
      n = d$n[[1L]],
      reps = nrow(d),
      n_success = sum(ok),
      failures = sum(!ok),
      truth = d$truth[[1L]],
      mean_estimate = if (length(est)) mean(est) else NA_real_,
      bias = if (length(bias)) mean(bias) else NA_real_,
      sd = if (length(est) > 1L) stats::sd(est) else NA_real_,
      rmse = if (length(bias)) sqrt(mean(bias^2)) else NA_real_,
      mc_se_bias = if (length(est) > 1L) stats::sd(est) / sqrt(length(est)) else NA_real_,
      mean_observed_fraction = mean(d$observed_fraction),
      mean_observed_raters = mean(d$mean_observed_raters),
      mean_zero_observed_fraction = mean(d$zero_observed_fraction),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[order(out$mechanism, out$coefficient, out$estimator), ]
}

summarise_category_shift <- function(df) {
  groups <- split(df, interaction(df$mechanism, df$category, drop = TRUE))
  rows <- lapply(groups, function(d) {
    data.frame(
      mechanism = d$mechanism[[1L]],
      category = d$category[[1L]],
      full_share = mean(d$full_share),
      observed_share = mean(d$observed_share),
      shift = mean(d$shift),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[order(out$mechanism, out$category), ]
}

dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)

truth <- do.call(rbind, lapply(dgps, estimate_truth))
replicate_rows <- list()
category_rows <- list()

for (dgp_name in names(dgps)) {
  dgp <- dgps[[dgp_name]]
  dgp_truth <- truth[truth$mechanism == dgp$mechanism, ]
  cat(sprintf("Running %s: n=%d, reps=%d\n", dgp$label, args$n, args$reps))
  for (b in seq_len(args$reps)) {
    if (b == 1L || b == args$reps || b %% max(1L, floor(args$reps / 5L)) == 0L) {
      cat(sprintf("  replicate %d/%d\n", b, args$reps))
    }
    out <- run_replicate(dgp, b, dgp_truth)
    replicate_rows[[length(replicate_rows) + 1L]] <- out$fits
    category_rows[[length(category_rows) + 1L]] <- out$category
  }
}

replicates <- do.call(rbind, replicate_rows)
summary <- summarise_replicates(replicates)
category_shift <- summarise_category_shift(do.call(rbind, category_rows))

mechanisms <- do.call(rbind, lapply(dgps, function(dgp) {
  data.frame(
    mechanism = dgp$mechanism,
    label = dgp$label,
    p_truth = paste(dgp$p_truth, collapse = " "),
    skill = paste(dgp$skill, collapse = " "),
    guess = paste(apply(dgp$guess, 1L, paste, collapse = " "), collapse = " | "),
    observe_prob_rater = paste(dgp$observe_prob_rater, collapse = " "),
    observe_prob_value = paste(dgp$observe_prob_value, collapse = " "),
    stringsAsFactors = FALSE
  )
}))

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) == 0L) {
  file.path(script_dir, "run_experiment.R")
} else {
  sub("^--file=", "", script_arg[[1L]])
}

metadata <- data.frame(
  key = c(
    "generated_at", "script", "smoke", "n", "reps", "truth_source", "seed_base",
    "C", "R", "weight", "em_max_iter", "em_tol", "misskappa_version",
    "r_version"
  ),
  value = c(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
    normalizePath(script_path, mustWork = FALSE),
    as.character(args$smoke),
    as.character(args$n),
    as.character(args$reps),
    "analytic",
    as.character(args$seed_base),
    as.character(C),
    as.character(R),
    args$weight,
    as.character(em_options$max_iter),
    as.character(em_options$tol),
    as.character(utils::packageVersion("misskappa")),
    R.version.string
  ),
  stringsAsFactors = FALSE
)

write.csv(metadata, file.path(args$out_dir, "metadata.csv"), row.names = FALSE)
write.csv(mechanisms, file.path(args$out_dir, "mechanisms.csv"), row.names = FALSE)
write.csv(truth, file.path(args$out_dir, "truth.csv"), row.names = FALSE)
write.csv(replicates, file.path(args$out_dir, "replicates.csv"), row.names = FALSE)
write.csv(summary, file.path(args$out_dir, "summary.csv"), row.names = FALSE)
write.csv(category_shift, file.path(args$out_dir, "category_shift.csv"), row.names = FALSE)

cat("Wrote:\n")
cat(sprintf("  %s\n", file.path(args$out_dir, "metadata.csv")))
cat(sprintf("  %s\n", file.path(args$out_dir, "mechanisms.csv")))
cat(sprintf("  %s\n", file.path(args$out_dir, "truth.csv")))
cat(sprintf("  %s\n", file.path(args$out_dir, "replicates.csv")))
cat(sprintf("  %s\n", file.path(args$out_dir, "summary.csv")))
cat(sprintf("  %s\n", file.path(args$out_dir, "category_shift.csv")))
