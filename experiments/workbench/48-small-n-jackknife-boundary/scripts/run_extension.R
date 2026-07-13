#!/usr/bin/env Rscript

# Experiment 48 extension: selected intervals under five DGPs with the same
# population quadratic Fleiss kappa but different disagreement geometry.

suppressPackageStartupMessages(library(misskappa))

usage <- function(status = 0L) {
  cat(paste0(
    "Usage: Rscript scripts/run_extension.R [options]\n\n",
    "Options:\n",
    "  --smoke                 Mechanical check: n=20, reps=2, bootstrap B=29.\n",
    "  --reps N                Replications per DGP-by-n cell. Default: 750.\n",
    "  --bootstrap-reps N      Subject bootstrap draws per replication. Default: 199.\n",
    "  --n-grid CSV            Sample sizes. Default: 20,30.\n",
    "  --seed-base N           Deterministic seed base. Default: 484800.\n",
    "  --cores N               Fork workers on non-Windows systems. Default: 1.\n",
    "  --out-dir PATH          Output directory. Default: experiment results/.\n",
    "  --progress              Print one line per DGP-by-n cell.\n",
    "  --help, -h              Show this help and exit.\n"
  ))
  quit(save = "no", status = status)
}

script_dir <- local({
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(arg) == 0L) getwd()
  else dirname(normalizePath(sub("^--file=", "", arg[[1L]]), mustWork = FALSE))
})
experiment_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
source(file.path(experiment_dir, "R", "interval_methods.R"), local = TRUE)

parse_int_csv <- function(x) {
  out <- as.integer(strsplit(x, ",", fixed = TRUE)[[1L]])
  if (any(is.na(out))) stop("--n-grid must contain integers.", call. = FALSE)
  unique(out)
}

parse_args <- function(argv) {
  out <- list(
    smoke = FALSE, reps = 750L, bootstrap_reps = 199L,
    n_grid = c(20L, 30L), seed_base = 484800L, cores = 1L,
    out_dir = file.path(experiment_dir, "results"), progress = FALSE
  )
  i <- 1L
  while (i <= length(argv)) {
    flag <- argv[[i]]
    if (flag %in% c("--help", "-h")) usage()
    if (flag %in% c("--smoke", "--progress")) {
      out[[sub("^--", "", flag)]] <- TRUE
      i <- i + 1L
      next
    }
    if (i == length(argv)) stop(flag, " needs a value.", call. = FALSE)
    value <- argv[[i + 1L]]
    if (flag == "--reps") out$reps <- as.integer(value)
    else if (flag == "--bootstrap-reps") out$bootstrap_reps <- as.integer(value)
    else if (flag == "--n-grid") out$n_grid <- parse_int_csv(value)
    else if (flag == "--seed-base") out$seed_base <- as.integer(value)
    else if (flag == "--cores") out$cores <- as.integer(value)
    else if (flag == "--out-dir") out$out_dir <- value
    else stop("Unknown argument: ", flag, call. = FALSE)
    i <- i + 2L
  }
  if (out$smoke) {
    out$reps <- 2L
    out$bootstrap_reps <- 29L
    out$n_grid <- 20L
    out$cores <- 1L
  }
  if (is.na(out$reps) || out$reps < 1L) stop("--reps must be >= 1.", call. = FALSE)
  if (is.na(out$bootstrap_reps) || out$bootstrap_reps < 20L) {
    stop("--bootstrap-reps must be >= 20.", call. = FALSE)
  }
  if (any(is.na(out$n_grid)) || any(out$n_grid < 8L)) {
    stop("--n-grid values must be >= 8.", call. = FALSE)
  }
  if (is.na(out$seed_base)) stop("--seed-base must be an integer.", call. = FALSE)
  if (is.na(out$cores) || out$cores < 1L) stop("--cores must be >= 1.", call. = FALSE)
  out
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))
C <- 5L
R <- 5L
target <- 0.90
values <- seq_len(C)
score_range <- max(values) - min(values)
selected_methods <- c(
  wald_if_fisher = "IF Wald (Fisher)",
  jack_fisher_direct_boundary = "Direct Fisher jackknife-t + boundary",
  bootstrap_t_jack_boundary = "Bootstrap-t, jackknife SE + boundary"
)

draw <- function(n, probability) {
  sample.int(length(probability), n, replace = TRUE, prob = probability)
}

pl_transition <- function(probability, knowledge) {
  knowledge * diag(C) +
    (1 - knowledge) * matrix(probability, C, C, byrow = TRUE)
}

transition_theory <- function(probability, transition) {
  margin <- as.numeric(probability %*% transition)
  chance <- 2 * (sum(values^2 * margin) - sum(values * margin)^2)
  conditional_mean <- as.numeric(transition %*% values)
  conditional_second <- as.numeric(transition %*% values^2)
  observed <- 2 * sum(
    probability * (conditional_second - conditional_mean^2)
  )
  list(
    kappa = 1 - observed / chance,
    chance = chance,
    observed = observed,
    event_probability = 1 - sum(probability * rowSums(transition^R))
  )
}

simulate_transition <- function(n, probability, transition) {
  truth <- draw(n, probability)
  x <- matrix(NA_integer_, n, R)
  for (j in seq_len(R)) {
    for (category in values) {
      index <- which(truth == category)
      if (length(index) > 0L) {
        x[index, j] <- draw(length(index), transition[category, ])
      }
    }
  }
  x
}

heterogeneous_theory <- function(probability, transitions) {
  margins <- vapply(
    transitions, function(P) as.numeric(probability %*% P), numeric(C)
  )
  margin <- rowMeans(margins)
  chance <- 2 * (sum(values^2 * margin) - sum(values * margin)^2)
  pairs <- utils::combn(R, 2L)
  observed <- mean(vapply(seq_len(ncol(pairs)), function(k) {
    P1 <- transitions[[pairs[1L, k]]]
    P2 <- transitions[[pairs[2L, k]]]
    m1 <- as.numeric(P1 %*% values)
    m2 <- as.numeric(P2 %*% values)
    q1 <- as.numeric(P1 %*% values^2)
    q2 <- as.numeric(P2 %*% values^2)
    sum(probability * (q1 + q2 - 2 * m1 * m2))
  }, numeric(1)))
  all_equal <- sum(vapply(values, function(category) {
    products <- lapply(transitions, function(P) P[, category])
    sum(probability * Reduce(function(x, y) x * y, products))
  }, numeric(1)))
  list(
    kappa = 1 - observed / chance,
    chance = chance,
    observed = observed,
    event_probability = 1 - all_equal
  )
}

make_dgps <- function() {
  uniform <- rep(1 / C, C)
  skewed <- c(0.55, 0.22, 0.12, 0.07, 0.04)
  knowledge <- sqrt(target)
  P_uniform <- pl_transition(uniform, knowledge)
  P_skewed <- pl_transition(skewed, knowledge)

  offsets <- seq(-0.8, 0.8, length.out = R)
  pair_kappa <- function(intercept) {
    p <- stats::plogis(intercept + offsets)
    mean(utils::combn(p, 2L, function(x) prod(x)))
  }
  intercept <- stats::uniroot(
    function(x) pair_kappa(x) - target, c(-10, 20)
  )$root
  rater_knowledge <- stats::plogis(intercept + offsets)
  rater_transitions <- lapply(
    rater_knowledge, function(p) pl_transition(uniform, p)
  )

  local_kernel <- matrix(0, C, C)
  for (category in values) {
    neighbor <- intersect(c(category - 1L, category + 1L), values)
    local_kernel[category, neighbor] <- 1 / length(neighbor)
  }
  local_kappa <- function(error) {
    P <- (1 - error) * diag(C) + error * local_kernel
    transition_theory(uniform, P)$kappa
  }
  local_error <- stats::uniroot(
    function(x) local_kappa(x) - target, c(0, 1)
  )$root
  P_local <- (1 - local_error) * diag(C) + local_error * local_kernel

  failure <- 1 - target
  chance_uniform <- 2 * (
    sum(values^2 * uniform) - sum(values * uniform)^2
  )
  failure_theory <- list(
    kappa = target,
    chance = chance_uniform,
    observed = failure * chance_uniform,
    event_probability = failure * (1 - sum(uniform^R))
  )

  list(
    pl_uniform = list(
      label = "PL: uniform margins",
      mechanism = "Independent knowledge/guessing; uniform truth and guesses",
      margin = uniform, tuning = sprintf("knowledge=%.8f", knowledge),
      theory = transition_theory(uniform, P_uniform),
      simulate = function(n) simulate_transition(n, uniform, P_uniform)
    ),
    pl_skewed = list(
      label = "PL: skewed margins",
      mechanism = "Independent knowledge/guessing; skewed truth and guesses",
      margin = skewed, tuning = sprintf("knowledge=%.8f", knowledge),
      theory = transition_theory(skewed, P_skewed),
      simulate = function(n) simulate_transition(n, skewed, P_skewed)
    ),
    pl_rater_heterogeneous = list(
      label = "PL: heterogeneous raters",
      mechanism = "Independent knowledge/guessing; rater-specific knowledge",
      margin = uniform,
      tuning = paste(
        "knowledge", paste(sprintf("%.6f", rater_knowledge), collapse = "|"),
        sep = "="
      ),
      theory = heterogeneous_theory(uniform, rater_transitions),
      simulate = function(n) {
        truth <- draw(n, uniform)
        x <- matrix(NA_integer_, n, R)
        for (j in seq_len(R)) {
          knows <- stats::runif(n) < rater_knowledge[[j]]
          x[, j] <- ifelse(knows, truth, draw(n, uniform))
        }
        x
      }
    ),
    rare_item_failure = list(
      label = "Rare item-wide failures",
      mechanism = "All raters correct, or an item-wide shock makes all guess",
      margin = uniform, tuning = sprintf("item_failure=%.8f", failure),
      theory = failure_theory,
      simulate = function(n) {
        truth <- draw(n, uniform)
        failed <- stats::runif(n) < failure
        x <- matrix(truth, n, R)
        for (j in seq_len(R)) {
          if (any(failed)) x[failed, j] <- draw(sum(failed), uniform)
        }
        x
      }
    ),
    local_near_miss = list(
      label = "Frequent local near-misses",
      mechanism = "Independent adjacent-category errors conditional on truth",
      margin = uniform, tuning = sprintf("local_error=%.8f", local_error),
      theory = transition_theory(uniform, P_local),
      simulate = function(n) simulate_transition(n, uniform, P_local)
    )
  )
}

dgps <- make_dgps()
dgp_table <- do.call(rbind, lapply(names(dgps), function(name) {
  d <- dgps[[name]]
  data.frame(
    dgp = name, dgp_label = d$label, mechanism = d$mechanism,
    margin = paste(sprintf("%.4f", d$margin), collapse = "|"),
    tuning = d$tuning, true_kappa = d$theory$kappa,
    chance_disagreement = d$theory$chance / score_range^2,
    observed_disagreement = d$theory$observed / score_range^2,
    disagreement_event_probability = d$theory$event_probability,
    positive_event_severity = (
      d$theory$observed / score_range^2
    ) / d$theory$event_probability,
    stringsAsFactors = FALSE
  )
}))
rownames(dgp_table) <- NULL
stopifnot(max(abs(dgp_table$true_kappa - target)) < 2e-6)

one_replicate <- function(dgp_name, dgp_index, n, rep) {
  seed <- opts$seed_base + 10000000L * dgp_index + 100000L * n + rep
  bootstrap_seed <- seed + 480000000L
  started <- proc.time()[["elapsed"]]
  tryCatch({
    set.seed(seed)
    x <- dgps[[dgp_name]]$simulate(n)
    fit <- kappa_counts(
      ratings_to_counts(x, categories = values),
      estimator = "fleiss_cuzick", weight = "quadratic", values = values
    )
    estimate <- unname(stats::coef(fit)[["Fleiss"]])
    se <- sqrt(unname(stats::vcov(fit)["Fleiss", "Fleiss"]))
    focused <- .focused_intervals(
      x, estimate, se, bootstrap_reps = opts$bootstrap_reps,
      min_valid_fraction = 0.70, bootstrap_seed = bootstrap_seed,
      score_range = score_range
    )
    intervals <- focused$intervals[names(selected_methods), , drop = FALSE]
    data.frame(
      dgp = dgp_name, dgp_label = dgps[[dgp_name]]$label,
      n = n, rep = rep, seed = seed, method = rownames(intervals),
      estimate = estimate, if_se = se, jackknife_se = focused$jackknife_se,
      lower = intervals[, "lower"], upper = intervals[, "upper"],
      truth = dgps[[dgp_name]]$theory$kappa,
      fallback_used = unname(focused$fallback_used[rownames(intervals)]),
      disagreement_events = focused$disagreement_events,
      exact_agreement = focused$exact_agreement,
      leave_one_out_boundary = focused$leave_one_out_boundary,
      valid_bootstrap_t_fraction = focused$valid_bootstrap_t_fraction,
      bootstrap_boundary_fraction = focused$bootstrap_boundary_fraction,
      package_delta_estimate = focused$package_delta_estimate,
      package_delta_se = focused$package_delta_se,
      elapsed_seconds = proc.time()[["elapsed"]] - started,
      status = "ok", stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      dgp = dgp_name, dgp_label = dgps[[dgp_name]]$label,
      n = n, rep = rep, seed = seed, method = names(selected_methods),
      estimate = NA_real_, if_se = NA_real_, jackknife_se = NA_real_,
      lower = NA_real_, upper = NA_real_,
      truth = dgps[[dgp_name]]$theory$kappa, fallback_used = FALSE,
      disagreement_events = NA_integer_, exact_agreement = NA,
      leave_one_out_boundary = NA, valid_bootstrap_t_fraction = NA_real_,
      bootstrap_boundary_fraction = NA_real_,
      package_delta_estimate = NA_real_, package_delta_se = NA_real_,
      elapsed_seconds = proc.time()[["elapsed"]] - started,
      status = conditionMessage(e), stringsAsFactors = FALSE
    )
  })
}

run_cell <- function(dgp_name, dgp_index, n) {
  if (opts$progress) {
    message(sprintf(
      "%s, n=%d, reps=%d, bootstrap B=%d, cores=%d",
      dgp_name, n, opts$reps, opts$bootstrap_reps, opts$cores
    ))
  }
  reps <- seq_len(opts$reps)
  worker <- function(rep) one_replicate(dgp_name, dgp_index, n, rep)
  pieces <- if (opts$cores > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(
      reps, worker, mc.cores = opts$cores, mc.preschedule = TRUE
    )
  } else {
    lapply(reps, worker)
  }
  do.call(rbind, pieces)
}

started_all <- proc.time()[["elapsed"]]
pieces <- list()
k <- 1L
for (dgp_index in seq_along(dgps)) {
  for (n in opts$n_grid) {
    pieces[[k]] <- run_cell(names(dgps)[[dgp_index]], dgp_index, n)
    k <- k + 1L
  }
}
replicates <- do.call(rbind, pieces)
elapsed_all <- proc.time()[["elapsed"]] - started_all
replicates$method_label <- unname(selected_methods[replicates$method])
replicates$valid <- replicates$status == "ok" &
  is.finite(replicates$lower) & is.finite(replicates$upper)
replicates$covered <- replicates$valid &
  replicates$lower <= replicates$truth & replicates$truth <= replicates$upper
replicates$miss_below <- replicates$valid & replicates$upper < replicates$truth
replicates$miss_above <- replicates$valid & replicates$lower > replicates$truth
replicates$length <- replicates$upper - replicates$lower

groups <- split(
  replicates,
  interaction(replicates$dgp, replicates$n, replicates$method, drop = TRUE)
)
summary <- do.call(rbind, lapply(groups, function(d) {
  valid <- d$valid
  ok <- d$status == "ok"
  coverage <- mean(d$covered)
  empirical_sd <- stats::sd(d$estimate[ok], na.rm = TRUE)
  data.frame(
    dgp = d$dgp[[1L]], dgp_label = d$dgp_label[[1L]], n = d$n[[1L]],
    method = d$method[[1L]], method_label = d$method_label[[1L]],
    reps = nrow(d), valid_reps = sum(valid), failures = sum(!valid),
    coverage = coverage,
    coverage_mcse = sqrt(coverage * (1 - coverage) / nrow(d)),
    coverage_valid = if (any(valid)) mean(d$covered[valid]) else NA_real_,
    miss_below = if (any(valid)) mean(d$miss_below[valid]) else NA_real_,
    miss_above = if (any(valid)) mean(d$miss_above[valid]) else NA_real_,
    mean_length = if (any(valid)) mean(d$length[valid]) else NA_real_,
    fallback_rate = mean(d$fallback_used[ok], na.rm = TRUE),
    bias = mean(d$estimate[ok], na.rm = TRUE) - d$truth[[1L]],
    empirical_sd = empirical_sd,
    mean_if_se = mean(d$if_se[ok], na.rm = TRUE),
    mean_jackknife_se = mean(d$jackknife_se[ok], na.rm = TRUE),
    if_se_ratio = mean(d$if_se[ok], na.rm = TRUE) / empirical_sd,
    jackknife_se_ratio = mean(d$jackknife_se[ok], na.rm = TRUE) / empirical_sd,
    stringsAsFactors = FALSE
  )
}))
summary <- summary[order(
  match(summary$dgp, names(dgps)), summary$n,
  match(summary$method, names(selected_methods))
), ]
rownames(summary) <- NULL

diag_groups <- split(
  replicates, interaction(replicates$dgp, replicates$n, drop = TRUE)
)
diagnostics <- do.call(rbind, lapply(diag_groups, function(d) {
  one <- d[d$method == "wald_if_fisher" & d$status == "ok", ]
  theory <- dgp_table[dgp_table$dgp == d$dgp[[1L]], ]
  data.frame(
    dgp = d$dgp[[1L]], dgp_label = d$dgp_label[[1L]], n = d$n[[1L]],
    expected_disagreement_event_rate = theory$disagreement_event_probability,
    mean_disagreement_event_rate = mean(one$disagreement_events / d$n[[1L]]),
    exact_agreement_rate = mean(one$exact_agreement),
    leave_one_out_boundary_rate = mean(one$leave_one_out_boundary),
    mean_valid_bootstrap_t_fraction = mean(one$valid_bootstrap_t_fraction),
    mean_bootstrap_boundary_fraction = mean(one$bootstrap_boundary_fraction),
    max_abs_package_delta_estimate = max(abs(one$package_delta_estimate)),
    max_abs_package_delta_se = max(abs(one$package_delta_se)),
    median_seconds_per_replicate = stats::median(one$elapsed_seconds),
    stringsAsFactors = FALSE
  )
}))
diagnostics <- diagnostics[order(
  match(diagnostics$dgp, names(dgps)), diagnostics$n
), ]
rownames(diagnostics) <- NULL

dir.create(opts$out_dir, recursive = TRUE, showWarnings = FALSE)
paths <- file.path(opts$out_dir, c(
  "extension_replicates.csv", "extension_summary.csv",
  "extension_diagnostics.csv", "extension_dgps.csv", "extension_metadata.csv"
))
write.csv(replicates, paths[[1L]], row.names = FALSE)
write.csv(summary, paths[[2L]], row.names = FALSE)
write.csv(diagnostics, paths[[3L]], row.names = FALSE)
write.csv(dgp_table, paths[[4L]], row.names = FALSE)
metadata <- data.frame(
  key = c(
    "experiment", "mode", "reps", "bootstrap_reps", "n_grid", "raters",
    "categories", "target", "seed_base", "cores", "dgps", "methods",
    "elapsed_seconds", "R_version", "misskappa_version",
    "misskappa_library", "generated_at"
  ),
  value = c(
    "48-small-n-jackknife-boundary-extension",
    if (opts$smoke) "smoke" else "extension",
    opts$reps, opts$bootstrap_reps, paste(opts$n_grid, collapse = ","),
    R, C, target, opts$seed_base, opts$cores, paste(names(dgps), collapse = ","),
    paste(names(selected_methods), collapse = ","), elapsed_all,
    R.version.string, as.character(utils::packageVersion("misskappa")),
    find.package("misskappa"), format(Sys.time(), tz = "UTC", usetz = TRUE)
  ),
  stringsAsFactors = FALSE
)
write.csv(metadata, paths[[5L]], row.names = FALSE)

cat("Wrote:\n")
for (path in paths) cat("  ", path, "\n", sep = "")
