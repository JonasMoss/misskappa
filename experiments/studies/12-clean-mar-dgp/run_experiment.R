#!/usr/bin/env Rscript
#
# 12-clean-mar-dgp
#
# Cleanly MAR replacement candidates for the old DGP C. The baseline
# mechanism always observes rater 1 and lets the remaining raters' observation
# probabilities depend on rater 1's observed category. Candidate search also
# includes sequential MAR mechanisms where rater j is observed according to
# the observed value/status of rater j - 1.
#
# Outputs (under results/):
#   baseline_summary.csv      n-sweep for the moderate anchor-MAR DGP
#   baseline_replicates.csv   per-replicate estimates for the baseline
#   candidate_summary.csv     candidate-by-method summaries at one n
#   candidate_ranking.csv     compact ranking of dramatic candidates
#   candidate_replicates.csv  per-replicate estimates for the search
#   mechanisms.csv            DGP/mechanism parameter table
#   metadata.csv              run metadata

suppressPackageStartupMessages({
  library(misskappa)
})

# ---- CLI ----------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
has_flag <- function(name) name %in% args
get_val <- function(name, default, parser = as.character) {
  i <- which(args == name)
  if (length(i) == 0L) return(default)
  parser(args[i + 1L])
}
parse_int_csv <- function(x) as.integer(strsplit(x, ",", fixed = TRUE)[[1L]])

if (has_flag("--help") || has_flag("-h")) {
  cat(
    "Usage: Rscript run_experiment.R [options]\n",
    " --smoke             Fast run (small reps/n_truth).\n",
    " --baseline-reps N   Replicates for baseline n-sweep (default 120).\n",
    " --search-reps N     Replicates for candidate search (default 80).\n",
    " --baseline-ns CSV   Baseline n grid (default 500,2000,8000).\n",
    " --search-n N        Candidate-search sample size (default 4000).\n",
    " --truth-n N         Complete-data truth sample size (default 200000).\n",
    " --seed-base K       Seed base (default 20260527).\n",
    " --help, -h          This help.\n",
    sep = ""
  )
  quit("no", status = 0)
}

smoke <- has_flag("--smoke")
seed_base <- get_val("--seed-base", 20260527L, as.integer)

baseline_reps <- if (smoke) 8L else 120L
search_reps <- if (smoke) 8L else 80L
baseline_ns <- if (smoke) c(300L, 1000L) else c(500L, 2000L, 8000L)
search_n <- if (smoke) 1000L else 4000L
truth_n <- if (smoke) 30000L else 200000L

baseline_reps <- get_val("--baseline-reps", baseline_reps, as.integer)
search_reps <- get_val("--search-reps", search_reps, as.integer)
baseline_ns <- get_val("--baseline-ns", baseline_ns, parse_int_csv)
search_n <- get_val("--search-n", search_n, as.integer)
truth_n <- get_val("--truth-n", truth_n, as.integer)

methods <- c("available", "ipw", "fiml", "gwet")

# ---- DGP specifications -------------------------------------------------
make_spec <- function(id, name, mechanism, C, R, p, rho_base, rho_truth_mult,
                      pi_anchor, pi_fallback = NA_real_) {
  list(
    id = id,
    name = name,
    mechanism = mechanism,
    C = C,
    R = R,
    p = p,
    rho_base = rho_base,
    rho_truth_mult = rho_truth_mult,
    pi_anchor = pi_anchor,
    pi_fallback = pi_fallback
  )
}

specs <- list(
  make_spec(
    1L, "anchor_moderate", "anchor", 3L, 4L,
    p = c(0.50, 0.30, 0.20),
    rho_base = c(0.94, 0.90, 0.84, 0.78),
    rho_truth_mult = c(1.05, 0.95, 0.75),
    pi_anchor = c(0.95, 0.60, 0.25)
  ),
  make_spec(
    2L, "anchor_severe", "anchor", 3L, 4L,
    p = c(0.45, 0.30, 0.25),
    rho_base = c(0.96, 0.90, 0.82, 0.74),
    rho_truth_mult = c(1.10, 0.90, 0.45),
    pi_anchor = c(0.98, 0.45, 0.06)
  ),
  make_spec(
    3L, "anchor_extreme_skew", "anchor", 3L, 4L,
    p = c(0.65, 0.25, 0.10),
    rho_base = c(0.97, 0.88, 0.78, 0.68),
    rho_truth_mult = c(1.05, 0.75, 0.30),
    pi_anchor = c(0.99, 0.35, 0.04)
  ),
  make_spec(
    4L, "sequential_severe", "sequential", 3L, 4L,
    p = c(0.45, 0.30, 0.25),
    rho_base = c(0.96, 0.90, 0.82, 0.74),
    rho_truth_mult = c(1.10, 0.90, 0.45),
    pi_anchor = c(0.98, 0.45, 0.06),
    pi_fallback = 0.12
  ),
  make_spec(
    5L, "anchor_four_cat", "anchor", 4L, 4L,
    p = c(0.40, 0.30, 0.20, 0.10),
    rho_base = c(0.96, 0.89, 0.80, 0.72),
    rho_truth_mult = c(1.10, 0.95, 0.65, 0.25),
    pi_anchor = c(0.99, 0.65, 0.22, 0.04)
  ),
  make_spec(
    6L, "sequential_four_cat", "sequential", 4L, 4L,
    p = c(0.40, 0.30, 0.20, 0.10),
    rho_base = c(0.96, 0.89, 0.80, 0.72),
    rho_truth_mult = c(1.10, 0.95, 0.65, 0.25),
    pi_anchor = c(0.99, 0.65, 0.22, 0.04),
    pi_fallback = 0.10
  )
)
names(specs) <- vapply(specs, `[[`, character(1), "name")

baseline_spec <- specs[["anchor_moderate"]]

# ---- Simulation helpers -------------------------------------------------
simulate_ratings <- function(n, spec) {
  C <- spec$C
  R <- spec$R
  truth <- sample.int(C, n, replace = TRUE, prob = spec$p)
  guess_mat <- matrix(rep(spec$p, each = R), nrow = R, byrow = TRUE)
  x_star <- matrix(NA_integer_, nrow = n, ncol = R)
  for (j in seq_len(R)) {
    rho_i <- pmin(pmax(spec$rho_base[j] * spec$rho_truth_mult[truth], 0), 1)
    correct <- stats::rbinom(n, 1L, rho_i) == 1L
    guessed <- sample.int(C, n, replace = TRUE, prob = guess_mat[j, ])
    x_star[, j] <- ifelse(correct, truth, guessed)
  }
  x_star
}

apply_missing_clean_mar <- function(x_star, spec) {
  R <- spec$R
  x <- x_star
  if (R == 1L) return(x)

  if (identical(spec$mechanism, "anchor")) {
    for (j in 2:R) {
      observed <- stats::runif(nrow(x_star)) < spec$pi_anchor[x_star[, 1L]]
      x[!observed, j] <- NA_integer_
    }
    return(x)
  }

  if (identical(spec$mechanism, "sequential")) {
    for (j in 2:R) {
      prev <- x[, j - 1L]
      p_obs <- ifelse(is.na(prev), spec$pi_fallback, spec$pi_anchor[prev])
      observed <- stats::runif(nrow(x_star)) < p_obs
      x[!observed, j] <- NA_integer_
    }
    return(x)
  }

  stop("Unknown missingness mechanism: ", spec$mechanism, call. = FALSE)
}

fit_one <- function(x, method) {
  args <- list(x = x, method = method, weight = "identity")
  if (identical(method, "fiml")) {
    args$em_options <- list(max_iter = 50000L, tol = 1e-8)
  }
  res <- try(do.call(misskappa::kappa, args), silent = TRUE)
  if (inherits(res, "try-error")) return(NA_real_)
  as.numeric(res$estimates[["Conger"]])
}

truth_for_spec <- function(spec) {
  set.seed(seed_base + spec$id * 100000L + 999L)
  x_star <- simulate_ratings(truth_n, spec)
  k <- misskappa::kappa(x_star, method = "available", weight = "identity")
  as.numeric(k$estimates[["Conger"]])
}

make_seed <- function(spec, n, b, phase_offset) {
  as.integer(seed_base) +
    as.integer(phase_offset) +
    as.integer(spec$id) * 100000L +
    as.integer(n) * 10L +
    as.integer(b)
}

run_cells <- function(spec, ns, reps, phase_offset) {
  truth <- truth_for_spec(spec)
  summary_rows <- list()
  rep_rows <- list()
  for (n in ns) {
    vals <- matrix(NA_real_, nrow = reps, ncol = length(methods),
                   dimnames = list(NULL, methods))
    obs_frac <- numeric(reps)
    for (b in seq_len(reps)) {
      set.seed(make_seed(spec, n, b, phase_offset))
      x_star <- simulate_ratings(n, spec)
      x <- apply_missing_clean_mar(x_star, spec)
      obs_frac[b] <- mean(!is.na(x))
      for (m in methods) vals[b, m] <- fit_one(x, m)
    }
    for (m in methods) {
      est <- vals[, m]
      ok <- is.finite(est)
      n_ok <- sum(ok)
      bias <- mean(est[ok]) - truth
      sd_est <- stats::sd(est[ok])
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(
        spec = spec$name,
        mechanism = spec$mechanism,
        C = spec$C,
        R = spec$R,
        n = n,
        reps = reps,
        n_ok = n_ok,
        method = m,
        truth = truth,
        mean_est = mean(est[ok]),
        bias = bias,
        abs_bias = abs(bias),
        sd_est = sd_est,
        mse = bias^2 + sd_est^2,
        mc_se_bias = sd_est / sqrt(n_ok),
        mean_observed_fraction = mean(obs_frac)
      )
      rep_rows[[length(rep_rows) + 1L]] <- data.frame(
        spec = spec$name,
        mechanism = spec$mechanism,
        n = n,
        b = seq_len(reps),
        method = m,
        truth = truth,
        est = est,
        observed_fraction = obs_frac
      )
    }
    cat(sprintf("[%s] %s n=%d reps=%d -- done\n",
                format(Sys.time(), "%H:%M:%S"), spec$name, n, reps))
  }
  list(
    summary = do.call(rbind, summary_rows),
    replicates = do.call(rbind, rep_rows)
  )
}

mechanism_rows <- do.call(rbind, lapply(specs, function(spec) {
  data.frame(
    spec = spec$name,
    id = spec$id,
    mechanism = spec$mechanism,
    C = spec$C,
    R = spec$R,
    p = paste(spec$p, collapse = ","),
    rho_base = paste(spec$rho_base, collapse = ","),
    rho_truth_mult = paste(spec$rho_truth_mult, collapse = ","),
    pi_anchor = paste(spec$pi_anchor, collapse = ","),
    pi_fallback = ifelse(is.finite(spec$pi_fallback), spec$pi_fallback, NA_real_),
    stringsAsFactors = FALSE
  )
}))

# ---- Run ----------------------------------------------------------------
t0 <- Sys.time()
baseline <- run_cells(baseline_spec, baseline_ns, baseline_reps, phase_offset = 0L)
candidate_pieces <- lapply(specs, run_cells, ns = search_n, reps = search_reps,
                           phase_offset = 50000000L)
candidate_summary <- do.call(rbind, lapply(candidate_pieces, `[[`, "summary"))
candidate_replicates <- do.call(rbind, lapply(candidate_pieces, `[[`, "replicates"))

candidate_ranking <- do.call(rbind, lapply(split(candidate_summary, candidate_summary$spec),
                                           function(g) {
  fiml <- g[g$method == "fiml", ]
  other <- g[g$method != "fiml", ]
  best_other <- other[which.max(other$abs_bias), ]
  data.frame(
    spec = unique(g$spec),
    mechanism = unique(g$mechanism),
    C = unique(g$C),
    R = unique(g$R),
    n = unique(g$n),
    reps = unique(g$reps),
    mean_observed_fraction = unique(g$mean_observed_fraction),
    truth = unique(g$truth),
    largest_non_fiml_method = best_other$method,
    largest_non_fiml_bias = best_other$bias,
    largest_non_fiml_abs_bias = best_other$abs_bias,
    fiml_bias = fiml$bias,
    fiml_abs_bias = fiml$abs_bias,
    abs_bias_gap = best_other$abs_bias - fiml$abs_bias,
    stringsAsFactors = FALSE
  )
}))
candidate_ranking <- candidate_ranking[order(-candidate_ranking$abs_bias_gap), ]

# ---- Write --------------------------------------------------------------
results_dir <- "results"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(baseline$summary, file.path(results_dir, "baseline_summary.csv"),
          row.names = FALSE)
write.csv(baseline$replicates, file.path(results_dir, "baseline_replicates.csv"),
          row.names = FALSE)
write.csv(candidate_summary, file.path(results_dir, "candidate_summary.csv"),
          row.names = FALSE)
write.csv(candidate_replicates, file.path(results_dir, "candidate_replicates.csv"),
          row.names = FALSE)
write.csv(candidate_ranking, file.path(results_dir, "candidate_ranking.csv"),
          row.names = FALSE)
write.csv(mechanism_rows, file.path(results_dir, "mechanisms.csv"),
          row.names = FALSE)

meta <- data.frame(
  key = c("seed_base", "smoke", "baseline_reps", "search_reps",
          "baseline_ns", "search_n", "truth_n", "methods", "weight",
          "kappa", "R_version", "misskappa_version", "started_at",
          "elapsed_s"),
  value = c(
    as.character(seed_base),
    as.character(smoke),
    as.character(baseline_reps),
    as.character(search_reps),
    paste(baseline_ns, collapse = ","),
    as.character(search_n),
    as.character(truth_n),
    paste(methods, collapse = ","),
    "identity",
    "Conger",
    as.character(getRversion()),
    as.character(utils::packageVersion("misskappa")),
    format(t0, "%Y-%m-%dT%H:%M:%S"),
    sprintf("%.1f", as.numeric(difftime(Sys.time(), t0, units = "secs")))
  ),
  stringsAsFactors = FALSE
)
write.csv(meta, file.path(results_dir, "metadata.csv"), row.names = FALSE)

cat(sprintf("Wrote %s/baseline_summary.csv (%d rows)\n",
            results_dir, nrow(baseline$summary)))
cat(sprintf("Wrote %s/candidate_summary.csv (%d rows)\n",
            results_dir, nrow(candidate_summary)))
cat(sprintf("Wrote %s/candidate_ranking.csv (%d rows)\n",
            results_dir, nrow(candidate_ranking)))
cat(sprintf("Wrote %s/metadata.csv\n", results_dir))
