#!/usr/bin/env Rscript
#
# Experiment 27: finite-sample MSE and Wald coverage for the current kappa
# estimator surface. The design deliberately mirrors the alpha-missing paper
# runner: several missingness mechanisms, a replicate-level CSV, and summaries
# for bias, MSE, SE calibration, coverage, failures, missingness diagnostics,
# and timing.
#
# The public estimator grid is:
#   ipw / cat_fiml       x nominal, absolute(linear), quadratic
#   pairwise / nt_fiml   x quadratic only
#
# The optional cat_fiml_jk5_* and cat_fiml_jk10_* rows are experiment-side
# grouped-jackknife diagnostics. They refit Cat-FIML after deleting
# deterministic folds and replace the point estimate by the delete-group
# jackknife bias correction.

suppressPackageStartupMessages({
  library(misskappa)
})

usage <- function(status = 0L) {
  cat(paste0(
    "Usage: Rscript run_experiment.R [mode/options]\n\n",
    "Modes:\n",
    "  --smoke            Minimal mechanical check: n=80, reps=1.\n",
    "  --small            Cheap pilot: n=120, reps=3. Default mode.\n",
    "  --big              Paper-grade starting grid: n=200,600, reps=500.\n\n",
    "Options:\n",
    "  --reps N            Override replications per cell.\n",
    "  --n-grid CSV        Override sample sizes.\n",
    "  --truth-n N         Complete-data Monte Carlo truth size.\n",
    "  --dgps CSV          DGP ids.\n",
    "  --mechanisms CSV    Missingness mechanism ids.\n",
    "  --methods CSV       Method ids.\n",
    "  --seed-base N       Base seed. Default: 272700.\n",
    "  --out-dir PATH      Output directory. Default: script-local results/.\n",
    "  --resume            Reuse existing per-cell checkpoints in --out-dir.\n",
    "  --progress          Print one line per design cell.\n",
    "  --help, -h          Show this help and exit.\n\n",
    "Examples:\n",
    "  Rscript run_experiment.R --smoke --progress\n",
    "  Rscript run_experiment.R --small --progress\n",
    "  Rscript run_experiment.R --big --reps 1000 --resume --progress\n",
    "  Rscript run_experiment.R --big --methods cat_fiml_quadratic,cat_fiml_jk5_quadratic --progress\n"
  ))
  quit(save = "no", status = status)
}

script_dir <- local({
  cmd <- commandArgs(FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) == 0L) getwd()
  else dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE))
})

parse_csv_chr <- function(x) trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
parse_csv_int <- function(x, arg) {
  out <- as.integer(strsplit(x, ",", fixed = TRUE)[[1L]])
  if (any(is.na(out))) stop(arg, " must be a comma-separated integer list.", call. = FALSE)
  unique(out)
}

mode_defaults <- function(mode) {
  if (mode == "smoke") {
    return(list(
      mode = mode,
      reps = 1L,
      n_grid = 80L,
      truth_n = 5000L,
      dgps = c("balanced4"),
      mechanisms = c("mcar30", "mar_anchor30"),
      methods = c("ipw_nominal", "cat_fiml_nominal",
                  "ipw_absolute", "cat_fiml_absolute",
                  "ipw_quadratic", "cat_fiml_quadratic",
                  "pairwise_quadratic", "nt_fiml_quadratic")
    ))
  }
  if (mode == "small") {
    return(list(
      mode = mode,
      reps = 3L,
      n_grid = 120L,
      truth_n = 15000L,
      dgps = c("balanced4", "biased4"),
      mechanisms = c("mcar30", "mcar_rater30", "mcar_subject30",
                     "mar_anchor15", "mar_anchor30",
                     "mar_anchor_nonlinear30", "mar_sequential30"),
      methods = c("ipw_nominal", "cat_fiml_nominal",
                  "ipw_absolute", "cat_fiml_absolute",
                  "ipw_quadratic", "cat_fiml_quadratic",
                  "pairwise_quadratic", "nt_fiml_quadratic")
    ))
  }
  if (mode == "big") {
    return(list(
      mode = mode,
      reps = 500L,
      n_grid = c(200L, 600L),
      truth_n = 300000L,
      dgps = c("balanced4", "biased4", "sparse5"),
      mechanisms = c("complete", "mcar15", "mcar30", "mcar_rater30",
                     "mcar_subject30", "mar_anchor15", "mar_anchor30",
                     "mar_anchor_nonlinear30", "mar_shifted30",
                     "mar_sequential30"),
      methods = c("ipw_nominal", "cat_fiml_nominal",
                  "ipw_absolute", "cat_fiml_absolute",
                  "ipw_quadratic", "cat_fiml_quadratic",
                  "pairwise_quadratic", "nt_fiml_quadratic")
    ))
  }
  stop("Unknown mode: ", mode, call. = FALSE)
}

parse_args <- function(argv) {
  mode <- "small"
  if ("--smoke" %in% argv) mode <- "smoke"
  if ("--small" %in% argv) mode <- "small"
  if ("--big" %in% argv) mode <- "big"
  opts <- mode_defaults(mode)
  opts$seed_base <- 272700L
  opts$out_dir <- file.path(script_dir, "results")
  opts$progress <- FALSE
  opts$resume <- FALSE

  i <- 1L
  while (i <= length(argv)) {
    arg <- argv[[i]]
    if (arg %in% c("--help", "-h")) usage(0L)
    if (arg %in% c("--smoke", "--small", "--big")) {
      i <- i + 1L
      next
    }
    if (arg == "--progress") {
      opts$progress <- TRUE
      i <- i + 1L
      next
    }
    if (arg == "--resume") {
      opts$resume <- TRUE
      i <- i + 1L
      next
    }
    needs_value <- c("--reps", "--n-grid", "--truth-n", "--dgps",
                     "--mechanisms", "--methods", "--seed-base", "--out-dir")
    if (arg %in% needs_value) {
      if (i == length(argv)) stop(arg, " needs a value.", call. = FALSE)
      val <- argv[[i + 1L]]
      if (arg == "--reps") opts$reps <- as.integer(val)
      else if (arg == "--n-grid") opts$n_grid <- parse_csv_int(val, arg)
      else if (arg == "--truth-n") opts$truth_n <- as.integer(val)
      else if (arg == "--dgps") opts$dgps <- parse_csv_chr(val)
      else if (arg == "--mechanisms") opts$mechanisms <- parse_csv_chr(val)
      else if (arg == "--methods") opts$methods <- parse_csv_chr(val)
      else if (arg == "--seed-base") opts$seed_base <- as.integer(val)
      else if (arg == "--out-dir") opts$out_dir <- val
      i <- i + 2L
      next
    }
    stop("Unknown argument: ", arg, call. = FALSE)
  }

  if (is.na(opts$reps) || opts$reps < 1L) stop("--reps must be >= 1.", call. = FALSE)
  if (any(is.na(opts$n_grid)) || any(opts$n_grid < 20L)) {
    stop("--n-grid entries must be integers >= 20.", call. = FALSE)
  }
  if (is.na(opts$truth_n) || opts$truth_n < 1000L) {
    stop("--truth-n must be an integer >= 1000.", call. = FALSE)
  }
  if (is.na(opts$seed_base)) stop("--seed-base must be an integer.", call. = FALSE)
  opts
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))

# ---- DGPs -------------------------------------------------------------------

centered_guess <- function(C, centers, spread = 0.9) {
  vals <- seq_len(C)
  G <- vapply(centers, function(ctr) {
    p <- exp(-0.5 * ((vals - ctr) / spread)^2)
    p / sum(p)
  }, numeric(C))
  t(G)
}

dgp_defs <- local({
  balanced4 <- list(
    id = "balanced4",
    label = "Balanced four-category, four-rater skill model",
    C = 4L,
    R = 4L,
    p = rep(0.25, 4L),
    skill = rep(0.72, 4L),
    guess = matrix(rep(0.25, 16L), nrow = 4L, byrow = TRUE),
    truth_skill_mult = rep(1, 4L)
  )
  biased4 <- list(
    id = "biased4",
    label = "Biased four-category, non-exchangeable raters",
    C = 4L,
    R = 4L,
    p = c(0.45, 0.30, 0.18, 0.07),
    skill = c(0.88, 0.78, 0.66, 0.58),
    guess = centered_guess(4L, c(1.2, 1.8, 3.0, 3.7), spread = 0.75),
    truth_skill_mult = c(1.05, 0.98, 0.86, 0.70)
  )
  sparse5 <- list(
    id = "sparse5",
    label = "Sparse five-category, five-rater high-agreement model",
    C = 5L,
    R = 5L,
    p = c(0.55, 0.22, 0.12, 0.07, 0.04),
    skill = c(0.90, 0.86, 0.82, 0.78, 0.74),
    guess = centered_guess(5L, c(1.1, 1.7, 2.6, 3.7, 4.5), spread = 0.8),
    truth_skill_mult = c(1.05, 0.98, 0.88, 0.74, 0.58)
  )
  list(balanced4 = balanced4, biased4 = biased4, sparse5 = sparse5)
})

simulate_ratings <- function(n, spec) {
  truth <- sample.int(spec$C, n, replace = TRUE, prob = spec$p)
  X <- matrix(NA_integer_, nrow = n, ncol = spec$R)
  for (j in seq_len(spec$R)) {
    p_correct <- pmin(pmax(spec$skill[j] * spec$truth_skill_mult[truth], 0), 1)
    correct <- stats::rbinom(n, 1L, p_correct) == 1L
    guessed <- sample.int(spec$C, n, replace = TRUE, prob = spec$guess[j, ])
    X[, j] <- ifelse(correct, truth, guessed)
  }
  colnames(X) <- paste0("r", seq_len(spec$R))
  X
}

category_z <- function(x, C) {
  z <- (as.numeric(x) - mean(seq_len(C))) / stats::sd(seq_len(C))
  z[!is.finite(z)] <- 0
  z
}

tune_intercept <- function(eta, target) {
  f <- function(a) mean(stats::plogis(a + eta)) - target
  stats::uniroot(f, c(-20, 20))$root
}

drop_with_prob <- function(X, P, keep_first = FALSE) {
  out <- X
  M <- matrix(stats::runif(length(X)), nrow = nrow(X), ncol = ncol(X)) < P
  if (keep_first && ncol(out) >= 1L) M[, 1L] <- FALSE
  out[M] <- NA_integer_
  out
}

apply_mcar <- function(X, rate) {
  drop_with_prob(X, matrix(rate, nrow(X), ncol(X)))
}

apply_mcar_rater <- function(X, target = 0.30) {
  R <- ncol(X)
  rates <- seq(max(0.04, target - 0.22), min(0.70, target + 0.22), length.out = R)
  rates <- pmin(pmax(rates - mean(rates) + target, 0.02), 0.80)
  P <- matrix(rates, nrow(X), R, byrow = TRUE)
  drop_with_prob(X, P)
}

apply_mcar_subject <- function(X, target = 0.30) {
  n <- nrow(X)
  p_i <- stats::rgamma(n, shape = 2.5, rate = 2.5)
  p_i <- pmin(0.85, p_i * target / mean(p_i))
  P <- matrix(p_i, n, ncol(X))
  drop_with_prob(X, P)
}

apply_anchor_mar <- function(X, C, target = 0.30, slope = 1.15,
                             nonlinear = FALSE, shifted = FALSE) {
  n <- nrow(X); R <- ncol(X)
  anchor <- category_z(X[, 1L], C)
  eta_base <- if (nonlinear) {
    centered <- anchor^2 - mean(anchor^2)
    1.35 * centered
  } else {
    slope * anchor
  }
  P <- matrix(0, n, R)
  for (j in seq_len(R)) {
    if (j == 1L) next
    shift <- if (shifted) (j - (R + 1) / 2) * 0.35 else 0
    eta <- eta_base + shift
    a <- tune_intercept(eta, target)
    P[, j] <- stats::plogis(a + eta)
  }
  drop_with_prob(X, P, keep_first = TRUE)
}

apply_sequential_mar <- function(X, C, target = 0.30) {
  out <- X
  R <- ncol(out)
  for (j in 2:R) {
    prev <- out[, j - 1L]
    prev_missing <- is.na(prev)
    z <- category_z(prev, C)
    eta <- 1.0 * z + ifelse(prev_missing, 1.25, 0)
    a <- tune_intercept(eta, target)
    p_miss <- stats::plogis(a + eta)
    out[stats::runif(nrow(out)) < p_miss, j] <- NA_integer_
  }
  out
}

apply_missing <- function(X, spec, mechanism) {
  if (mechanism == "complete") return(X)
  if (mechanism == "mcar15") return(apply_mcar(X, 0.15))
  if (mechanism == "mcar30") return(apply_mcar(X, 0.30))
  if (mechanism == "mcar_rater30") return(apply_mcar_rater(X, 0.30))
  if (mechanism == "mcar_subject30") return(apply_mcar_subject(X, 0.30))
  if (mechanism == "mar_anchor15") return(apply_anchor_mar(X, spec$C, 0.15))
  if (mechanism == "mar_anchor30") return(apply_anchor_mar(X, spec$C, 0.30))
  if (mechanism == "mar_anchor_nonlinear30") {
    return(apply_anchor_mar(X, spec$C, 0.30, nonlinear = TRUE))
  }
  if (mechanism == "mar_shifted30") {
    return(apply_anchor_mar(X, spec$C, 0.30, shifted = TRUE))
  }
  if (mechanism == "mar_sequential30") return(apply_sequential_mar(X, spec$C, 0.30))
  stop("Unknown mechanism: ", mechanism, call. = FALSE)
}

# ---- estimators -------------------------------------------------------------

method_defs <- data.frame(
  method = c("ipw_nominal", "cat_fiml_nominal",
             "cat_fiml_jk5_nominal", "cat_fiml_jk10_nominal",
             "ipw_absolute", "cat_fiml_absolute",
             "cat_fiml_jk5_absolute", "cat_fiml_jk10_absolute",
             "ipw_quadratic", "cat_fiml_quadratic",
             "cat_fiml_jk5_quadratic", "cat_fiml_jk10_quadratic",
             "pairwise_quadratic", "nt_fiml_quadratic"),
  estimator = c("ipw", "cat_fiml", "cat_fiml_jk5", "cat_fiml_jk10",
                "ipw", "cat_fiml", "cat_fiml_jk5", "cat_fiml_jk10",
                "ipw", "cat_fiml", "cat_fiml_jk5", "cat_fiml_jk10",
                "pairwise", "nt_fiml"),
  base_estimator = c("ipw", "cat_fiml", "cat_fiml", "cat_fiml",
                     "ipw", "cat_fiml", "cat_fiml", "cat_fiml",
                     "ipw", "cat_fiml", "cat_fiml", "cat_fiml",
                     "pairwise", "nt_fiml"),
  weight = c("nominal", "nominal", "nominal", "nominal",
             "linear", "linear", "linear", "linear",
             "quadratic", "quadratic", "quadratic", "quadratic",
             "quadratic", "quadratic"),
  weight_label = c("nominal", "nominal", "nominal", "nominal",
                   "absolute", "absolute", "absolute", "absolute",
                   "quadratic", "quadratic", "quadratic", "quadratic",
                   "quadratic", "quadratic"),
  jackknife_groups = c(0L, 0L, 5L, 10L, 0L, 0L, 5L, 10L,
                       0L, 0L, 5L, 10L, 0L, 0L),
  stringsAsFactors = FALSE
)

valid_dgps <- names(dgp_defs)
valid_mechanisms <- c("complete", "mcar15", "mcar30", "mcar_rater30",
                      "mcar_subject30", "mar_anchor15", "mar_anchor30",
                      "mar_anchor_nonlinear30", "mar_shifted30",
                      "mar_sequential30")
valid_methods <- method_defs$method

if (!all(opts$dgps %in% valid_dgps)) {
  stop("--dgps must contain only: ", paste(valid_dgps, collapse = ","), call. = FALSE)
}
if (!all(opts$mechanisms %in% valid_mechanisms)) {
  stop("--mechanisms must contain only: ", paste(valid_mechanisms, collapse = ","), call. = FALSE)
}
if (!all(opts$methods %in% valid_methods)) {
  stop("--methods must contain only: ", paste(valid_methods, collapse = ","), call. = FALSE)
}
method_defs <- method_defs[match(opts$methods, method_defs$method), , drop = FALSE]

coef_names <- c("Conger", "Fleiss")
z975 <- stats::qnorm(0.975)

se_vector <- function(fit, coefs) {
  V <- stats::vcov(fit)
  out <- rep(NA_real_, length(coefs)); names(out) <- coefs
  have <- intersect(coefs, rownames(V))
  vals <- diag(V)[have]
  vals[!is.finite(vals) | vals < 0] <- NA_real_
  out[have] <- sqrt(vals)
  out
}

fisher_interval <- function(est, se) {
  if (!is.finite(est) || !is.finite(se) || se < 0 || abs(est) >= 1) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  se_z <- se / (1 - est^2)
  g <- atanh(est)
  c(lower = tanh(g - z975 * se_z), upper = tanh(g + z975 * se_z))
}

em_options_for <- function(estimator) {
  if (estimator == "cat_fiml") {
    return(list(tol = 1e-7, max_iter = 12000L, prune_tol = 1e-10,
                start_alpha = 0.1, info_rcond = 1e-4))
  }
  if (estimator == "nt_fiml") {
    return(list(tol = 1e-8, max_iter = 12000L, fd_h = 1e-5))
  }
  list()
}

fit_kappa <- function(X, spec, estimator, weight) {
  misskappa::kappa(
    X,
    estimator = estimator,
    weight = weight,
    values = seq(0, spec$C - 1L),
    em_options = em_options_for(estimator)
  )
}

fit_one <- function(X, spec, method_row) {
  start <- proc.time()[["elapsed"]]
  out <- tryCatch({
    fit <- fit_kappa(X, spec, method_row$base_estimator, method_row$weight)
    est_raw <- stats::coef(fit)
    se <- se_vector(fit, coef_names)

    est <- est_raw
    jk_bias <- rep(NA_real_, length(coef_names)); names(jk_bias) <- coef_names
    jk_refits <- 0L
    jk_groups <- as.integer(method_row$jackknife_groups)
    if (jk_groups > 0L) {
      jk_groups <- min(jk_groups, nrow(X) - 1L)
      if (jk_groups < 2L) stop("Grouped jackknife needs at least two folds.", call. = FALSE)

      fold_id <- ((seq_len(nrow(X)) - 1L) %% jk_groups) + 1L
      fold_rows <- split(seq_len(nrow(X)), fold_id)
      delete_est <- matrix(
        NA_real_,
        nrow = length(fold_rows),
        ncol = length(coef_names),
        dimnames = list(NULL, coef_names)
      )
      for (g in seq_along(fold_rows)) {
        fit_g <- fit_kappa(X[-fold_rows[[g]], , drop = FALSE],
                           spec, method_row$base_estimator, method_row$weight)
        delete_est[g, ] <- stats::coef(fit_g)[coef_names]
      }
      jk_refits <- length(fold_rows)
      delete_mean <- colMeans(delete_est)
      jk_bias <- (jk_refits - 1) * (delete_mean - est_raw[coef_names])
      est[coef_names] <- est_raw[coef_names] - jk_bias
    }

    rows <- lapply(coef_names, function(coef) {
      fci <- fisher_interval(est[[coef]], se[[coef]])
      data.frame(
        coefficient = coef,
        estimate = unname(est[[coef]]),
        estimate_raw = unname(est_raw[[coef]]),
        jackknife_bias = unname(jk_bias[[coef]]),
        jackknife_groups = jk_groups,
        jackknife_refits = jk_refits,
        se = unname(se[[coef]]),
        lower = unname(est[[coef]] - z975 * se[[coef]]),
        upper = unname(est[[coef]] + z975 * se[[coef]]),
        fisher_lower = fci[["lower"]],
        fisher_upper = fci[["upper"]],
        error = "",
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, rows)
  }, error = function(e) {
    data.frame(
      coefficient = coef_names,
      estimate = NA_real_,
      estimate_raw = NA_real_,
      jackknife_bias = NA_real_,
      jackknife_groups = as.integer(method_row$jackknife_groups),
      jackknife_refits = NA_integer_,
      se = NA_real_,
      lower = NA_real_,
      upper = NA_real_,
      fisher_lower = NA_real_,
      fisher_upper = NA_real_,
      error = conditionMessage(e),
      stringsAsFactors = FALSE
    )
  })
  out$elapsed_ms <- 1000 * (proc.time()[["elapsed"]] - start)
  out
}

min_pair_count <- function(X) {
  R <- ncol(X)
  counts <- integer(R * (R - 1L) / 2L)
  idx <- 0L
  for (j in seq_len(R - 1L)) {
    for (k in (j + 1L):R) {
      idx <- idx + 1L
      counts[idx] <- sum(!is.na(X[, j]) & !is.na(X[, k]))
    }
  }
  min(counts)
}

truth_for_dgp <- function(spec, seed, truth_n) {
  set.seed(seed)
  X <- simulate_ratings(truth_n, spec)
  rows <- list()
  weights <- unique(method_defs[, c("weight", "weight_label")])
  for (i in seq_len(nrow(weights))) {
    w <- weights$weight[i]
    wl <- weights$weight_label[i]
    fit <- misskappa::kappa(
      X,
      estimator = "ipw",
      weight = w,
      values = seq(0, spec$C - 1L)
    )
    est <- stats::coef(fit)
    rows[[length(rows) + 1L]] <- data.frame(
      dgp = spec$id,
      dgp_label = spec$label,
      C = spec$C,
      R = spec$R,
      weight = w,
      weight_label = wl,
      coefficient = coef_names,
      truth = unname(est[coef_names]),
      truth_n = truth_n,
      truth_method = "complete-data Monte Carlo using kappa(estimator='ipw')",
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

split_keys <- function(data, keys) {
  interaction(data[, keys], drop = TRUE, lex.order = TRUE)
}

summarize_replicates <- function(df) {
  keys <- c("dgp", "dgp_label", "C", "R", "mechanism", "n", "method",
            "estimator", "weight_label", "coefficient")
  pieces <- lapply(split(df, split_keys(df, keys)), function(g) {
    ok <- is.finite(g$estimate)
    se_ok <- ok & is.finite(g$se) & g$se > 0
    raw_est <- if ("estimate_raw" %in% names(g)) g$estimate_raw else g$estimate
    raw_ok <- is.finite(raw_est)
    jk_bias <- if ("jackknife_bias" %in% names(g)) g$jackknife_bias else rep(NA_real_, nrow(g))
    jk_refits <- if ("jackknife_refits" %in% names(g)) g$jackknife_refits else rep(NA_integer_, nrow(g))
    jk_groups <- if ("jackknife_groups" %in% names(g)) g$jackknife_groups else rep(0L, nrow(g))
    err <- g$estimate - g$truth
    raw_err <- raw_est - g$truth
    nat_cover <- se_ok & is.finite(g$lower) & is.finite(g$upper) &
      g$lower <= g$truth & g$truth <= g$upper
    fish_cover <- se_ok & is.finite(g$fisher_lower) & is.finite(g$fisher_upper) &
      g$fisher_lower <= g$truth & g$truth <= g$fisher_upper
    data.frame(
      dgp = g$dgp[[1L]],
      dgp_label = g$dgp_label[[1L]],
      C = g$C[[1L]],
      R = g$R[[1L]],
      mechanism = g$mechanism[[1L]],
      n = g$n[[1L]],
      method = g$method[[1L]],
      estimator = g$estimator[[1L]],
      weight_label = g$weight_label[[1L]],
      coefficient = g$coefficient[[1L]],
      reps = length(unique(g$rep)),
      n_valid = sum(ok),
      failures = sum(!ok),
      truth = g$truth[[1L]],
      mean_estimate = if (any(ok)) mean(g$estimate[ok]) else NA_real_,
      mean_raw_estimate = if (any(raw_ok)) mean(raw_est[raw_ok]) else NA_real_,
      bias = if (any(ok)) mean(err[ok]) else NA_real_,
      raw_bias = if (any(raw_ok)) mean(raw_err[raw_ok]) else NA_real_,
      mean_jackknife_bias = if (any(is.finite(jk_bias))) mean(jk_bias, na.rm = TRUE) else NA_real_,
      sd_estimate = if (sum(ok) > 1L) stats::sd(g$estimate[ok]) else NA_real_,
      mc_se_bias = if (sum(ok) > 1L) stats::sd(g$estimate[ok]) / sqrt(sum(ok)) else NA_real_,
      mse = if (any(ok)) mean(err[ok]^2) else NA_real_,
      rmse = if (any(ok)) sqrt(mean(err[ok]^2)) else NA_real_,
      mean_se = if (any(se_ok)) mean(g$se[se_ok]) else NA_real_,
      se_over_sd = if (sum(se_ok) > 1L && stats::sd(g$estimate[se_ok]) > 1e-12) {
        mean(g$se[se_ok]) / stats::sd(g$estimate[se_ok])
      } else NA_real_,
      coverage95 = if (any(se_ok)) mean(nat_cover[se_ok]) else NA_real_,
      fisher_coverage95 = if (any(se_ok)) mean(fish_cover[se_ok]) else NA_real_,
      mean_ci_length = if (any(se_ok)) mean(g$upper[se_ok] - g$lower[se_ok]) else NA_real_,
      mean_fisher_length = if (any(se_ok)) {
        mean(g$fisher_upper[se_ok] - g$fisher_lower[se_ok], na.rm = TRUE)
      } else NA_real_,
      mean_elapsed_ms = mean(g$elapsed_ms, na.rm = TRUE),
      median_elapsed_ms = stats::median(g$elapsed_ms, na.rm = TRUE),
      mean_jackknife_refits = if (any(is.finite(jk_refits))) mean(jk_refits, na.rm = TRUE) else NA_real_,
      mean_jackknife_groups = if (any(is.finite(jk_groups))) mean(jk_groups, na.rm = TRUE) else NA_real_,
      mean_observed_fraction = mean(g$observed_fraction, na.rm = TRUE),
      mean_subjects_used = mean(g$subjects_used, na.rm = TRUE),
      mean_empty_rows = mean(g$empty_rows, na.rm = TRUE),
      mean_complete_rows = mean(g$complete_rows, na.rm = TRUE),
      min_pair_count_min = min(g$min_pair_count, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  ans <- do.call(rbind, pieces)
  ans <- ans[order(ans$dgp, ans$mechanism, ans$n, ans$weight_label,
                   ans$estimator, ans$coefficient), ]
  rownames(ans) <- NULL
  ans
}

log_progress <- function(...) {
  if (opts$progress) message(format(Sys.time(), "%H:%M:%S"), " ", sprintf(...))
}

write_csv_atomic <- function(x, path) {
  tmp <- paste0(path, ".tmp")
  write.csv(x, tmp, row.names = FALSE)
  if (file.exists(path)) unlink(path)
  ok <- file.rename(tmp, path)
  if (!ok) stop("Failed to move temporary file into place: ", path, call. = FALSE)
}

checkpoint_file <- function(cell, grid, dir) {
  file.path(
    dir,
    sprintf(
      "cell-%03d-%s-%s-n%d.csv",
      cell,
      grid$dgp[[cell]],
      grid$mechanism[[cell]],
      grid$n[[cell]]
    )
  )
}

truth_rows <- do.call(rbind, lapply(seq_along(opts$dgps), function(i) {
  spec <- dgp_defs[[opts$dgps[[i]]]]
  truth_for_dgp(spec, opts$seed_base + 900000L + i, opts$truth_n)
}))

grid <- expand.grid(
  dgp = opts$dgps,
  mechanism = opts$mechanisms,
  n = opts$n_grid,
  stringsAsFactors = FALSE
)
grid <- grid[order(grid$dgp, grid$mechanism, grid$n), , drop = FALSE]
grid$cell <- seq_len(nrow(grid))

mechanisms_out <- data.frame(
  mechanism = valid_mechanisms,
  family = c("sanity", "mcar", "mcar", "mcar", "mcar",
             "mar", "mar", "mar", "mar", "mar"),
  description = c(
    "No deletion; complete-data sanity row.",
    "Independent per-cell MCAR with 15 percent deletion.",
    "Independent per-cell MCAR with 30 percent deletion.",
    "Rater-specific MCAR deletion rates with mean near 30 percent.",
    "Subject-clustered MCAR deletion propensities with mean near 30 percent.",
    "Clean anchor MAR: rater 1 always observed; later deletion depends linearly on rater 1 category, target 15 percent.",
    "Clean anchor MAR: rater 1 always observed; later deletion depends linearly on rater 1 category, target 30 percent.",
    "Clean nonlinear anchor MAR: later deletion is highest at extreme rater 1 categories, target 30 percent.",
    "Clean shifted anchor MAR: anchor dependence plus rater-specific deletion shifts, target 30 percent on average.",
    "Sequential MAR: later deletion depends on the previous observed rating or its missingness, target 30 percent."
  ),
  stringsAsFactors = FALSE
)

dgp_rows <- do.call(rbind, lapply(dgp_defs[opts$dgps], function(spec) {
  data.frame(
    dgp = spec$id,
    dgp_label = spec$label,
    C = spec$C,
    R = spec$R,
    category_prob = paste(signif(spec$p, 4), collapse = ","),
    skill = paste(signif(spec$skill, 4), collapse = ","),
    truth_skill_mult = paste(signif(spec$truth_skill_mult, 4), collapse = ","),
    support_size = spec$C^spec$R,
    stringsAsFactors = FALSE
  )
}))

dir.create(opts$out_dir, recursive = TRUE, showWarnings = FALSE)
checkpoint_dir <- file.path(opts$out_dir, "checkpoints")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
write_csv_atomic(truth_rows, file.path(opts$out_dir, "truth.csv"))
write_csv_atomic(grid[, c("cell", "dgp", "mechanism", "n")],
                 file.path(opts$out_dir, "cell_plan.csv"))
write_csv_atomic(mechanisms_out, file.path(opts$out_dir, "mechanisms.csv"))
write_csv_atomic(dgp_rows, file.path(opts$out_dir, "dgps.csv"))

checkpoint_files <- character(nrow(grid))
t0 <- Sys.time()
for (cell in seq_len(nrow(grid))) {
  spec <- dgp_defs[[grid$dgp[[cell]]]]
  checkpoint_files[[cell]] <- checkpoint_file(cell, grid, checkpoint_dir)
  if (opts$resume && file.exists(checkpoint_files[[cell]])) {
    log_progress(
      "cell %d/%d: reusing %s",
      cell, nrow(grid), basename(checkpoint_files[[cell]])
    )
    next
  }

  cell_rows <- list()
  log_progress(
    "cell %d/%d: %s %s n=%d reps=%d",
    cell, nrow(grid), grid$dgp[[cell]], grid$mechanism[[cell]],
    grid$n[[cell]], opts$reps
  )
  for (rep_id in seq_len(opts$reps)) {
    seed <- opts$seed_base + 1000000L * cell + rep_id
    set.seed(seed)
    X_star <- simulate_ratings(grid$n[[cell]], spec)
    X <- apply_missing(X_star, spec, grid$mechanism[[cell]])
    observed_fraction <- mean(!is.na(X))
    empty_rows <- sum(rowSums(!is.na(X)) == 0L)
    keep <- rowSums(!is.na(X)) > 0L
    X <- X[keep, , drop = FALSE]
    subjects_used <- nrow(X)
    complete_rows <- sum(stats::complete.cases(X))
    min_pair <- min_pair_count(X)

    for (m in seq_len(nrow(method_defs))) {
      method_row <- method_defs[m, ]
      fit <- fit_one(X, spec, method_row)
      truth_key <- truth_rows$dgp == spec$id &
        truth_rows$weight_label == method_row$weight_label
      truth_now <- truth_rows[truth_key, ]
      fit$truth <- truth_now$truth[match(fit$coefficient, truth_now$coefficient)]
      fit$covered <- fit$lower <= fit$truth & fit$truth <= fit$upper
      fit$fisher_covered <- fit$fisher_lower <= fit$truth & fit$truth <= fit$fisher_upper

      cell_rows[[length(cell_rows) + 1L]] <- data.frame(
        dgp = spec$id,
        dgp_label = spec$label,
        C = spec$C,
        R = spec$R,
        mechanism = grid$mechanism[[cell]],
        n = grid$n[[cell]],
        rep = rep_id,
        seed = seed,
        method = method_row$method,
        estimator = method_row$estimator,
        weight = method_row$weight,
        weight_label = method_row$weight_label,
        fit,
        observed_fraction = observed_fraction,
        subjects_used = subjects_used,
        empty_rows = empty_rows,
        complete_rows = complete_rows,
        min_pair_count = min_pair,
        stringsAsFactors = FALSE
      )
    }
  }
  cell_replicates <- if (length(cell_rows)) do.call(rbind, cell_rows) else data.frame()
  write_csv_atomic(cell_replicates, checkpoint_files[[cell]])
  log_progress(
    "cell %d/%d: wrote %s (%d rows)",
    cell, nrow(grid), basename(checkpoint_files[[cell]]), nrow(cell_replicates)
  )
}

missing_checkpoints <- checkpoint_files[!file.exists(checkpoint_files)]
if (length(missing_checkpoints)) {
  stop("Missing checkpoint(s): ", paste(missing_checkpoints, collapse = ", "),
       call. = FALSE)
}
replicates <- do.call(rbind, lapply(checkpoint_files, function(path) {
  read.csv(path, stringsAsFactors = FALSE)
}))
summary <- summarize_replicates(replicates)

write_csv_atomic(replicates, file.path(opts$out_dir, "replicates.csv"))
write_csv_atomic(summary, file.path(opts$out_dir, "summary.csv"))

cmd <- paste(commandArgs(FALSE), collapse = " ")
metadata <- data.frame(
  key = c("generated_at", "script", "command", "mode", "reps", "n_grid",
          "truth_n", "dgps", "mechanisms", "methods", "seed_base",
          "resume", "checkpoint_dir", "elapsed_seconds",
          "misskappa_version", "r_version"),
  value = c(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
    normalizePath(file.path(script_dir, "run_experiment.R"), mustWork = FALSE),
    cmd,
    opts$mode,
    as.character(opts$reps),
    paste(opts$n_grid, collapse = ","),
    as.character(opts$truth_n),
    paste(opts$dgps, collapse = ","),
    paste(opts$mechanisms, collapse = ","),
    paste(opts$methods, collapse = ","),
    as.character(opts$seed_base),
    as.character(opts$resume),
    normalizePath(checkpoint_dir, mustWork = FALSE),
    as.character(round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 3)),
    as.character(utils::packageVersion("misskappa")),
    R.version.string
  ),
  stringsAsFactors = FALSE
)
write_csv_atomic(metadata, file.path(opts$out_dir, "metadata.csv"))

cat("Wrote:\n")
for (name in c("truth.csv", "cell_plan.csv", "replicates.csv", "summary.csv",
               "mechanisms.csv", "dgps.csv", "metadata.csv")) {
  cat("  ", normalizePath(file.path(opts$out_dir, name), mustWork = FALSE), "\n", sep = "")
}
cat("  ", normalizePath(checkpoint_dir, mustWork = FALSE), "/\n", sep = "")

headline <- summary[summary$coefficient == "Conger" &
                      summary$n == min(summary$n) &
                      summary$mechanism %in% c("mcar30", "mar_anchor30"),
                    c("dgp", "mechanism", "method", "weight_label",
                      "bias", "rmse", "coverage95", "n_valid")]
if (nrow(headline)) {
  cat("\nConger pilot rows:\n")
  print(format(headline, digits = 3), row.names = FALSE)
}
