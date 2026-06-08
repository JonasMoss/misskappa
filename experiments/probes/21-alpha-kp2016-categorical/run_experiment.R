#!/usr/bin/env Rscript
#
# Experiment 21: Kelley & Pornprasertmanit (2016)-style design for the
# alpha-missing paper, adapted to ORDERED-CATEGORICAL items.
#
# Kelley & Pornprasertmanit (2016, Psychological Methods 21:69-92) is the
# modern complete-data evaluation of confidence-interval methods for
# composite reliability. Their continuous Study 1 crosses sample size,
# number of items, factor-loading pattern (equal = tau-equivalent vs.
# unequal = congeneric), population reliability, and item distribution.
#
# This experiment keeps that design skeleton but generates ordered
# categorical items (Cat-FIML in misskappa needs them; pairwise and
# normal FIML do not, listwise is a sanity check). A single-factor
# standardized latent model is calibrated so the LATENT coefficient alpha
# hits a target in {.7,.8,.9}; the items are then categorised, and the
# reported population truth is the exact post-categorisation alpha (the
# Maydeu-Olivares, Coffman & Hartmann 2007 convention). Per the paper's
# preference the congeneric ("graded") loading pattern is the primary
# baseline; the equal/tau pattern is retained as a contrast.
#
# Cat-FIML is run only where the saturated C^J table is tractable
# (<= --cat-em-max-cells). For the longer K&P scales it is infeasible;
# those cells record a skip rather than a fit. That infeasibility is part
# of the estimator-selection story, not a bug.
#
# Outputs (stable column names):
#   results/truth.csv
#   results/replicates.csv
#   results/summary.csv
#   results/metadata.csv

suppressPackageStartupMessages({
  library(misskappa)
})

usage <- function(status = 0L) {
  cat(paste0(
    "Usage: Rscript run_experiment.R [options]\n\n",
    "Design factors (K&P 2016, categorical adaptation):\n",
    "  --items CSV         Items per scale J. Default: 4,8,20.\n",
    "  --alphas CSV        Target LATENT alpha. Default: 0.7,0.8,0.9.\n",
    "  --loadings CSV      Loading patterns: equal,graded. Default: both.\n",
    "  --categories CSV    Categories per item C. Default: 2,5.\n",
    "  --shapes CSV        Margin shapes: symmetric,skew. Default: both.\n",
    "  --n-grid CSV        Sample sizes. Default: 200,500.\n",
    "  --mechanisms CSV    Missingness. Default: complete,mcar30,mar_anchor30.\n",
    "  --methods CSV       Methods. Default: pairwise,cat_em,normal_fiml,listwise.\n\n",
    "Run control:\n",
    "  --reps N            Replicates per cell. Default: 40.\n",
    "  --seed-base N       Base seed. Default: 210000.\n",
    "  --cat-em-max-cells N  Skip Cat-FIML when C^J exceeds this. Default: 4096.\n",
    "  --out-dir PATH      Output directory. Default: script-local results/.\n",
    "  --progress          Print one line per design cell.\n",
    "  --smoke             Cheap check: tiny grid, reps=4.\n",
    "  --help, -h          Show this help and exit.\n\n",
    "Examples:\n",
    "  Rscript run_experiment.R --smoke --progress\n",
    "  Rscript run_experiment.R --items 4,8,12,16,20 --n-grid 50,100,200,400,1000 \\\n",
    "    --reps 500 --progress\n"
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
parse_csv_num <- function(x, arg) {
  out <- as.numeric(strsplit(x, ",", fixed = TRUE)[[1L]])
  if (any(is.na(out))) stop(arg, " must be a comma-separated numeric list.", call. = FALSE)
  unique(out)
}

parse_args <- function(argv) {
  opts <- list(
    items = c(4L, 8L, 20L),
    alphas = c(0.7, 0.8, 0.9),
    loadings = c("equal", "graded"),
    categories = c(2L, 5L),
    shapes = c("symmetric", "skew"),
    n_grid = c(200L, 500L),
    mechanisms = c("complete", "mcar30", "mar_anchor30"),
    methods = c("pairwise", "cat_em", "normal_fiml", "listwise"),
    reps = 40L,
    seed_base = 210000L,
    cat_em_max_cells = 4096,
    out_dir = file.path(script_dir, "results"),
    progress = FALSE,
    smoke = FALSE
  )
  explicit <- as.list(setNames(rep(FALSE, 8L),
    c("items", "alphas", "loadings", "categories", "shapes",
      "n_grid", "mechanisms", "methods")))

  i <- 1L
  while (i <= length(argv)) {
    arg <- argv[[i]]
    if (arg %in% c("--help", "-h")) usage(0L)
    if (arg == "--progress") { opts$progress <- TRUE; i <- i + 1L; next }
    if (arg == "--smoke") { opts$smoke <- TRUE; i <- i + 1L; next }
    needs_value <- c("--items", "--alphas", "--loadings", "--categories",
                     "--shapes", "--n-grid", "--mechanisms", "--methods",
                     "--reps", "--seed-base", "--cat-em-max-cells", "--out-dir")
    if (arg %in% needs_value) {
      if (i == length(argv)) stop(arg, " needs a value.", call. = FALSE)
      val <- argv[[i + 1L]]
      if (arg == "--items") { opts$items <- parse_csv_int(val, arg); explicit$items <- TRUE }
      else if (arg == "--alphas") { opts$alphas <- parse_csv_num(val, arg); explicit$alphas <- TRUE }
      else if (arg == "--loadings") { opts$loadings <- parse_csv_chr(val); explicit$loadings <- TRUE }
      else if (arg == "--categories") { opts$categories <- parse_csv_int(val, arg); explicit$categories <- TRUE }
      else if (arg == "--shapes") { opts$shapes <- parse_csv_chr(val); explicit$shapes <- TRUE }
      else if (arg == "--n-grid") { opts$n_grid <- parse_csv_int(val, arg); explicit$n_grid <- TRUE }
      else if (arg == "--mechanisms") { opts$mechanisms <- parse_csv_chr(val); explicit$mechanisms <- TRUE }
      else if (arg == "--methods") { opts$methods <- parse_csv_chr(val); explicit$methods <- TRUE }
      else if (arg == "--reps") opts$reps <- as.integer(val)
      else if (arg == "--seed-base") opts$seed_base <- as.integer(val)
      else if (arg == "--cat-em-max-cells") opts$cat_em_max_cells <- as.numeric(val)
      else if (arg == "--out-dir") opts$out_dir <- val
      i <- i + 2L
      next
    }
    stop("Unknown argument: ", arg, call. = FALSE)
  }

  if (opts$smoke) {
    if (!explicit$items) opts$items <- 4L
    if (!explicit$alphas) opts$alphas <- 0.8
    if (!explicit$loadings) opts$loadings <- "graded"
    if (!explicit$categories) opts$categories <- 5L
    if (!explicit$shapes) opts$shapes <- "symmetric"
    if (!explicit$n_grid) opts$n_grid <- 100L
    if (!explicit$mechanisms) opts$mechanisms <- c("complete", "mcar30")
    if (!explicit$methods) opts$methods <- c("pairwise", "cat_em", "normal_fiml", "listwise")
    opts$reps <- if (any(grepl("^--reps$", argv))) opts$reps else 4L
  }

  valid_loadings <- c("equal", "graded")
  valid_shapes <- c("symmetric", "skew")
  valid_mechanisms <- c("complete", "mcar15", "mcar30", "mar_anchor15", "mar_anchor30")
  valid_methods <- c("pairwise", "cat_em", "normal_fiml", "listwise")
  if (!all(opts$loadings %in% valid_loadings))
    stop("--loadings must be in: ", paste(valid_loadings, collapse = ","), call. = FALSE)
  if (!all(opts$shapes %in% valid_shapes))
    stop("--shapes must be in: ", paste(valid_shapes, collapse = ","), call. = FALSE)
  if (!all(opts$mechanisms %in% valid_mechanisms))
    stop("--mechanisms must be in: ", paste(valid_mechanisms, collapse = ","), call. = FALSE)
  if (!all(opts$methods %in% valid_methods))
    stop("--methods must be in: ", paste(valid_methods, collapse = ","), call. = FALSE)
  if (any(opts$items < 3L)) stop("--items entries must be >= 3.", call. = FALSE)
  if (any(opts$categories < 2L)) stop("--categories entries must be >= 2.", call. = FALSE)
  if (is.na(opts$seed_base)) stop("--seed-base must be an integer.", call. = FALSE)
  opts
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))

if (!exists("alpha", asNamespace("misskappa"), mode = "function")) {
  stop("misskappa::alpha() unavailable; install misskappa first.", call. = FALSE)
}
if (!"type" %in% names(formals(misskappa::alpha))) {
  stop("misskappa::alpha() has no 'type' argument; update to a build that ",
       "dispatches normal vs categorical FIML via type=.", call. = FALSE)
}

# ---------------------------------------------------------------------
# Margins: scores, category probabilities, latent thresholds per
# (categories, shape). "symmetric" margins are near-normal (D1-like);
# "skew" margins push mass to one end to enter the high-skew regime
# (MOCH 2007 / K&P D2-D4) where normal-theory inference is stressed.
# ---------------------------------------------------------------------

margin_spec <- function(C, shape) {
  scores <- 0:(C - 1L)
  if (shape == "symmetric") {
    if (C == 2L) {
      probs <- c(0.5, 0.5)
    } else {
      # discretised symmetric bell over C categories
      mids <- (seq_len(C) - 0.5) / C
      dens <- stats::dnorm(stats::qnorm(mids))
      probs <- dens / sum(dens)
    }
  } else if (shape == "skew") {
    if (C == 2L) {
      probs <- c(0.85, 0.15)          # |skewness| ~ 2
    } else {
      # right-skewed, rare upper categories
      probs <- (C:1)^2
      probs <- probs / sum(probs)
    }
  } else {
    stop("Unknown shape: ", shape, call. = FALSE)
  }
  cdf <- cumsum(probs)[-C]
  thresholds <- stats::qnorm(pmin(pmax(cdf, 1e-6), 1 - 1e-6))
  score_mean <- sum(probs * scores)
  score_var <- sum(probs * (scores - score_mean)^2)
  sk <- if (score_var > 0) sum(probs * (scores - score_mean)^3) / score_var^1.5 else 0
  list(C = C, shape = shape, scores = scores, probs = probs,
       thresholds = thresholds, score_mean = score_mean,
       score_var = score_var, skewness = sk)
}

# ---------------------------------------------------------------------
# Loadings + calibration to a target LATENT alpha.
# Single-factor standardised model: lambda_i loadings, latent item
# variance 1, latent correlation matrix R[i,j] = lambda_i lambda_j.
# ---------------------------------------------------------------------

base_loadings <- function(J, pattern) {
  if (pattern == "equal") rep(1, J)
  else if (pattern == "graded") seq(0.2, 0.8, length.out = J)
  else stop("Unknown loading pattern: ", pattern, call. = FALSE)
}

latent_alpha <- function(lambda) {
  J <- length(lambda)
  s1 <- sum(lambda); s2 <- sum(lambda^2)
  denom <- J + s1^2 - s2          # = sum of latent covariance (diag 1)
  (J / (J - 1)) * (1 - J / denom)
}

CAP_LAMBDA <- 0.95
calibrate_loadings <- function(J, pattern, target_alpha) {
  L0 <- base_loadings(J, pattern)
  c_max <- CAP_LAMBDA / max(L0)
  a_max <- latent_alpha(c_max * L0)
  clamped <- FALSE
  if (a_max <= target_alpha) {
    cc <- c_max
    clamped <- TRUE
  } else {
    cc <- stats::uniroot(
      function(c) latent_alpha(c * L0) - target_alpha,
      interval = c(1e-6, c_max), tol = 1e-10
    )$root
  }
  lambda <- cc * L0
  list(lambda = lambda, latent_alpha = latent_alpha(lambda),
       clamped = clamped, scale = cc)
}

# ---------------------------------------------------------------------
# Exact post-categorisation covariance from the latent Gaussian copula.
# ---------------------------------------------------------------------

bvn_rect_prob <- function(xlo, xhi, ylo, yhi, rho) {
  if (abs(rho) < 1e-12) {
    return((stats::pnorm(xhi) - stats::pnorm(xlo)) *
             (stats::pnorm(yhi) - stats::pnorm(ylo)))
  }
  s <- sqrt(1 - rho * rho)
  integrand <- function(x) {
    stats::dnorm(x) *
      (stats::pnorm((yhi - rho * x) / s) - stats::pnorm((ylo - rho * x) / s))
  }
  stats::integrate(integrand, lower = xlo, upper = xhi,
                   rel.tol = 1e-10, abs.tol = 1e-12, subdivisions = 200L)$value
}

ordinal_cov_cache <- new.env(parent = emptyenv())
ordinal_cov_from_latent <- function(rho, marg) {
  key <- sprintf("%s|%d|%.12f", marg$shape, marg$C, rho)
  if (exists(key, ordinal_cov_cache, inherits = FALSE))
    return(get(key, ordinal_cov_cache, inherits = FALSE))
  lower <- c(-Inf, marg$thresholds)
  upper <- c(marg$thresholds, Inf)
  S <- length(marg$scores)
  joint <- matrix(0, S, S)
  for (a in seq_len(S)) for (b in seq_len(S))
    joint[a, b] <- bvn_rect_prob(lower[a], upper[a], lower[b], upper[b], rho)
  cov <- sum(outer(marg$scores, marg$scores) * joint) - marg$score_mean^2
  assign(key, cov, ordinal_cov_cache)
  cov
}

sigma_for_dgp <- function(lambda, marg) {
  J <- length(lambda)
  R <- tcrossprod(lambda); diag(R) <- 1
  eig <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  if (min(eig) <= 1e-8) stop("Latent correlation matrix not positive definite.")
  Sigma <- matrix(0, J, J)
  diag(Sigma) <- marg$score_var
  for (j in seq_len(J - 1L)) for (k in (j + 1L):J) {
    Sigma[j, k] <- ordinal_cov_from_latent(R[j, k], marg)
    Sigma[k, j] <- Sigma[j, k]
  }
  Sigma
}

alpha_from_sigma <- function(Sigma) {
  p <- ncol(Sigma)
  (p / (p - 1)) * (1 - sum(diag(Sigma)) / sum(Sigma))
}

simulate_complete <- function(n, lambda, marg) {
  J <- length(lambda)
  R <- tcrossprod(lambda); diag(R) <- 1
  z <- matrix(stats::rnorm(n * J), n, J) %*% chol(R)
  x <- apply(z, 2L, function(col) findInterval(col, marg$thresholds))
  x <- matrix(as.integer(x), nrow = n, ncol = J)
  colnames(x) <- paste0("y", seq_len(J))
  x
}

# ---------------------------------------------------------------------
# Missingness (item 1 is the fully observed anchor; matches the paper).
# ---------------------------------------------------------------------

target_items <- function(p) seq.int(2L, p)

apply_mcar <- function(X, rate) {
  out <- X
  for (j in target_items(ncol(X)))
    out[stats::runif(nrow(X)) < rate, j] <- NA_integer_
  out
}

tune_intercept <- function(z, slope, rate, offset = 0.0) {
  f <- function(a) mean(stats::plogis(a + offset + slope * z)) - rate
  stats::uniroot(f, c(-12, 12))$root
}

apply_anchor_mar <- function(X, rate) {
  out <- X
  z <- as.numeric(scale(X[, 1L]))
  z[!is.finite(z)] <- 0
  slope <- -1.35
  p <- ncol(X)
  for (j in target_items(p)) {
    shift <- (j - (p + 2) / 2) * 0.08
    a <- tune_intercept(z, slope, rate, shift)
    p_miss <- stats::plogis(a + shift + slope * z)
    out[stats::runif(nrow(X)) < p_miss, j] <- NA_integer_
  }
  out
}

apply_missing <- function(X, mechanism) {
  switch(mechanism,
    complete = X,
    mcar15 = apply_mcar(X, 0.15),
    mcar30 = apply_mcar(X, 0.30),
    mar_anchor15 = apply_anchor_mar(X, 0.15),
    mar_anchor30 = apply_anchor_mar(X, 0.30),
    stop("Unknown mechanism: ", mechanism, call. = FALSE)
  )
}

# ---------------------------------------------------------------------
# Estimators (identical wiring to experiment 18).
# ---------------------------------------------------------------------

se_from_fit <- function(fit) {
  v <- as.numeric(stats::vcov(fit)[1L, 1L])
  if (!is.finite(v) || v < 0) return(NA_real_)
  sqrt(v)
}

fit_pairwise <- function(X) {
  fit <- misskappa::alpha(X, method = "available")
  list(estimate = as.numeric(stats::coef(fit)[["alpha"]]), se = se_from_fit(fit), iter = NA_integer_)
}

fit_cat_em <- function(X) {
  fit <- misskappa::alpha(
    X, method = "fiml", type = "categorical",
    em_options = list(tol = 1e-7, max_iter = 20000L, prune_tol = 1e-10,
                      start_alpha = 0.1, info_rcond = 1e-4)
  )
  list(estimate = as.numeric(stats::coef(fit)[["alpha"]]), se = se_from_fit(fit), iter = NA_integer_)
}

fit_normal_fiml <- function(X) {
  X_num <- matrix(as.numeric(X), nrow = nrow(X), ncol = ncol(X))
  fit <- misskappa::alpha(
    X_num, method = "fiml", type = "normal", se_type = "sandwich",
    em_options = list(tol = 1e-8, max_iter = 10000L, fd_h = 1e-5)
  )
  it <- tryCatch(fit$moments$iterations, error = function(e) NA_integer_)
  if (is.null(it)) it <- NA_integer_
  list(estimate = as.numeric(stats::coef(fit)[["alpha"]]), se = se_from_fit(fit), iter = it)
}

fit_listwise <- function(X) {
  cc <- stats::complete.cases(X)
  if (sum(cc) < ncol(X)) stop("too few complete rows for listwise alpha.")
  fit <- misskappa::alpha(X[cc, , drop = FALSE], method = "available")
  list(estimate = as.numeric(stats::coef(fit)[["alpha"]]), se = se_from_fit(fit), iter = NA_integer_)
}

fit_one <- function(X, method, cat_em_ok) {
  if (method == "cat_em" && !cat_em_ok) {
    return(list(estimate = NA_real_, se = NA_real_, iter = NA_integer_,
                elapsed_ms = NA_real_, skipped = TRUE,
                error = "cat_em skipped: C^J exceeds --cat-em-max-cells"))
  }
  start <- proc.time()[["elapsed"]]
  out <- tryCatch({
    ans <- switch(method,
      pairwise = fit_pairwise(X),
      cat_em = fit_cat_em(X),
      normal_fiml = fit_normal_fiml(X),
      listwise = fit_listwise(X))
    ans$error <- ""; ans$skipped <- FALSE; ans
  }, error = function(e) {
    list(estimate = NA_real_, se = NA_real_, iter = NA_integer_,
         skipped = FALSE, error = conditionMessage(e))
  })
  out$elapsed_ms <- 1000 * (proc.time()[["elapsed"]] - start)
  out
}

min_pair_count <- function(X) {
  p <- ncol(X); counts <- integer(0L)
  for (j in seq_len(p - 1L)) for (k in (j + 1L):p)
    counts <- c(counts, sum(!is.na(X[, j]) & !is.na(X[, k])))
  min(counts)
}

log_progress <- function(...) {
  if (opts$progress) message(format(Sys.time(), "%H:%M:%S"), " ", sprintf(...))
}

# ---------------------------------------------------------------------
# Enumerate the DGP grid and compute population truth per DGP.
# ---------------------------------------------------------------------

dgp_grid <- expand.grid(
  J = opts$items, target_alpha = opts$alphas, loading = opts$loadings,
  C = opts$categories, shape = opts$shapes,
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
dgp_grid$dgp <- sprintf("J%d_a%02d_%s_C%d_%s",
  dgp_grid$J, round(100 * dgp_grid$target_alpha), dgp_grid$loading,
  dgp_grid$C, dgp_grid$shape)

dgp_defs <- list()
truth_rows <- vector("list", nrow(dgp_grid))
for (r in seq_len(nrow(dgp_grid))) {
  g <- dgp_grid[r, ]
  marg <- margin_spec(g$C, g$shape)
  cal <- calibrate_loadings(g$J, g$loading, g$target_alpha)
  Sigma <- sigma_for_dgp(cal$lambda, marg)
  cells <- g$C^g$J
  dgp_defs[[g$dgp]] <- list(J = g$J, lambda = cal$lambda, marg = marg,
                            cat_em_ok = is.finite(cells) && cells <= opts$cat_em_max_cells)
  truth_rows[[r]] <- data.frame(
    dgp = g$dgp, n_items = g$J, loading = g$loading,
    categories = g$C, shape = g$shape, item_skewness = round(marg$skewness, 4),
    target_latent_alpha = g$target_alpha, latent_alpha = round(cal$latent_alpha, 6),
    alpha = alpha_from_sigma(Sigma), loading_clamped = cal$clamped,
    min_loading = round(min(cal$lambda), 4), max_loading = round(max(cal$lambda), 4),
    cells = cells, cat_em_feasible = dgp_defs[[g$dgp]]$cat_em_ok,
    truth_method = "exact post-categorisation alpha from Gaussian copula",
    stringsAsFactors = FALSE)
}
truth <- do.call(rbind, truth_rows)
truth_for <- function(dgp) truth$alpha[match(dgp, truth$dgp)]

# ---------------------------------------------------------------------
# Run grid: DGP x mechanism x n.
# ---------------------------------------------------------------------

grid <- expand.grid(dgp = dgp_grid$dgp, mechanism = opts$mechanisms,
                    n = opts$n_grid, stringsAsFactors = FALSE)
grid <- grid[order(grid$dgp, grid$mechanism, grid$n), , drop = FALSE]

row_seed <- function(cell, rep_id) opts$seed_base + 1000000L * cell + rep_id

rep_rows <- list()
t0 <- Sys.time()
for (cell in seq_len(nrow(grid))) {
  dgp_name <- grid$dgp[cell]
  spec <- dgp_defs[[dgp_name]]
  tr <- truth_for(dgp_name)
  log_progress("cell %d/%d: %s %s n=%d reps=%d (cat_em=%s)",
    cell, nrow(grid), dgp_name, grid$mechanism[cell], grid$n[cell],
    opts$reps, if (spec$cat_em_ok) "on" else "skip")
  for (rep_id in seq_len(opts$reps)) {
    set.seed(row_seed(cell, rep_id))
    X_star <- simulate_complete(grid$n[cell], spec$lambda, spec$marg)
    X <- apply_missing(X_star, grid$mechanism[cell])
    observed_fraction <- mean(!is.na(X))
    complete_rows <- sum(stats::complete.cases(X))
    min_pair <- min_pair_count(X)
    for (method in opts$methods) {
      fit <- fit_one(X, method, spec$cat_em_ok)
      rep_rows[[length(rep_rows) + 1L]] <- data.frame(
        dgp = dgp_name, n_items = spec$J, mechanism = grid$mechanism[cell],
        n = grid$n[cell], rep = rep_id, seed = row_seed(cell, rep_id),
        method = method, truth = tr, estimate = fit$estimate,
        bias = fit$estimate - tr, se = fit$se, elapsed_ms = fit$elapsed_ms,
        em_iter = fit$iter, skipped = isTRUE(fit$skipped),
        observed_fraction = observed_fraction, complete_rows = complete_rows,
        min_pair_count = min_pair, error = fit$error, stringsAsFactors = FALSE)
    }
  }
}
replicates <- if (length(rep_rows)) do.call(rbind, rep_rows) else data.frame()

# ---------------------------------------------------------------------
# Summaries.
# ---------------------------------------------------------------------

split_keys <- function(data, keys) interaction(data[, keys], drop = TRUE, lex.order = TRUE)

summarize_replicates <- function(df) {
  if (!nrow(df)) return(df)
  keys <- c("dgp", "n_items", "mechanism", "n", "method")
  out <- lapply(split(df, split_keys(df, keys)), function(g) {
    ok <- is.finite(g$estimate)
    se_ok <- ok & is.finite(g$se) & g$se > 0
    cov95 <- abs(g$estimate[se_ok] - g$truth[se_ok]) <= stats::qnorm(0.975) * g$se[se_ok]
    data.frame(
      dgp = g$dgp[[1L]], n_items = g$n_items[[1L]], mechanism = g$mechanism[[1L]],
      n = g$n[[1L]], method = g$method[[1L]], reps = length(unique(g$rep)),
      n_valid = sum(ok), failures = sum(!ok & !g$skipped), skipped = sum(g$skipped),
      truth = g$truth[[1L]],
      mean_estimate = if (any(ok)) mean(g$estimate[ok]) else NA_real_,
      bias = if (any(ok)) mean(g$estimate[ok] - g$truth[ok]) else NA_real_,
      sd_estimate = if (sum(ok) > 1L) stats::sd(g$estimate[ok]) else NA_real_,
      rmse = if (any(ok)) sqrt(mean((g$estimate[ok] - g$truth[ok])^2)) else NA_real_,
      mean_se = if (any(se_ok)) mean(g$se[se_ok]) else NA_real_,
      se_over_sd = if (sum(se_ok) > 1L && stats::sd(g$estimate[se_ok]) > 1e-12)
        mean(g$se[se_ok]) / stats::sd(g$estimate[se_ok]) else NA_real_,
      coverage95 = if (any(se_ok)) mean(cov95) else NA_real_,
      mean_ci_length = if (any(se_ok)) mean(2 * stats::qnorm(0.975) * g$se[se_ok]) else NA_real_,
      median_elapsed_ms = stats::median(g$elapsed_ms, na.rm = TRUE),
      mean_observed_fraction = mean(g$observed_fraction, na.rm = TRUE),
      stringsAsFactors = FALSE)
  })
  ans <- do.call(rbind, out)
  ans <- merge(ans, truth[, c("dgp", "loading", "categories", "shape",
                              "target_latent_alpha", "item_skewness")],
               by = "dgp", all.x = TRUE)
  ans <- ans[order(ans$dgp, ans$mechanism, ans$n, ans$method), ]
  rownames(ans) <- NULL
  ans
}
summary <- summarize_replicates(replicates)

# ---------------------------------------------------------------------
# Write outputs.
# ---------------------------------------------------------------------

dir.create(opts$out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(truth, file.path(opts$out_dir, "truth.csv"), row.names = FALSE)
write.csv(replicates, file.path(opts$out_dir, "replicates.csv"), row.names = FALSE)
write.csv(summary, file.path(opts$out_dir, "summary.csv"), row.names = FALSE)

metadata <- data.frame(
  key = c("generated_at", "script", "smoke", "reps", "items", "alphas",
          "loadings", "categories", "shapes", "n_grid", "mechanisms",
          "methods", "seed_base", "cat_em_max_cells", "n_dgps",
          "elapsed_seconds", "design", "misskappa_version", "r_version"),
  value = c(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
    normalizePath(file.path(script_dir, "run_experiment.R"), mustWork = FALSE),
    as.character(opts$smoke), as.character(opts$reps),
    paste(opts$items, collapse = ","), paste(opts$alphas, collapse = ","),
    paste(opts$loadings, collapse = ","), paste(opts$categories, collapse = ","),
    paste(opts$shapes, collapse = ","), paste(opts$n_grid, collapse = ","),
    paste(opts$mechanisms, collapse = ","), paste(opts$methods, collapse = ","),
    as.character(opts$seed_base), as.character(opts$cat_em_max_cells),
    as.character(nrow(dgp_grid)),
    as.character(round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 3)),
    "Kelley & Pornprasertmanit (2016) skeleton, ordered-categorical items, latent-alpha calibration, exact post-categorisation truth",
    as.character(utils::packageVersion("misskappa")), R.version.string),
  stringsAsFactors = FALSE)
write.csv(metadata, file.path(opts$out_dir, "metadata.csv"), row.names = FALSE)

message("Wrote ", nrow(dgp_grid), " DGPs x ", length(opts$mechanisms),
        " mechanisms x ", length(opts$n_grid), " sample sizes to ", opts$out_dir)
for (f in c("truth.csv", "replicates.csv", "summary.csv", "metadata.csv"))
  message("  ", file.path(opts$out_dir, f))
