#!/usr/bin/env Rscript
#
# Experiment 18: paper-facing continuous alpha simulation.
#
# First pass for the alpha-missing paper simulation. It compares pairwise
# covariance alpha, normal-FIML alpha with sandwich and normal-theory SEs, and a
# listwise strawman under paper-backed measurement models and clean MCAR/MAR
# item missingness. The pairwise point estimate is included now; pairwise SE /
# coverage stays NA until the package exposes the overlapping-subsample SE.
#
# Outputs:
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
    "Options:\n",
    "  --help, -h          Show this help and exit.\n",
    "  --smoke             Cheap check: n=160, reps=4, core DGP/mechanism grid.\n",
    "  --reps N            Replicates per cell. Default: 40.\n",
    "  --n-grid CSV        Sample sizes. Default: 200,500.\n",
    "  --dgps CSV          DGPs. Default: all continuous DGPs.\n",
    "  --dists CSV         Distributions: normal,nonnormal. Default: normal.\n",
    "  --mechanisms CSV    Missingness mechanisms. Default: complete,mcar15,mcar30,mar_zy_light,mar_anchor15,mar_anchor30.\n",
    "  --methods CSV       Methods. Default: pairwise,fiml_sandwich,fiml_normal,listwise.\n",
    "  --seed-base N       Base seed. Default: 181800.\n",
    "  --out-dir PATH      Output directory. Default: script-local results/.\n",
    "  --progress          Print one line per design cell.\n\n",
    "Examples:\n",
    "  Rscript run_experiment.R --smoke\n",
    "  Rscript run_experiment.R --reps 500 --n-grid 200,500,1000 --progress\n"
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

parse_args <- function(argv) {
  opts <- list(
    smoke = FALSE,
    reps = 40L,
    n_grid = c(200L, 500L),
    dgps = c("tau6_essential", "parallel6_zy", "congeneric6_zy",
             "congeneric8_gradient", "twofactor8"),
    dists = "normal",
    mechanisms = c("complete", "mcar15", "mcar30", "mar_zy_light",
                   "mar_anchor15", "mar_anchor30"),
    methods = c("pairwise", "fiml_sandwich", "fiml_normal", "listwise"),
    seed_base = 181800L,
    out_dir = file.path(script_dir, "results"),
    progress = FALSE
  )
  explicit <- list(reps = FALSE, n_grid = FALSE, dgps = FALSE,
                   dists = FALSE, mechanisms = FALSE, methods = FALSE)

  i <- 1L
  while (i <= length(argv)) {
    arg <- argv[[i]]
    if (arg %in% c("--help", "-h")) usage(0L)
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
    needs_value <- c("--reps", "--n-grid", "--dgps", "--dists",
                     "--mechanisms", "--methods", "--seed-base", "--out-dir")
    if (arg %in% needs_value) {
      if (i == length(argv)) stop(arg, " needs a value.", call. = FALSE)
      val <- argv[[i + 1L]]
      if (arg == "--reps") {
        opts$reps <- as.integer(val); explicit$reps <- TRUE
      } else if (arg == "--n-grid") {
        opts$n_grid <- parse_csv_int(val, arg); explicit$n_grid <- TRUE
      } else if (arg == "--dgps") {
        opts$dgps <- parse_csv_chr(val); explicit$dgps <- TRUE
      } else if (arg == "--dists") {
        opts$dists <- parse_csv_chr(val); explicit$dists <- TRUE
      } else if (arg == "--mechanisms") {
        opts$mechanisms <- parse_csv_chr(val); explicit$mechanisms <- TRUE
      } else if (arg == "--methods") {
        opts$methods <- parse_csv_chr(val); explicit$methods <- TRUE
      } else if (arg == "--seed-base") {
        opts$seed_base <- as.integer(val)
      } else if (arg == "--out-dir") {
        opts$out_dir <- val
      }
      i <- i + 2L
      next
    }
    stop("Unknown argument: ", arg, call. = FALSE)
  }

  if (opts$smoke) {
    if (!explicit$reps) opts$reps <- 4L
    if (!explicit$n_grid) opts$n_grid <- 160L
    if (!explicit$dgps) opts$dgps <- c("tau6_essential", "congeneric6_zy", "twofactor8")
    if (!explicit$dists) opts$dists <- "normal"
    if (!explicit$mechanisms) opts$mechanisms <- c("complete", "mcar30", "mar_anchor30")
    if (!explicit$methods) opts$methods <- c("pairwise", "fiml_sandwich",
                                             "fiml_normal", "listwise")
  }

  if (is.na(opts$reps) || opts$reps < 1L) stop("--reps must be >= 1.", call. = FALSE)
  if (any(is.na(opts$n_grid)) || any(opts$n_grid < 20L)) {
    stop("--n-grid entries must be integers >= 20.", call. = FALSE)
  }
  if (is.na(opts$seed_base)) stop("--seed-base must be an integer.", call. = FALSE)
  opts
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))

dgp_defs <- list(
  tau6_essential = list(
    label = "Essential-tau, 6 items",
    loadings = rep(sqrt(0.60), 6L),
    residual_var = c(0.30, 0.34, 0.38, 0.42, 0.46, 0.50),
    intercepts = c(-0.35, -0.20, -0.05, 0.05, 0.20, 0.35),
    family = "essential_tau"
  ),
  parallel6_zy = list(
    label = "Parallel Zhang-Yuan, 6 items",
    loadings = rep(sqrt(0.60), 6L),
    residual_var = rep(0.40, 6L),
    intercepts = rep(0.0, 6L),
    family = "parallel"
  ),
  congeneric6_zy = list(
    label = "Congeneric Zhang-Yuan, 6 items",
    loadings = sqrt(c(0.20, 0.20, 0.20, 0.60, 0.60, 0.60)),
    residual_var = c(0.80, 0.80, 0.80, 0.40, 0.40, 0.40),
    intercepts = rep(0.0, 6L),
    family = "congeneric"
  ),
  congeneric8_gradient = list(
    label = "Congeneric gradient, 8 items",
    loadings = seq(0.45, 0.85, length.out = 8L),
    residual_var = 1 - seq(0.45, 0.85, length.out = 8L)^2,
    intercepts = seq(-0.35, 0.35, length.out = 8L),
    family = "congeneric"
  ),
  twofactor8 = list(
    label = "Two-factor, 8 items",
    loadings = rep(0.75, 8L),
    residual_var = rep(1 - 0.75^2, 8L),
    intercepts = seq(-0.20, 0.20, length.out = 8L),
    family = "twofactor",
    factor_cor = 0.35
  )
)

valid_dgps <- names(dgp_defs)
valid_dists <- c("normal", "nonnormal")
valid_mechanisms <- c("complete", "mcar15", "mcar30", "mar_zy_light",
                      "mar_anchor15", "mar_anchor30", "mar_anchor_nonlinear30")
valid_methods <- c("pairwise", "fiml_sandwich", "fiml_normal", "listwise")

if (!all(opts$dgps %in% valid_dgps)) {
  stop("--dgps must contain only: ", paste(valid_dgps, collapse = ","), call. = FALSE)
}
if (!all(opts$dists %in% valid_dists)) {
  stop("--dists must contain only: ", paste(valid_dists, collapse = ","), call. = FALSE)
}
if (!all(opts$mechanisms %in% valid_mechanisms)) {
  stop("--mechanisms must contain only: ", paste(valid_mechanisms, collapse = ","), call. = FALSE)
}
if (!all(opts$methods %in% valid_methods)) {
  stop("--methods must contain only: ", paste(valid_methods, collapse = ","), call. = FALSE)
}
if (!exists("alpha_continuous", asNamespace("misskappa"), mode = "function")) {
  stop("misskappa::alpha_continuous() is not available; install the normal-FIML alpha package build first.",
       call. = FALSE)
}

alpha_from_sigma <- function(Sigma) {
  p <- ncol(Sigma)
  (p / (p - 1)) * (1 - sum(diag(Sigma)) / sum(Sigma))
}

sigma_for_dgp <- function(spec) {
  p <- length(spec$loadings)
  if (identical(spec$family, "twofactor")) {
    Lambda <- matrix(0, p, 2L)
    Lambda[1:(p / 2L), 1L] <- spec$loadings[1:(p / 2L)]
    Lambda[(p / 2L + 1L):p, 2L] <- spec$loadings[(p / 2L + 1L):p]
    Phi <- matrix(c(1, spec$factor_cor, spec$factor_cor, 1), 2L, 2L)
    Lambda %*% Phi %*% t(Lambda) + diag(spec$residual_var)
  } else {
    tcrossprod(spec$loadings) + diag(spec$residual_var)
  }
}

gen_factor_scores <- function(n, dist, q) {
  if (dist == "normal") {
    matrix(stats::rnorm(n * q), n, q)
  } else {
    z <- matrix((stats::rchisq(n * q, df = 3) - 3) / sqrt(6), n, q)
    z
  }
}

simulate_complete <- function(n, spec, dist) {
  p <- length(spec$loadings)
  if (identical(spec$family, "twofactor")) {
    F <- gen_factor_scores(n, dist, 2L)
    Phi_chol <- chol(matrix(c(1, spec$factor_cor, spec$factor_cor, 1), 2L, 2L))
    F <- F %*% Phi_chol
    X <- matrix(0, n, p)
    half <- p / 2L
    X[, 1:half] <- F[, 1L] %o% spec$loadings[1:half]
    X[, (half + 1L):p] <- F[, 2L] %o% spec$loadings[(half + 1L):p]
  } else {
    F <- gen_factor_scores(n, dist, 1L)
    X <- F[, 1L] %o% spec$loadings
  }
  errors <- matrix(stats::rnorm(n * p), n, p)
  X <- X + errors * rep(sqrt(spec$residual_var), each = n)
  X <- sweep(X, 2L, spec$intercepts, "+")
  colnames(X) <- paste0("y", seq_len(p))
  X
}

target_items <- function(p) seq.int(floor(p / 2L) + 1L, p)

apply_mcar <- function(X, rate) {
  out <- X
  targets <- target_items(ncol(X))
  miss <- matrix(stats::runif(nrow(X) * length(targets)) < rate,
                 nrow(X), length(targets))
  out[, targets][miss] <- NA_real_
  out
}

tune_intercept <- function(z, slope, rate) {
  f <- function(a) mean(stats::plogis(a + slope * z)) - rate
  stats::uniroot(f, c(-12, 12))$root
}

apply_anchor_mar <- function(X, rate, nonlinear = FALSE) {
  out <- X
  z <- as.numeric(scale(X[, 1L]))
  if (nonlinear) z <- abs(z) - mean(abs(z))
  slope <- if (nonlinear) 1.35 else -1.35
  a <- tune_intercept(z, slope, rate)
  base <- stats::plogis(a + slope * z)
  for (j in 2:ncol(X)) {
    shift <- (j - (ncol(X) + 1) / 2) * 0.08
    p_miss <- pmin(pmax(stats::plogis(stats::qlogis(base) + shift), 0.001), 0.999)
    out[stats::runif(nrow(X)) < p_miss, j] <- NA_real_
  }
  out
}

apply_zy_mar <- function(X) {
  if (ncol(X) != 6L) return(NULL)
  out <- X
  q1 <- stats::quantile(X[, 1L], probs = c(0.1, 0.2), type = 8)
  q4 <- stats::quantile(X[, 4L], probs = c(0.8, 0.9), type = 8)
  out[X[, 1L] <= q1[1L], 5L] <- NA_real_
  out[X[, 1L] > q1[1L] & X[, 1L] <= q1[2L], 6L] <- NA_real_
  out[X[, 4L] >= q4[2L], 2L] <- NA_real_
  out[X[, 4L] >= q4[1L] & X[, 4L] < q4[2L], 3L] <- NA_real_
  out
}

apply_missing <- function(X, mechanism) {
  if (mechanism == "complete") return(X)
  if (mechanism == "mcar15") return(apply_mcar(X, 0.15))
  if (mechanism == "mcar30") return(apply_mcar(X, 0.30))
  if (mechanism == "mar_zy_light") return(apply_zy_mar(X))
  if (mechanism == "mar_anchor15") return(apply_anchor_mar(X, 0.15))
  if (mechanism == "mar_anchor30") return(apply_anchor_mar(X, 0.30))
  if (mechanism == "mar_anchor_nonlinear30") return(apply_anchor_mar(X, 0.30, nonlinear = TRUE))
  stop("Unknown mechanism: ", mechanism, call. = FALSE)
}

fit_pairwise <- function(X) {
  S <- stats::cov(X, use = "pairwise.complete.obs")
  if (any(!is.finite(S))) stop("pairwise covariance is not finite.")
  list(estimate = alpha_from_sigma(S), se = NA_real_, iter = NA_integer_)
}

fit_fiml <- function(X, se_type) {
  fit <- misskappa::alpha_continuous(
    X,
    se_type = se_type,
    em_options = list(tol = 1e-8, max_iter = 10000L, fd_h = 1e-5)
  )
  list(
    estimate = as.numeric(stats::coef(fit)[["alpha"]]),
    se = sqrt(as.numeric(stats::vcov(fit)[1L, 1L])),
    iter = fit$moments$iterations
  )
}

fit_listwise <- function(X) {
  cc <- stats::complete.cases(X)
  if (sum(cc) <= ncol(X) + 2L) stop("too few complete rows for listwise alpha.")
  fit <- misskappa::alpha_continuous(
    X[cc, , drop = FALSE],
    se_type = "sandwich",
    em_options = list(tol = 1e-8, max_iter = 10000L, fd_h = 1e-5)
  )
  list(
    estimate = as.numeric(stats::coef(fit)[["alpha"]]),
    se = sqrt(as.numeric(stats::vcov(fit)[1L, 1L])),
    iter = fit$moments$iterations
  )
}

fit_one <- function(X, method) {
  start <- proc.time()[["elapsed"]]
  out <- tryCatch({
    ans <- switch(
      method,
      pairwise = fit_pairwise(X),
      fiml_sandwich = fit_fiml(X, "sandwich"),
      fiml_normal = fit_fiml(X, "normal"),
      listwise = fit_listwise(X)
    )
    ans$error <- ""
    ans
  }, error = function(e) {
    list(estimate = NA_real_, se = NA_real_, iter = NA_integer_,
         error = conditionMessage(e))
  })
  out$elapsed_ms <- 1000 * (proc.time()[["elapsed"]] - start)
  out
}

min_pair_count <- function(X) {
  p <- ncol(X)
  counts <- integer(p * (p - 1L) / 2L)
  k <- 0L
  for (j in seq_len(p - 1L)) {
    for (l in (j + 1L):p) {
      k <- k + 1L
      counts[k] <- sum(is.finite(X[, j]) & is.finite(X[, l]))
    }
  }
  min(counts)
}

row_seed <- function(cell, rep_id) opts$seed_base + 1000000L * cell + rep_id

truth_rows <- do.call(rbind, lapply(opts$dgps, function(name) {
  spec <- dgp_defs[[name]]
  S <- sigma_for_dgp(spec)
  data.frame(
    dgp = name,
    dgp_label = spec$label,
    family = spec$family,
    p = ncol(S),
    truth_alpha = alpha_from_sigma(S),
    min_loading = min(spec$loadings),
    max_loading = max(spec$loadings),
    min_residual_var = min(spec$residual_var),
    max_residual_var = max(spec$residual_var),
    truth_method = "analytic covariance",
    stringsAsFactors = FALSE
  )
}))

truth_for <- function(dgp) truth_rows$truth_alpha[match(dgp, truth_rows$dgp)]

split_keys <- function(data, keys) interaction(data[, keys], drop = TRUE, lex.order = TRUE)

summarize_replicates <- function(df) {
  keys <- c("dgp", "dgp_label", "family", "p", "dist", "mechanism", "n", "method")
  out <- lapply(split(df, split_keys(df, keys)), function(g) {
    ok <- is.finite(g$estimate)
    se_ok <- ok & is.finite(g$se) & g$se > 0
    z <- (g$estimate[se_ok] - g$truth[se_ok]) / g$se[se_ok]
    cov95 <- abs(g$estimate[se_ok] - g$truth[se_ok]) <= 1.96 * g$se[se_ok]
    data.frame(
      dgp = g$dgp[[1L]],
      dgp_label = g$dgp_label[[1L]],
      family = g$family[[1L]],
      p = g$p[[1L]],
      dist = g$dist[[1L]],
      mechanism = g$mechanism[[1L]],
      n = g$n[[1L]],
      method = g$method[[1L]],
      reps = length(unique(g$rep)),
      n_valid = sum(ok),
      failures = sum(!ok),
      truth = g$truth[[1L]],
      mean_estimate = if (any(ok)) mean(g$estimate[ok]) else NA_real_,
      bias = if (any(ok)) mean(g$estimate[ok] - g$truth[ok]) else NA_real_,
      sd_estimate = if (sum(ok) > 1L) stats::sd(g$estimate[ok]) else NA_real_,
      mc_se_bias = if (sum(ok) > 1L) stats::sd(g$estimate[ok]) / sqrt(sum(ok)) else NA_real_,
      rmse = if (any(ok)) sqrt(mean((g$estimate[ok] - g$truth[ok])^2)) else NA_real_,
      mean_se = if (any(se_ok)) mean(g$se[se_ok]) else NA_real_,
      se_over_sd = if (sum(se_ok) > 1L && stats::sd(g$estimate[se_ok]) > 1e-12) {
        mean(g$se[se_ok]) / stats::sd(g$estimate[se_ok])
      } else NA_real_,
      coverage95 = if (any(se_ok)) mean(cov95) else NA_real_,
      mean_z = if (length(z)) mean(z) else NA_real_,
      sd_z = if (length(z) > 1L) stats::sd(z) else NA_real_,
      mean_elapsed_ms = mean(g$elapsed_ms, na.rm = TRUE),
      median_elapsed_ms = stats::median(g$elapsed_ms, na.rm = TRUE),
      mean_observed_fraction = mean(g$observed_fraction, na.rm = TRUE),
      mean_complete_rows = mean(g$complete_rows, na.rm = TRUE),
      min_pair_count_min = min(g$min_pair_count, na.rm = TRUE),
      mean_em_iter = mean(g$em_iter, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  ans <- do.call(rbind, out)
  ans <- ans[order(ans$dgp, ans$dist, ans$mechanism, ans$n, ans$method), ]
  rownames(ans) <- NULL
  ans
}

log_progress <- function(...) {
  if (opts$progress) message(format(Sys.time(), "%H:%M:%S"), " ", sprintf(...))
}

grid <- expand.grid(
  dgp = opts$dgps,
  dist = opts$dists,
  mechanism = opts$mechanisms,
  n = opts$n_grid,
  stringsAsFactors = FALSE
)
grid <- grid[vapply(seq_len(nrow(grid)), function(i) {
  p <- length(dgp_defs[[grid$dgp[i]]]$loadings)
  !(grid$mechanism[i] == "mar_zy_light" && p != 6L)
}, logical(1)), , drop = FALSE]
grid <- grid[order(grid$dgp, grid$dist, grid$mechanism, grid$n), , drop = FALSE]

rep_rows <- list()
t0 <- Sys.time()
for (cell in seq_len(nrow(grid))) {
  dgp_name <- grid$dgp[cell]
  spec <- dgp_defs[[dgp_name]]
  truth <- truth_for(dgp_name)
  log_progress(
    "cell %d/%d: %s %s %s n=%d reps=%d",
    cell, nrow(grid), dgp_name, grid$dist[cell], grid$mechanism[cell],
    grid$n[cell], opts$reps
  )
  for (rep_id in seq_len(opts$reps)) {
    set.seed(row_seed(cell, rep_id))
    X_star <- simulate_complete(grid$n[cell], spec, grid$dist[cell])
    X <- apply_missing(X_star, grid$mechanism[cell])
    if (is.null(X)) next
    observed_fraction <- mean(is.finite(X))
    complete_rows <- sum(stats::complete.cases(X))
    min_pair <- min_pair_count(X)
    for (method in opts$methods) {
      fit <- fit_one(X, method)
      rep_rows[[length(rep_rows) + 1L]] <- data.frame(
        dgp = dgp_name,
        dgp_label = spec$label,
        family = spec$family,
        p = length(spec$loadings),
        dist = grid$dist[cell],
        mechanism = grid$mechanism[cell],
        n = grid$n[cell],
        rep = rep_id,
        seed = row_seed(cell, rep_id),
        method = method,
        truth = truth,
        estimate = fit$estimate,
        bias = fit$estimate - truth,
        se = fit$se,
        elapsed_ms = fit$elapsed_ms,
        em_iter = fit$iter,
        observed_fraction = observed_fraction,
        complete_rows = complete_rows,
        min_pair_count = min_pair,
        error = fit$error,
        stringsAsFactors = FALSE
      )
    }
  }
}

replicates <- if (length(rep_rows)) do.call(rbind, rep_rows) else data.frame()
summary <- summarize_replicates(replicates)

dir.create(opts$out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(truth_rows, file.path(opts$out_dir, "truth.csv"), row.names = FALSE)
write.csv(replicates, file.path(opts$out_dir, "replicates.csv"), row.names = FALSE)
write.csv(summary, file.path(opts$out_dir, "summary.csv"), row.names = FALSE)

metadata <- data.frame(
  key = c(
    "generated_at", "script", "smoke", "reps", "n_grid", "dgps", "dists",
    "mechanisms", "methods", "seed_base", "elapsed_seconds",
    "misskappa_version", "r_version", "pairwise_se_status"
  ),
  value = c(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
    normalizePath(file.path(script_dir, "run_experiment.R"), mustWork = FALSE),
    as.character(opts$smoke),
    as.character(opts$reps),
    paste(opts$n_grid, collapse = ","),
    paste(opts$dgps, collapse = ","),
    paste(opts$dists, collapse = ","),
    paste(opts$mechanisms, collapse = ","),
    paste(opts$methods, collapse = ","),
    as.character(opts$seed_base),
    as.character(round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 3)),
    as.character(utils::packageVersion("misskappa")),
    R.version.string,
    "pending package port of overlapping-subsample pairwise SE"
  ),
  stringsAsFactors = FALSE
)
write.csv(metadata, file.path(opts$out_dir, "metadata.csv"), row.names = FALSE)

cat("Wrote:\n")
cat("  ", normalizePath(file.path(opts$out_dir, "truth.csv"), mustWork = FALSE), "\n", sep = "")
cat("  ", normalizePath(file.path(opts$out_dir, "replicates.csv"), mustWork = FALSE), "\n", sep = "")
cat("  ", normalizePath(file.path(opts$out_dir, "summary.csv"), mustWork = FALSE), "\n", sep = "")
cat("  ", normalizePath(file.path(opts$out_dir, "metadata.csv"), mustWork = FALSE), "\n", sep = "")
