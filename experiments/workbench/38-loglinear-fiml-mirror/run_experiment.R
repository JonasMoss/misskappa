#!/usr/bin/env Rscript
#
# Experiment 38: is the all-2-way-interaction log-linear FIML the categorical
# mirror of nt-fiml?
#
# nt-fiml is consistent under MCAR for any distribution and efficient under MAR
# when normality holds, because the Gaussian is saturated in exactly the moments
# quadratic kappa needs (means + covariances) and max-entropy beyond. The
# discrete analogue for nominal-weight Conger/Fleiss -- which depend only on the
# univariate + bivariate margins -- is the all-2-way log-linear model
# (max-entropy given 1-way and 2-way margins). This probe tests, by simulation,
# whether that mirror holds: model-free consistency under MCAR (even when the
# truth has genuine 3-way+ rater interactions), and consistency + efficiency
# under MAR when the restriction holds. Point estimates, delta-method SEs, and
# exact truth.
#
# Hypotheses:
#   H1  MCAR model-free consistency: on a 3-way-interaction truth, loglin2 bias
#       -> 0 as n grows, unlike ds_fiml (nonzero misspecification asymptote).
#   H2  MAR efficiency under the restriction: on an all-2-way truth under MAR,
#       loglin2 is consistent and lower-RMSE than the MCAR-only baselines.
#   H3  MAR robustness boundary: on a 3-way-interaction truth under MAR,
#       loglin2 carries the fill-in misspecification bias; characterise it.
#   H4  SE calibration: loglin2's delta-method SE (Louis information in the
#       sufficient-statistic coordinates, added in 1187547) tracks the sampling
#       SD (se/sd -> 1) and gives nominal Wald coverage wherever the point
#       estimate is consistent; where it is biased (H3), coverage about the
#       truth drops while se/sd stays ~1 -- a bias effect, not an SE defect.

suppressPackageStartupMessages(library(misskappa))

usage <- function(status = 0L) {
  cat(paste0(
    "Usage: Rscript run_experiment.R [mode/options]\n\n",
    "Modes:\n",
    "  --smoke         Mechanical check: n=40, reps=5, one DGP, nominal.\n",
    "  --main          Full grid (default).\n\n",
    "Options:\n",
    "  --reps N        Replications per cell.\n",
    "  --n-grid CSV    Sample sizes.\n",
    "  --dgps CSV      DGP ids.\n",
    "  --mechanisms CSV  Missingness mechanism ids.\n",
    "  --weights CSV   Weight schemes (nominal, quadratic).\n",
    "  --seed-base N   Base seed. Default: 380000.\n",
    "  --cores N       Worker processes. Default: min(12, detected - 1).\n",
    "  --out-dir PATH  Output directory. Default: script-local results/.\n",
    "  --help, -h      Show this help and exit.\n"
  ))
  quit(save = "no", status = status)
}

script_dir <- local({
  cmd <- commandArgs(FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) == 0L) getwd()
  else dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE))
})

source(file.path(script_dir, "R", "loglin2_fiml.R"))
source(file.path(script_dir, "..", "..", "studies",
                 "29-kappa-planned-missingness", "R", "ds_fiml.R"))

parse_csv_chr <- function(x) trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
parse_csv_int <- function(x) as.integer(strsplit(x, ",", fixed = TRUE)[[1L]])

opts <- list(
  mode = "main", reps = NULL, n_grid = NULL, dgps = NULL, mechanisms = NULL,
  weights = NULL, seed_base = 380000L,
  cores = max(1L, min(12L, parallel::detectCores() - 1L)),
  out_dir = file.path(script_dir, "results")
)
args <- commandArgs(TRUE)
i <- 1L
while (i <= length(args)) {
  arg <- args[[i]]
  if (arg %in% c("--help", "-h")) usage()
  if (arg == "--smoke") { opts$mode <- "smoke"; i <- i + 1L; next }
  if (arg == "--main")  { opts$mode <- "main";  i <- i + 1L; next }
  if (!arg %in% c("--reps", "--n-grid", "--dgps", "--mechanisms", "--weights",
                  "--seed-base", "--cores", "--out-dir")) {
    stop("Unknown argument: ", arg, " (see --help)", call. = FALSE)
  }
  if (i == length(args)) stop(arg, " needs a value.", call. = FALSE)
  val <- args[[i + 1L]]
  if (arg == "--reps") opts$reps <- as.integer(val)
  if (arg == "--n-grid") opts$n_grid <- parse_csv_int(val)
  if (arg == "--dgps") opts$dgps <- parse_csv_chr(val)
  if (arg == "--mechanisms") opts$mechanisms <- parse_csv_chr(val)
  if (arg == "--weights") opts$weights <- parse_csv_chr(val)
  if (arg == "--seed-base") opts$seed_base <- as.integer(val)
  if (arg == "--cores") opts$cores <- as.integer(val)
  if (arg == "--out-dir") opts$out_dir <- val
  i <- i + 2L
}
smoke <- opts$mode == "smoke"
dir.create(opts$out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- DGPs --------------------------------------------------------------------
# Each spec carries C, R, a simulate(n) function, and exact true per-rater
# marginals `marg` (C x R) + pairwise joints `P[[j]][[l]]` (C x C) so the
# closed-form truth is the same fixed-rater estimand for every estimator.

ds_row_prob <- function(C, truth, accuracy, severity) {
  cats <- seq_len(C)
  center <- truth + severity
  error <- exp(-0.5 * ((cats - center) / 0.85)^2)
  error[truth] <- 0
  if (sum(error) <= 0) error <- rep(1, C)
  error <- error / sum(error)
  p <- (1 - accuracy) * error
  p[truth] <- p[truth] + accuracy
  p / sum(p)
}

make_ds_spec <- function(id, C, R, pi_true, accuracy, severity) {
  # True per-class emission matrices theta_j (C x C: row = true class).
  theta <- lapply(seq_len(R), function(j) {
    t(vapply(seq_len(C), function(k) ds_row_prob(C, k, accuracy[j], severity[j]),
             numeric(C)))
  })
  marg <- vapply(seq_len(R), function(j) as.numeric(colSums(pi_true * theta[[j]])),
                 numeric(C))
  P <- lapply(seq_len(R), function(j) {
    lapply(seq_len(R), function(l) {
      if (l <= j) return(NULL)
      t(theta[[j]]) %*% (pi_true * theta[[l]])     # C x C bivariate joint
    })
  })
  simulate <- function(n) {
    truth <- sample.int(C, n, replace = TRUE, prob = pi_true)
    X <- matrix(NA_integer_, n, R)
    for (j in seq_len(R)) for (k in seq_len(C)) {
      idx <- which(truth == k)
      if (length(idx)) X[idx, j] <- sample.int(C, length(idx), replace = TRUE,
                                               prob = theta[[j]][k, ])
    }
    X
  }
  list(id = id, C = C, R = R, kind = "threeway", marg = marg, P = P,
       simulate = simulate)
}

# All-2-way (no 3-way+) truth: log p = sum_j u_j(y_j) + beta * sum_{j<l} 1{y_j=y_l}.
make_loglin2_spec <- function(id, C, R, beta, main_sd = 0.25, seed = 11L) {
  set.seed(seed)
  u <- matrix(stats::rnorm(C * R, sd = main_sd), C, R)   # main effects, C x R
  grid <- as.matrix(expand.grid(rep(list(seq_len(C)), R)))  # C^R x R
  logp <- as.numeric(apply(grid, 1L, function(y) {
    s <- sum(u[cbind(y, seq_len(R))])
    for (j in seq_len(R - 1L)) for (l in (j + 1L):R) if (y[j] == y[l]) s <- s + beta
    s
  }))
  p <- exp(logp - max(logp)); p <- p / sum(p)
  tab <- array(p, dim = rep(C, R))
  marg <- vapply(seq_len(R), function(j) as.numeric(apply(tab, j, sum)), numeric(C))
  P <- lapply(seq_len(R), function(j) {
    lapply(seq_len(R), function(l) if (l <= j) NULL else apply(tab, c(j, l), sum))
  })
  cells <- seq_along(p)
  simulate <- function(n) {
    draws <- sample(cells, n, replace = TRUE, prob = p)
    idx <- arrayInd(draws, .dim = rep(C, R))
    matrix(as.integer(idx), n, R)
  }
  list(id = id, C = C, R = R, kind = "loglin2", marg = marg, P = P,
       table = tab, simulate = simulate)
}

dgp_defs <- list(
  loglin2_4x4 = make_loglin2_spec("loglin2_4x4", 4L, 4L, beta = 1.1),
  ds_4x4 = make_ds_spec("ds_4x4", 4L, 4L,
                        pi_true = c(0.40, 0.30, 0.20, 0.10),
                        accuracy = c(0.85, 0.80, 0.75, 0.70),
                        severity = c(-0.45, -0.15, 0.15, 0.45)),
  ds_5x6 = make_ds_spec("ds_5x6", 5L, 6L,
                        pi_true = c(0.38, 0.27, 0.18, 0.11, 0.06),
                        accuracy = c(0.86, 0.82, 0.78, 0.74, 0.70, 0.66),
                        severity = c(-0.65, -0.35, -0.10, 0.15, 0.40, 0.70))
)

dgp_truth <- function(spec, weight) {
  L <- 1 - loglin2_weight_matrix(weight, seq_len(spec$C) - 1L)
  loglin2_kappa_from_margins(spec$marg, spec$P, L)
}

# ---- missingness mechanisms (mirroring study 29) -----------------------------

category_z <- function(x, C) {
  z <- (as.numeric(x) - mean(seq_len(C))) / stats::sd(seq_len(C))
  z[!is.finite(z)] <- 0
  z
}
tune_intercept <- function(eta, target) {
  stats::uniroot(function(a) mean(stats::plogis(a + eta)) - target, c(-20, 20))$root
}
drop_with_prob <- function(X, P, keep_first = FALSE) {
  out <- X
  M <- matrix(stats::runif(length(X)), nrow(X), ncol(X)) < P
  if (keep_first && ncol(out) >= 1L) M[, 1L] <- FALSE
  out[M] <- NA_integer_
  out
}
apply_mcar <- function(X, rate) drop_with_prob(X, matrix(rate, nrow(X), ncol(X)))
apply_anchor_mar <- function(X, C, target = 0.30, slope = 1.15) {
  R <- ncol(X)
  eta_base <- slope * category_z(X[, 1L], C)
  P <- matrix(0, nrow(X), R)
  for (j in seq_len(R)) if (j > 1L) {
    a <- tune_intercept(eta_base, target)
    P[, j] <- stats::plogis(a + eta_base)
  }
  drop_with_prob(X, P, keep_first = TRUE)
}
apply_missing <- function(X, spec, mechanism) {
  if (mechanism == "complete") return(X)
  if (mechanism == "mcar30") return(apply_mcar(X, 0.30))
  if (mechanism == "mar_anchor30") return(apply_anchor_mar(X, spec$C, 0.30))
  stop("Unknown mechanism: ", mechanism, call. = FALSE)
}

# ---- estimators --------------------------------------------------------------
# Each returns list(est, se), both named (Conger, Fleiss), with NA pairs on
# failure. loglin2_fiml uses the C++ estimate_kappa_loglin2 (the R EM in
# R/loglin2_fiml.R is kept only for the closed-form truth helpers): it matches
# the prototype to ~1e-5, runs in milliseconds rather than seconds, and carries
# the delta-method SEs added in 1187547. ds_fiml stays point-only here (its SEs
# are validated separately in workbench 37). SEs are read from each fit's vcov()
# when it exposes Conger/Fleiss; coverage and se/sd are computed from them.

coefs <- c("Conger", "Fleiss")
zcrit <- stats::qnorm(0.975)
na_pair <- function() c(Conger = NA_real_, Fleiss = NA_real_)

# Conger/Fleiss SEs from any fit with a named vcov(); NA where absent or non-PSD.
se_from_vcov <- function(fit) {
  out <- na_pair()
  V <- tryCatch(stats::vcov(fit), error = function(e) NULL)
  if (is.null(V) || is.null(rownames(V))) return(out)
  have <- intersect(coefs, rownames(V))
  vals <- diag(V)[have]
  vals[!is.finite(vals) | vals < 0] <- NA_real_
  out[have] <- sqrt(vals)
  out
}

est_pair <- function(fit) {
  cf <- stats::coef(fit)
  c(Conger = unname(cf["Conger"]), Fleiss = unname(cf["Fleiss"]))
}

fit_estimator <- function(estimator, X, spec, weight) {
  vals <- seq_len(spec$C) - 1L
  tryCatch({
    if (estimator == "loglin2_fiml") {
      # The C++ EM converges on the full C^R table to opts$tol; at the default
      # 1e-8 it grinds for seconds on sparse small-n tables resolving tiny cells
      # the estimand never uses. tol=1e-6 matches the R prototype's effective
      # precision (which converged the *margins* at 1e-6) -- ~5x faster at n=40,
      # estimate within ~3e-3 and SE within ~6e-4 of the tight default, no
      # non-convergence; the larger cap covers the sparsest tables.
      f <- misskappa:::estimate_kappa_loglin2(
        X, weight = weight, values = as.numeric(vals),
        em_options = list(tol = 1e-6, max_iter = 5000L))
      list(est = est_pair(f), se = se_from_vcov(f))
    } else if (estimator == "ds_fiml") {
      f <- fit_ds_fiml(X, spec = NULL, weight = weight, se = FALSE)
      list(est = f$coefficients[coefs], se = na_pair())
    } else {
      cpp_weight <- if (weight == "nominal") "nominal" else weight
      f <- suppressWarnings(misskappa::kappa(
        X, estimator = estimator, weight = cpp_weight, values = vals,
        em_options = if (estimator == "cat_fiml")
          list(tol = 1e-8, max_iter = 12000L) else list()))
      list(est = est_pair(f), se = se_from_vcov(f))
    }
  }, error = function(e) list(est = na_pair(), se = na_pair()))
}

estimators_for <- function(weight) {
  if (weight == "nominal") c("ipw", "cat_fiml", "ds_fiml", "loglin2_fiml")
  else c("ipw", "pairwise", "cat_fiml", "ds_fiml", "loglin2_fiml", "nt_fiml")
}

# ---- design grid -------------------------------------------------------------

build_grid <- function() {
  if (smoke) {
    return(expand.grid(dgp = "ds_4x4", mechanism = c("complete", "mcar30"),
                       n = 40L, weight = "nominal", stringsAsFactors = FALSE))
  }
  nominal <- expand.grid(
    dgp = c("loglin2_4x4", "ds_4x4"),
    mechanism = c("complete", "mcar30", "mar_anchor30"),
    n = c(40L, 100L, 400L, 1600L), weight = "nominal",
    stringsAsFactors = FALSE)
  quad <- expand.grid(
    dgp = c("loglin2_4x4", "ds_4x4"),
    mechanism = c("complete", "mcar30", "mar_anchor30"),
    n = c(100L, 400L), weight = "quadratic", stringsAsFactors = FALSE)
  # ds_5x6 (15625-cell table) stays out of the default grid: the prototype EM is
  # ~2 min/fit and slow to converge there. Reachable via --dgps ds_5x6 for a
  # one-off look; larger C,R belongs to the eventual C++ port.
  rbind(nominal, quad)
}

grid <- build_grid()
if (!is.null(opts$dgps)) grid <- grid[grid$dgp %in% opts$dgps, , drop = FALSE]
if (!is.null(opts$mechanisms)) grid <- grid[grid$mechanism %in% opts$mechanisms, , drop = FALSE]
if (!is.null(opts$weights)) grid <- grid[grid$weight %in% opts$weights, , drop = FALSE]
if (!is.null(opts$n_grid)) grid <- grid[grid$n %in% opts$n_grid, , drop = FALSE]
if (is.null(opts$reps)) opts$reps <- if (smoke) 5L else 500L
if (nrow(grid) == 0L) stop("Empty design grid after filtering.", call. = FALSE)
rownames(grid) <- NULL

# ---- run ---------------------------------------------------------------------

rep_rows <- vector("list", nrow(grid))
for (cell in seq_len(nrow(grid))) {
  dgp <- grid$dgp[cell]; mech <- grid$mechanism[cell]
  n <- grid$n[cell]; weight <- grid$weight[cell]
  spec <- dgp_defs[[dgp]]
  ests <- estimators_for(weight)
  t0 <- proc.time()[["elapsed"]]
  reps <- parallel::mclapply(seq_len(opts$reps), function(r) {
    set.seed(opts$seed_base + cell * 1000000L + r)
    X <- apply_missing(spec$simulate(n), spec, mech)
    do.call(rbind, lapply(ests, function(est) {
      est_t0 <- proc.time()[["elapsed"]]
      k <- fit_estimator(est, X, spec, weight)
      data.frame(dgp = dgp, mechanism = mech, n = n, weight = weight, rep = r,
                 estimator = est, coefficient = coefs,
                 estimate = unname(k$est[coefs]), se = unname(k$se[coefs]),
                 elapsed_ms = 1000 * (proc.time()[["elapsed"]] - est_t0),
                 stringsAsFactors = FALSE)
    }))
  }, mc.cores = opts$cores)
  ok <- reps[!vapply(reps, is.null, logical(1)) &
               vapply(reps, is.data.frame, logical(1))]
  rep_rows[[cell]] <- do.call(rbind, ok)
  cat(sprintf("cell %d/%d  %s %s n=%d %s  reps=%d  %.1fs\n",
              cell, nrow(grid), dgp, mech, n, weight, opts$reps,
              proc.time()[["elapsed"]] - t0))
}
replicates <- do.call(rbind, rep_rows)
path_r <- file.path(opts$out_dir, "replicates.csv")
utils::write.csv(replicates, path_r, row.names = FALSE)

# Truth per (dgp, weight, coefficient).
truth_rows <- list()
for (dgp in unique(grid$dgp)) for (weight in unique(grid$weight[grid$dgp == dgp])) {
  k <- dgp_truth(dgp_defs[[dgp]], weight)
  truth_rows[[length(truth_rows) + 1L]] <- data.frame(
    dgp = dgp, weight = weight, coefficient = coefs, truth = unname(k[coefs]),
    stringsAsFactors = FALSE)
}
truths <- do.call(rbind, truth_rows)
path_t <- file.path(opts$out_dir, "truths.csv")
utils::write.csv(truths, path_t, row.names = FALSE)

# Summary: bias / sd / rmse plus SE calibration (mean SE, se/sd ratio, Wald
# coverage about the exact truth) per cell x estimator x coefficient. se/sd
# isolates SE calibration; coverage-about-truth additionally absorbs point bias,
# so low coverage with se/sd ~ 1 is a bias effect (small-n or the H3 MAR
# boundary), not an SE defect. Estimators carried point-only (ds_fiml) get NA SE
# columns.
ok <- replicates[is.finite(replicates$estimate), ]
ok <- merge(ok, truths, by = c("dgp", "weight", "coefficient"))
keys <- c("dgp", "mechanism", "n", "weight", "estimator", "coefficient")
summary <- do.call(rbind, by(ok, ok[keys], function(d) {
  bias <- mean(d$estimate) - d$truth[1L]
  sdv <- stats::sd(d$estimate)
  se_ok <- is.finite(d$se)
  mean_se <- if (any(se_ok)) mean(d$se[se_ok]) else NA_real_
  cover <- if (any(se_ok))
    mean(abs(d$estimate[se_ok] - d$truth[1L]) <= zcrit * d$se[se_ok]) else NA_real_
  data.frame(d[1L, keys, drop = FALSE], truth = d$truth[1L], n_ok = nrow(d),
             bias = bias, sd = sdv, rmse = sqrt(bias^2 + sdv^2),
             mean_se = mean_se, se_sd = mean_se / sdv, coverage = cover,
             n_se = sum(se_ok), med_ms = stats::median(d$elapsed_ms),
             stringsAsFactors = FALSE)
}))
rownames(summary) <- NULL
summary <- summary[order(summary$dgp, summary$weight, summary$mechanism,
                         summary$coefficient, summary$estimator, summary$n), ]
path_s <- file.path(opts$out_dir, "summary.csv")
utils::write.csv(summary, path_s, row.names = FALSE)

n_fail <- sum(!is.finite(replicates$estimate))
meta <- data.frame(
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  args = paste(args, collapse = " "), seed_base = opts$seed_base,
  reps = opts$reps, cores = opts$cores, cells = nrow(grid),
  n_failed_fits = n_fail, r_version = as.character(getRversion()),
  stringsAsFactors = FALSE)
utils::write.csv(meta, file.path(opts$out_dir, "metadata.csv"), row.names = FALSE)

cat("Wrote", path_r, "\n"); cat("Wrote", path_t, "\n")
cat("Wrote", path_s, "\n")
if (n_fail > 0L) cat("Non-finite fits:", n_fail, "\n")
final <- summary[summary$weight == "nominal" & summary$coefficient == "Conger",
                 c("dgp", "mechanism", "n", "estimator", "bias", "sd",
                   "mean_se", "se_sd", "coverage", "med_ms")]
print(final, digits = 3, row.names = FALSE)
