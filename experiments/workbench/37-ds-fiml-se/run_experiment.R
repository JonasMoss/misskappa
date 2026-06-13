#!/usr/bin/env Rscript
#
# Experiment 36: are the new Dawid-Skene FIML standard errors correct?
#
# The DS-FIML prototype (study 29, R/ds_fiml.R) gained delta-method inference:
# analytic observed-data scores in softmax coordinates, FD observed
# information, rank-truncated pseudo-inverse, and ML / robust-sandwich vcovs
# for the model-implied Conger and Fleiss kappas. This probe validates the
# machinery (score vs numerical gradient, information vs numerical Hessian,
# bootstrap cross-check) and calibrates the SEs by simulation: se/sd ratios
# and Wald coverage under the correctly specified ds5x6_severity DGP
# (complete / MCAR / anchor-MAR) and under copula misspecification, where
# coverage is taken around the mechanism-specific pseudo-true value.

usage <- function(status = 0L) {
  cat(paste0(
    "Usage: Rscript run_experiment.R [mode/options]\n\n",
    "Modes:\n",
    "  --checks      Correctness checks only (score, Hessian, bootstrap).\n",
    "  --calibrate   Simulation calibration of the SEs (default reps: 500).\n",
    "  --smoke       Cheap end-to-end pass: checks with a small bootstrap,\n",
    "                calibration with reps=3 on one cell.\n\n",
    "Options:\n",
    "  --reps N        Replications per calibration cell.\n",
    "  --boot-b N      Bootstrap replicates for the cross-check.\n",
    "  --seed-base N   Base seed. Default: 360000.\n",
    "  --cores N       Worker processes. Default: min(8, detected - 1).\n",
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

source(file.path(script_dir, "..", "..", "studies",
                 "29-kappa-planned-missingness", "R", "ds_fiml.R"))

opts <- list(
  mode = NULL, reps = NULL, boot_b = NULL, seed_base = 360000L,
  cores = max(1L, min(8L, parallel::detectCores() - 1L)),
  out_dir = file.path(script_dir, "results")
)
args <- commandArgs(TRUE)
i <- 1L
while (i <= length(args)) {
  arg <- args[[i]]
  if (arg %in% c("--help", "-h")) usage()
  if (arg %in% c("--checks", "--calibrate", "--smoke")) {
    opts$mode <- c(opts$mode, sub("^--", "", arg))
    i <- i + 1L
    next
  }
  if (!arg %in% c("--reps", "--boot-b", "--seed-base", "--cores", "--out-dir")) {
    stop("Unknown argument: ", arg, " (see --help)", call. = FALSE)
  }
  if (i == length(args)) stop(arg, " needs a value.", call. = FALSE)
  val <- args[[i + 1L]]
  if (arg == "--reps") opts$reps <- as.integer(val)
  if (arg == "--boot-b") opts$boot_b <- as.integer(val)
  if (arg == "--seed-base") opts$seed_base <- as.integer(val)
  if (arg == "--cores") opts$cores <- as.integer(val)
  if (arg == "--out-dir") opts$out_dir <- val
  i <- i + 2L
}
if (is.null(opts$mode)) opts$mode <- c("checks", "calibrate")
smoke <- "smoke" %in% opts$mode
if (smoke) opts$mode <- c("checks", "calibrate")
if (is.null(opts$reps)) opts$reps <- if (smoke) 3L else 500L
if (is.null(opts$boot_b)) opts$boot_b <- if (smoke) 20L else 300L
dir.create(opts$out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- DGPs and mechanisms (mirroring study 29 run_experiment.R) ---------------

ds5x6_severity <- list(
  id = "ds5x6_severity", C = 5L, R = 6L,
  p = c(0.38, 0.27, 0.18, 0.11, 0.06),
  accuracy = c(0.86, 0.82, 0.78, 0.74, 0.70, 0.66),
  severity = c(-0.65, -0.35, -0.10, 0.15, 0.40, 0.70)
)

copula5x6_high <- list(
  id = "copula5x6_high", C = 5L, R = 6L,
  p = rep(0.2, 5L), latent_corr = 0.72
)

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

ds_true_params <- function(spec) {
  theta <- lapply(seq_len(spec$R), function(j) {
    t(vapply(seq_len(spec$C), function(k) {
      ds_row_prob(spec$C, k, spec$accuracy[j], spec$severity[j])
    }, numeric(spec$C)))
  })
  list(pi = spec$p, theta = theta)
}

simulate_ds <- function(n, spec) {
  par <- ds_true_params(spec)
  truth <- sample.int(spec$C, n, replace = TRUE, prob = spec$p)
  X <- matrix(NA_integer_, nrow = n, ncol = spec$R)
  for (j in seq_len(spec$R)) {
    for (k in seq_len(spec$C)) {
      idx <- which(truth == k)
      if (length(idx)) {
        X[idx, j] <- sample.int(spec$C, length(idx), replace = TRUE,
                                prob = par$theta[[j]][k, ])
      }
    }
  }
  X
}

simulate_copula <- function(n, spec) {
  Sigma <- matrix(spec$latent_corr, spec$R, spec$R)
  diag(Sigma) <- 1
  Z <- matrix(stats::rnorm(n * spec$R), nrow = n, ncol = spec$R) %*% chol(Sigma)
  cuts <- c(-Inf, stats::qnorm(cumsum(spec$p)[-spec$C]), Inf)
  X <- apply(Z, 2L, function(z) as.integer(cut(z, breaks = cuts, labels = FALSE)))
  matrix(X, nrow = n, ncol = spec$R)
}

simulate_ratings <- function(n, spec) {
  if (spec$id == "copula5x6_high") simulate_copula(n, spec) else simulate_ds(n, spec)
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

apply_mcar <- function(X, rate) drop_with_prob(X, matrix(rate, nrow(X), ncol(X)))

apply_anchor_mar <- function(X, C, target = 0.30, slope = 1.15) {
  n <- nrow(X)
  R <- ncol(X)
  eta_base <- slope * category_z(X[, 1L], C)
  P <- matrix(0, n, R)
  for (j in seq_len(R)) {
    if (j == 1L) next
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

weights <- c("nominal", "quadratic")
coefs <- c("Conger", "Fleiss")
zcrit <- stats::qnorm(0.975)

ds_L <- function(weight, C) 1 - ds_weight_matrix(weight, seq_len(C) - 1)

# One EM fit, both weights served from the shared likelihood cache.
fit_both_weights <- function(X, C, se = TRUE) {
  Y <- matrix(match(X, seq_len(C)), nrow(X), ncol(X))
  if (length(unique(Y[!is.na(Y)])) < C) {
    stop("A category is unobserved in this sample.", call. = FALSE)
  }
  fit <- ds_em_fit(Y, C = C)
  cache <- if (se) ds_vcov_cache(fit, Y) else NULL
  out <- lapply(weights, function(w) {
    L <- ds_L(w, C)
    est <- ds_implied_kappa(fit, L)
    if (!se) return(list(est = est))
    vc <- ds_kappa_vcov(cache, L)
    list(est = est, se_ml = sqrt(pmax(diag(vc$vcov_ml), 0)),
         se_rb = sqrt(pmax(diag(vc$vcov_robust), 0)), null_frac = vc$null_frac)
  })
  names(out) <- weights
  out
}

# ---- correctness checks ------------------------------------------------------

run_checks <- function() {
  rows <- list()
  add <- function(check, statistic, value, tol, pass) {
    rows[[length(rows) + 1L]] <<- data.frame(
      check = check, statistic = statistic, value = value,
      tolerance = tol, pass = pass, stringsAsFactors = FALSE)
  }

  # 1. Analytic score equals the numerical gradient of the observed-data
  #    log-likelihood, at a non-MLE point of a small holed dataset.
  set.seed(opts$seed_base + 1L)
  spec_small <- list(id = "ds_small", C = 3L, R = 4L,
                     p = c(0.5, 0.3, 0.2),
                     accuracy = c(0.8, 0.75, 0.7, 0.65),
                     severity = c(-0.3, 0, 0.2, 0.5))
  Xs <- apply_mcar(simulate_ds(30L, spec_small), 0.3)
  Ys <- matrix(match(Xs, seq_len(3L)), nrow(Xs), ncol(Xs))
  K <- 3L; C <- 3L; R <- 4L
  p_dim <- K + R * K * C
  eta0 <- stats::rnorm(p_dim, sd = 0.5)
  s_analytic <- ds_score_total_eta(eta0, Ys, K, C, R)
  h <- 1e-6
  s_fd <- vapply(seq_len(p_dim), function(r) {
    e <- numeric(p_dim); e[r] <- h
    (ds_loglik_eta(eta0 + e, Ys, K, C, R) -
       ds_loglik_eta(eta0 - e, Ys, K, C, R)) / (2 * h)
  }, numeric(1))
  err <- max(abs(s_analytic - s_fd)) / (1 + max(abs(s_fd)))
  add("score_vs_fd_gradient", "max_rel_err", err, 1e-5, err < 1e-5)

  # 2. FD-of-score observed information agrees with second differences of the
  #    log-likelihood on a random coordinate block.
  info <- ds_observed_info(eta0, Ys, K, C, R)
  set.seed(opts$seed_base + 2L)
  idx <- sample.int(p_dim, 6L)
  h2 <- 1e-4
  worst <- 0
  for (r in idx) {
    for (s in idx) {
      er <- numeric(p_dim); er[r] <- h2
      es <- numeric(p_dim); es[s] <- h2
      hrs <- -(ds_loglik_eta(eta0 + er + es, Ys, K, C, R) -
                 ds_loglik_eta(eta0 + er - es, Ys, K, C, R) -
                 ds_loglik_eta(eta0 - er + es, Ys, K, C, R) +
                 ds_loglik_eta(eta0 - er - es, Ys, K, C, R)) / (4 * h2^2)
      worst <- max(worst, abs(hrs - info[r, s]) / (1 + abs(hrs)))
    }
  }
  add("info_vs_fd2_loglik", "max_rel_err", worst, 1e-3, worst < 1e-3)

  # 3. At an MLE on the target DGP: total score ~ 0, the truncated null space
  #    is at least as large as the softmax gauge (1 + R*K), and the kappa
  #    gradient barely overlaps it.
  set.seed(opts$seed_base + 3L)
  X <- apply_mcar(simulate_ds(150L, ds5x6_severity), 0.3)
  Y <- matrix(match(X, seq_len(5L)), nrow(X), ncol(X))
  fit <- ds_em_fit(Y, C = 5L)
  cache <- ds_vcov_cache(fit, Y)
  score_norm <- max(abs(colSums(cache$S))) / nrow(Y)
  add("mle_score_zero", "max_abs_mean_score", score_norm, 1e-4, score_norm < 1e-4)
  n_null <- sum(!cache$keep)
  gauge <- 1L + 6L * 5L
  add("gauge_null_space", "n_truncated_minus_gauge", n_null - gauge, NA,
      n_null >= gauge)
  vc <- ds_kappa_vcov(cache, ds_L("quadratic", 5L))
  nf <- max(vc$null_frac)
  add("kappa_gradient_null_frac", "max_null_frac", nf, 1e-4, nf < 1e-4)

  # 4. Bootstrap cross-check: delta-method SEs against the nonparametric
  #    bootstrap SD on one MCAR dataset. Loose band: both estimate the same
  #    asymptotic quantity; B and n keep either side from being exact.
  set.seed(opts$seed_base + 4L)
  Xb <- apply_mcar(simulate_ds(200L, ds5x6_severity), 0.3)
  fit_se <- fit_both_weights(Xb, 5L, se = TRUE)
  boot <- parallel::mclapply(seq_len(opts$boot_b), function(b) {
    set.seed(opts$seed_base + 1000L + b)
    Xr <- Xb[sample.int(nrow(Xb), replace = TRUE), , drop = FALSE]
    tryCatch(
      vapply(fit_both_weights(Xr, 5L, se = FALSE),
             function(f) f$est, numeric(2)),
      error = function(e) matrix(NA_real_, 2, length(weights),
                                 dimnames = list(coefs, weights)))
  }, mc.cores = opts$cores)
  for (w in weights) {
    bmat <- vapply(boot, function(b) b[, w], numeric(2))
    bsd <- apply(bmat, 1L, stats::sd, na.rm = TRUE)
    for (k in seq_along(coefs)) {
      for (flavour in c("se_ml", "se_rb")) {
        ratio <- fit_se[[w]][[flavour]][k] / bsd[k]
        add(paste("bootstrap", w, coefs[k], flavour, sep = "_"),
            "se_over_boot_sd", ratio, 0.25, abs(ratio - 1) < 0.25)
      }
    }
  }

  out <- do.call(rbind, rows)
  path <- file.path(opts$out_dir, "checks.csv")
  utils::write.csv(out, path, row.names = FALSE)
  cat("Wrote", path, "\n")
  print(out, digits = 3)
  invisible(out)
}

# ---- calibration -------------------------------------------------------------

calibration_cells <- function() {
  cells <- expand.grid(
    dgp = "ds5x6_severity", mechanism = c("complete", "mcar30", "mar_anchor30"),
    n = c(40L, 100L, 400L), stringsAsFactors = FALSE)
  mis <- expand.grid(
    dgp = "copula5x6_high", mechanism = "mcar30",
    n = c(100L, 400L), stringsAsFactors = FALSE)
  out <- rbind(cells, mis)
  if (smoke) out <- out[out$dgp == "ds5x6_severity" & out$mechanism == "mcar30" &
                          out$n == 40L, , drop = FALSE]
  out
}

spec_for <- function(id) if (id == "copula5x6_high") copula5x6_high else ds5x6_severity

# Truth: exact implied kappa of the true parameters on the correctly specified
# DS DGP; mechanism-specific pseudo-true value (one large fit on data passed
# through the same mechanism) under copula misspecification.
truth_table <- function(cells) {
  out <- list()
  for (dgp in unique(cells$dgp)) {
    spec <- spec_for(dgp)
    if (dgp == "ds5x6_severity") {
      par <- ds_true_params(spec)
      for (w in weights) {
        kap <- ds_implied_kappa(par, ds_L(w, spec$C))
        out[[length(out) + 1L]] <- data.frame(
          dgp = dgp, mechanism = unique(cells$mechanism[cells$dgp == dgp]),
          weight = w, Conger = kap["Conger"], Fleiss = kap["Fleiss"],
          truth_kind = "exact", stringsAsFactors = FALSE)
      }
    } else {
      for (mech in unique(cells$mechanism[cells$dgp == dgp])) {
        set.seed(opts$seed_base + 9L)
        n_big <- if (smoke) 5000L else 300000L
        Xb <- apply_missing(simulate_ratings(n_big, spec), spec, mech)
        Yb <- matrix(match(Xb, seq_len(spec$C)), nrow(Xb), ncol(Xb))
        fit <- ds_em_fit(Yb, C = spec$C, n_starts = 2L)
        for (w in weights) {
          kap <- ds_implied_kappa(fit, ds_L(w, spec$C))
          out[[length(out) + 1L]] <- data.frame(
            dgp = dgp, mechanism = mech, weight = w,
            Conger = kap["Conger"], Fleiss = kap["Fleiss"],
            truth_kind = "pseudo_true_mc", stringsAsFactors = FALSE)
        }
      }
    }
  }
  do.call(rbind, out)
}

run_calibration <- function() {
  cells <- calibration_cells()
  truths <- truth_table(cells)
  path_t <- file.path(opts$out_dir, "truths.csv")
  utils::write.csv(truths, path_t, row.names = FALSE)

  rep_rows <- list()
  for (cell in seq_len(nrow(cells))) {
    dgp <- cells$dgp[cell]
    mech <- cells$mechanism[cell]
    n <- cells$n[cell]
    spec <- spec_for(dgp)
    t0 <- proc.time()[["elapsed"]]
    reps <- parallel::mclapply(seq_len(opts$reps), function(r) {
      set.seed(opts$seed_base + cell * 1000000L + r)
      X <- apply_missing(simulate_ratings(n, spec), spec, mech)
      res <- tryCatch(fit_both_weights(X, spec$C, se = TRUE),
                      error = function(e) conditionMessage(e))
      if (is.character(res)) {
        return(data.frame(dgp = dgp, mechanism = mech, n = n, rep = r,
                          weight = NA_character_, coefficient = NA_character_,
                          estimate = NA_real_, se_ml = NA_real_,
                          se_robust = NA_real_, null_frac = NA_real_,
                          error = res, stringsAsFactors = FALSE))
      }
      do.call(rbind, lapply(weights, function(w) {
        data.frame(dgp = dgp, mechanism = mech, n = n, rep = r, weight = w,
                   coefficient = coefs, estimate = unname(res[[w]]$est),
                   se_ml = unname(res[[w]]$se_ml),
                   se_robust = unname(res[[w]]$se_rb),
                   null_frac = unname(res[[w]]$null_frac),
                   error = "", stringsAsFactors = FALSE)
      }))
    }, mc.cores = opts$cores)
    rep_rows[[cell]] <- do.call(rbind, reps)
    cat(sprintf("cell %d/%d  %s %s n=%d  reps=%d  %.1fs\n",
                cell, nrow(cells), dgp, mech, n, opts$reps,
                proc.time()[["elapsed"]] - t0))
  }
  replicates <- do.call(rbind, rep_rows)
  path_r <- file.path(opts$out_dir, "replicates.csv")
  utils::write.csv(replicates, path_r, row.names = FALSE)

  ok <- replicates[replicates$error == "" & !is.na(replicates$estimate), ]
  ok <- merge(ok, truths, by = c("dgp", "mechanism", "weight"))
  ok$truth <- ifelse(ok$coefficient == "Conger", ok$Conger, ok$Fleiss)
  ok$cover_ml <- abs(ok$estimate - ok$truth) <= zcrit * ok$se_ml
  ok$cover_rb <- abs(ok$estimate - ok$truth) <= zcrit * ok$se_robust
  agg <- function(d) {
    data.frame(
      n_ok = nrow(d), truth = d$truth[1L],
      bias = mean(d$estimate) - d$truth[1L],
      sd = stats::sd(d$estimate),
      mean_se_ml = mean(d$se_ml), mean_se_robust = mean(d$se_robust),
      ratio_ml = mean(d$se_ml) / stats::sd(d$estimate),
      ratio_robust = mean(d$se_robust) / stats::sd(d$estimate),
      cover_ml = mean(d$cover_ml), cover_robust = mean(d$cover_rb),
      mean_null_frac = mean(d$null_frac), stringsAsFactors = FALSE)
  }
  keys <- c("dgp", "mechanism", "n", "weight", "coefficient")
  summary <- do.call(rbind, by(ok, ok[keys], function(d) {
    cbind(d[1L, keys, drop = FALSE], agg(d))
  }))
  rownames(summary) <- NULL
  summary <- summary[order(summary$dgp, summary$mechanism, summary$n,
                           summary$weight, summary$coefficient), ]
  path_s <- file.path(opts$out_dir, "summary.csv")
  utils::write.csv(summary, path_s, row.names = FALSE)

  n_err <- sum(replicates$error != "")
  meta <- data.frame(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    args = paste(args, collapse = " "),
    seed_base = opts$seed_base, reps = opts$reps, boot_b = opts$boot_b,
    cores = opts$cores, n_error_reps = n_err,
    r_version = as.character(getRversion()), stringsAsFactors = FALSE)
  path_m <- file.path(opts$out_dir, "metadata.csv")
  utils::write.csv(meta, path_m, row.names = FALSE)
  cat("Wrote", path_t, "\n")
  cat("Wrote", path_r, "\n")
  cat("Wrote", path_s, "\n")
  cat("Wrote", path_m, "\n")
  if (n_err > 0L) cat("Replicates with errors:", n_err, "\n")
  print(summary, digits = 3)
  invisible(summary)
}

if ("checks" %in% opts$mode) run_checks()
if ("calibrate" %in% opts$mode) run_calibration()
