#!/usr/bin/env Rscript

# Probe 46: cat-FIML feasibility frontier.
#
# Maps where saturated categorical FIML (`alpha(estimator = "cat_fiml")`)
# becomes infeasible, separating two distinct walls that earlier notes
# conflated:
#
#   1. SINGLE-FIT wall  -- one EM fit takes too long / degenerates. Driven by
#      the cell count C^p and (super-linearly) by n.
#   2. SIM-THROUGHPUT wall -- single-fit-time x B replications. ~1000x lower
#      than the single-fit wall, and the one that actually bit probe 16
#      (its 5^6, n=4000 "stall" was B=200 x ~12 s/fit ~= 40 min for the cell,
#      not one fit).
#
# Per fit we record wall-time, the null-space degeneracy fraction `null_frac`
# (the statistical feasibility axis: start-dependence + pseudo-inverse SE),
# whether the SE is finite, and the bias against a large-sample reference
# alpha. A per-fit elapsed guard turns "too slow" into a recorded datum
# (status = "timeout") rather than a hang.

suppressMessages(library(misskappa))
options(warn = -1)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[[i + 1L]]
}

this_file <- sub("^--file=", "",
                 commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1L])
script_dir <- dirname(normalizePath(this_file))

reps      <- as.integer(get_arg("--reps", "3"))
seed_base <- as.integer(get_arg("--seed-base", "20260627"))
max_cells <- as.numeric(get_arg("--max-cells", "200000"))
fit_tl    <- as.numeric(get_arg("--fit-timeout", "60"))
sim_B     <- as.numeric(get_arg("--sim-B", "1000"))
out_dir   <- get_arg("--out-dir", file.path(script_dir, "results"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

Cs    <- as.integer(strsplit(get_arg("--cats", "3,5,7"), ",")[[1L]])
ps    <- as.integer(strsplit(get_arg("--items", "4,6,8"), ",")[[1L]])
ns    <- as.integer(strsplit(get_arg("--n", "100,250,500,1000,2000,4000"), ",")[[1L]])
mechs <- strsplit(get_arg("--mechanisms", "complete,mar30"), ",")[[1L]]

# One-factor latent Gaussian thresholded to C equiprobable categories; item 1
# is the fully observed MAR anchor. lambda is the common loading.
gen <- function(n, p, C, mech, seed) {
  set.seed(seed)
  lambda <- 0.65
  z <- matrix(stats::rnorm(n), n, 1L)
  Y <- z %*% matrix(lambda, 1L, p) +
    matrix(stats::rnorm(n * p, sd = sqrt(1 - lambda^2)), n, p)
  th <- stats::qnorm(cumsum(rep(1, C) / C)[-C])
  X <- matrix(as.integer(apply(Y, 2L, function(col) findInterval(col, th))), n, p)
  if (mech == "mar30") {
    a1 <- as.numeric(scale(X[, 1L]))
    a1[!is.finite(a1)] <- 0
    for (j in 2:p) {
      intercept <- stats::uniroot(
        function(a) mean(stats::plogis(a - 1.35 * a1)) - 0.30,
        c(-12, 12))$root
      X[stats::runif(n) < stats::plogis(intercept - 1.35 * a1), j] <- NA_integer_
    }
  }
  X
}

# Large-sample reference alpha per (p, C), via the (cheap, structure-free)
# pairwise estimator on a big complete sample.
truth_cache <- new.env(parent = emptyenv())
truth_alpha <- function(p, C) {
  key <- paste(p, C, sep = "x")
  if (exists(key, truth_cache, inherits = FALSE)) {
    return(get(key, truth_cache, inherits = FALSE))
  }
  X <- gen(200000L, p, C, "complete", 987654L)
  a <- as.numeric(stats::coef(alpha(X, estimator = "pairwise"))[["alpha"]])
  assign(key, a, truth_cache)
  a
}

rows <- list()
idx <- 0L
total <- length(Cs) * length(ps) * length(ns) * length(mechs) * reps
for (C in Cs) for (p in ps) {
  cells <- C^p
  for (n in ns) for (mech in mechs) for (r in seq_len(reps)) {
    idx <- idx + 1L
    base <- data.frame(C = C, p = p, cells = cells, n = n, mechanism = mech,
                       rep = r, obs_patterns = NA_integer_, elapsed_s = NA_real_,
                       null_frac = NA_real_, se_finite = NA, est = NA_real_,
                       bias = NA_real_, status = NA_character_,
                       stringsAsFactors = FALSE)
    if (cells > max_cells) {
      base$status <- "skipped_too_large"
      rows[[idx]] <- base
      next
    }
    X <- gen(n, p, C, mech, seed_base + 7919L * idx)
    base$obs_patterns <- nrow(unique(X))
    tru <- truth_alpha(p, C)
    setTimeLimit(elapsed = fit_tl, transient = TRUE)
    t0 <- proc.time()[["elapsed"]]
    fit <- tryCatch(
      alpha(X, estimator = "cat_fiml",
            em_options = list(tol = 1e-7, max_iter = 20000L, info_rcond = 1e-4)),
      error = function(e) e)
    base$elapsed_s <- proc.time()[["elapsed"]] - t0
    setTimeLimit()
    if (inherits(fit, "error")) {
      base$status <- if (grepl("elapsed time limit", conditionMessage(fit))) {
        "timeout"
      } else {
        "error"
      }
      base$se_finite <- FALSE
    } else {
      est <- as.numeric(stats::coef(fit)[["alpha"]])
      se <- sqrt(diag(stats::vcov(fit)))[1L]
      base$est <- est
      base$bias <- est - tru
      base$se_finite <- is.finite(se) && se > 0
      base$null_frac <- if (!is.null(fit$null_frac)) as.numeric(fit$null_frac) else NA_real_
      base$status <- "ok"
    }
    rows[[idx]] <- base
    cat(sprintf("[%3d/%3d] C=%d p=%d cells=%d n=%d %-8s rep=%d  %.2fs  %s\n",
                idx, total, C, p, cells, n, mech, r, base$elapsed_s, base$status))
    flush.console()
  }
}

reps_df <- do.call(rbind, rows)
write.csv(reps_df, file.path(out_dir, "replicates.csv"), row.names = FALSE)

# Per-cell summary + derived feasibility verdicts.
key <- with(reps_df, paste(C, p, cells, n, mechanism, sep = "|"))
summ <- do.call(rbind, lapply(split(reps_df, key), function(g) {
  ok <- g$status == "ok"
  med_t <- stats::median(g$elapsed_s, na.rm = TRUE)
  data.frame(
    C = g$C[1L], p = g$p[1L], cells = g$cells[1L], n = g$n[1L],
    mechanism = g$mechanism[1L],
    reps = nrow(g), n_ok = sum(ok),
    median_s = med_t,
    max_null_frac = if (any(ok)) max(g$null_frac[ok], na.rm = TRUE) else NA_real_,
    frac_se_finite = if (any(ok)) mean(g$se_finite[ok]) else 0,
    mean_bias = if (any(ok)) mean(g$bias[ok], na.rm = TRUE) else NA_real_,
    status = if (all(g$status == "skipped_too_large")) "skipped_too_large"
             else if (any(g$status == "timeout")) "timeout"
             else if (any(g$status == "error")) "error" else "ok",
    single_fit_feasible = isTRUE(med_t < fit_tl) && any(ok),
    sim_feasible = isTRUE(med_t * sim_B < 1800) && any(ok),  # B reps < 30 min/cell
    stringsAsFactors = FALSE)
}))
summ <- summ[order(summ$C, summ$p, summ$n, summ$mechanism), ]
write.csv(summ, file.path(out_dir, "summary.csv"), row.names = FALSE)

cat("\nWrote ", file.path(out_dir, "replicates.csv"), " and summary.csv\n", sep = "")
