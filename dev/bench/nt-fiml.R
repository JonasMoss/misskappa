#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(misskappa))

reps <- as.integer(Sys.getenv("MISSKAPPA_NTFIML_BENCH_REPS", "3"))
if (!is.finite(reps) || reps < 1L) reps <- 3L

make_matrix <- function(n, p, miss, seed, rho = 0.35) {
  set.seed(seed)
  Sigma <- matrix(rho, p, p)
  diag(Sigma) <- 1
  X <- matrix(stats::rnorm(n * p), n, p) %*% chol(Sigma)
  if (miss > 0) X[matrix(stats::runif(n * p) < miss, n, p)] <- NA_real_
  X
}

make_array <- function(n, raters, features, miss, seed) {
  set.seed(seed)
  X <- array(stats::rnorm(n * raters * features), c(n, raters, features))
  for (r in seq_len(raters)) X[, r, ] <- X[, r, ] + (r - 1) * 0.15
  if (miss > 0) X[array(stats::runif(length(X)) < miss, dim(X))] <- NA_real_
  X
}

bench_one <- function(target, n, p, miss, make_fit) {
  invisible(make_fit())
  timings <- numeric(reps)
  iterations <- integer(reps)
  converged <- logical(reps)
  for (i in seq_len(reps)) {
    gc()
    fit <- NULL
    timing <- system.time(fit <- make_fit())[["elapsed"]]
    timings[[i]] <- timing
    iterations[[i]] <- fit$moments$iterations
    converged[[i]] <- isTRUE(fit$moments$converged)
  }
  data.frame(
    target = target,
    n = n,
    p = p,
    missing = miss,
    reps = reps,
    median_seconds = stats::median(timings),
    min_seconds = min(timings),
    max_seconds = max(timings),
    median_iterations = stats::median(iterations),
    all_converged = all(converged),
    stringsAsFactors = FALSE
  )
}

rows <- list()
case_id <- 1L
for (case in list(
  list(n = 1000L, p = 8L, miss = 0.00),
  list(n = 1000L, p = 8L, miss = 0.15),
  list(n = 1000L, p = 8L, miss = 0.30),
  list(n = 1000L, p = 12L, miss = 0.15),
  list(n = 1000L, p = 12L, miss = 0.30)
)) {
  X <- make_matrix(case$n, case$p, case$miss, seed = 100L + case_id)
  rows[[length(rows) + 1L]] <- bench_one(
    "kappa_quadratic_fiml", case$n, case$p, case$miss,
    function() misskappa:::kappa_quadratic_fiml(X)
  )
  rows[[length(rows) + 1L]] <- bench_one(
    "alpha_continuous", case$n, case$p, case$miss,
    function() misskappa:::alpha_continuous(X)
  )
  case_id <- case_id + 1L
}

Xv <- make_array(1000L, raters = 4L, features = 3L, miss = 0.15, seed = 400L)
W <- matrix(c(
  1.0, 0.2, 0.1,
  0.2, 1.4, 0.3,
  0.1, 0.3, 0.9
), 3L, 3L, byrow = TRUE)
rows[[length(rows) + 1L]] <- bench_one(
  "kappa_vector_quadratic_nt_fiml", 1000L, 12L, 0.15,
  function() misskappa:::kappa_vector_quadratic(Xv, method = "nt_fiml", W = W)
)

out <- do.call(rbind, rows)
print(out, row.names = FALSE)
