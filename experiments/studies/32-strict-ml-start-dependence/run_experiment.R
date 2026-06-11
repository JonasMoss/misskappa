# Study 32: start-dependence of strict-ML cat_fiml on sparse tables.
#
# Question: with flatten = 0 (strict ML), how much do the kappa point
# estimate, its SE, and the null_frac diagnostic move when the EM start
# changes, per dataset? The design guard already pins kappa on the
# population MLE face (note §3 of cat-fiml-mle-face-flattening.md); the
# only leak is finite-sample width through zero-count cells (§4). Study 31
# showed face-point *selection* is statistically free in aggregate (RMSE
# ratio <= 1.005 vs the flattened analytic center); this study measures the
# per-dataset spread directly, including the SE, which study 30 did not.
#
# Start perturbation: start_alpha in {0.01, 0.1, 1, 10}. This is the only
# start perturbation the C++ supports — initialise_theta() is deterministic
# (theta_j propto start_alpha + sum_g n_subs_g / n_comps_g), with no RNG.
#
# Grid: study-31 DGPs and mechanisms, n in {20, 40} only (face width is
# ~1e-12 by n = 160; study 30 script 04).
#
# Pre-registered reading (in README.md): if the fraction of datasets where
# the start moves the estimate by > 0.1 SE or the SE by > 10% is negligible
# (< 1-2%) in every cell, start-dependence does not justify changing the
# package default away from strict ML.
#
# Usage: Rscript run_experiment.R [reps] [cores]

suppressMessages(library(misskappa))

args <- commandArgs(trailingOnly = TRUE)
REPS  <- if (length(args) >= 1) as.integer(args[1]) else 1000L
CORES <- if (length(args) >= 2) as.integer(args[2]) else max(1L, parallel::detectCores() - 2L)

STARTS <- c(0.01, 0.1, 1, 10)
NS     <- c(20L, 40L)
MECHS  <- c("pairs", "anchor")
DGPS <- list(
  lc3x3 = list(C = 3L, R = 3L, acc = 0.75),
  lc4x4 = list(C = 4L, R = 4L, acc = 0.70),
  lc5x4 = list(C = 5L, R = 4L, acc = 0.70)
)

true_kappa <- function(C, acc) {
  po <- (acc + (1 - acc) / C)^2 + (C - 1) * ((1 - acc) / C)^2
  (po - 1 / C) / (1 - 1 / C)
}

gen_data <- function(n, C, R, acc, mech, seed) {
  set.seed(seed)
  z <- sample.int(C, n, replace = TRUE)
  x <- vapply(seq_len(R), function(j) {
    ifelse(runif(n) < acc, z, sample.int(C, n, replace = TRUE))
  }, integer(n))
  if (mech == "pairs") {
    prs <- t(combn(R, 2))
    for (i in seq_len(n)) {
      keep <- prs[sample.int(nrow(prs), 1L), ]
      x[i, setdiff(seq_len(R), keep)] <- NA_integer_
    }
  } else {
    p_obs <- ifelse(x[, 1] <= ceiling(C / 2), 0.9, 0.35)
    for (j in 2:R) x[runif(n) >= p_obs, j] <- NA_integer_
  }
  x
}

cells <- expand.grid(dgp = names(DGPS), mech = MECHS, n = NS,
                     rep = seq_len(REPS), stringsAsFactors = FALSE)

one_cell <- function(k) {
  cell <- cells[k, ]
  cfg <- DGPS[[cell$dgp]]
  seed <- 32000000L + k
  x <- gen_data(cell$n, cfg$C, cfg$R, cfg$acc, cell$mech, seed)
  truth <- true_kappa(cfg$C, cfg$acc)

  out <- vector("list", length(STARTS))
  for (si in seq_along(STARTS)) {
    sa <- STARTS[si]
    fit <- tryCatch(
      suppressWarnings(misskappa::kappa(
        x, estimator = "cat_fiml", weight = "nominal",
        em_options = list(flatten = 0, start_alpha = sa)
      )),
      error = function(e) conditionMessage(e)
    )
    if (is.character(fit)) {
      out[[si]] <- data.frame(
        dgp = cell$dgp, mech = cell$mech, n = cell$n, rep = cell$rep,
        start_alpha = sa, coef = NA_character_, est = NA_real_, se = NA_real_,
        null_frac = NA_real_, truth = truth, error = substr(fit, 1, 60),
        stringsAsFactors = FALSE
      )
    } else {
      nf <- fit$null_frac
      if (is.null(nf) || length(nf) != 3L) nf <- rep(NA_real_, 3L)
      out[[si]] <- data.frame(
        dgp = cell$dgp, mech = cell$mech, n = cell$n, rep = cell$rep,
        start_alpha = sa, coef = names(fit$estimates),
        est = unname(fit$estimates), se = sqrt(diag(fit$vcov)),
        null_frac = unname(nf), truth = truth, error = NA_character_,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, out)
}

message(sprintf("study 32: %d datasets x %d starts on %d cores",
                nrow(cells), length(STARTS), CORES))
t0 <- Sys.time()
res <- parallel::mclapply(seq_len(nrow(cells)), one_cell,
                          mc.cores = CORES, mc.preschedule = TRUE)
bad <- vapply(res, function(r) inherits(r, "try-error") || is.null(r), logical(1))
if (any(bad)) message(sprintf("dropped %d worker failures", sum(bad)))
res <- do.call(rbind, res[!bad])
message(sprintf("done in %.1f min", as.numeric(Sys.time() - t0, units = "mins")))

dir.create("results", showWarnings = FALSE)
write.csv(res, "results/start-spread-raw.csv", row.names = FALSE)
source("summarize.R")
summ <- summarize_starts(res)
write.csv(summ, "results/start-spread-summary.csv", row.names = FALSE)
print(summ[summ$coef %in% c(NA, "Conger"), ], digits = 3, row.names = FALSE)
