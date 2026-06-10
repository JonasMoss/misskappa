# Study 31: calibration of the cat_fiml flattening constant.
#
# Question: how does the total Dirichlet pseudo-mass `flatten` (c) affect
# bias / SD / RMSE / CI coverage of raw cat_fiml kappa, across the sparse
# regimes where the old identification guard fired? Expectation from
# papers/00_wip/kappa-missing/dev/notes/cat-fiml-mle-face-flattening.md:
# everything c <= 1 is indistinguishable at the n^{-1/2} scale, with c = 0.1
# a safe default; c = 0 (strict ML, deterministic-start face point) is the
# control arm.
#
# Design: latent-class DGP with exact closed-form truth. z ~ Uniform(1..C);
# given z each rater independently reports z with prob `acc`, otherwise a
# uniform draw from 1..C. All raters share the uniform margin, so Conger,
# Fleiss, and BP all have the same population value
#   kappa = (po - 1/C) / (1 - 1/C),
#   po    = (acc + (1-acc)/C)^2 + (C-1) ((1-acc)/C)^2.
#
# Mechanisms (both leave every rater pair co-observed w.h.p., so the design
# guard passes and any failure/wobble is the saturated-nuisance story):
#   pairs  - planned missingness, MCAR: each subject keeps one random rater
#            pair (the regime where the old guard misfired).
#   anchor - MAR: rater 1 always observed; each other rater observed with
#            prob 0.9 when x1 <= ceil(C/2), else 0.35.
#
# Common random numbers: one dataset per (dgp, mech, n, rep), all flatten
# values fit on the same data.
#
# Usage: Rscript run_calibration.R [reps] [cores]

suppressMessages(library(misskappa))

args <- commandArgs(trailingOnly = TRUE)
REPS  <- if (length(args) >= 1) as.integer(args[1]) else 1000L
CORES <- if (length(args) >= 2) as.integer(args[2]) else max(1L, parallel::detectCores() - 2L)

FLATTEN <- c(0, 0.01, 0.1, 0.5, 1)
NS      <- c(20L, 40L, 100L)
MECHS   <- c("pairs", "anchor")
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
  seed <- 31000000L + k
  x <- gen_data(cell$n, cfg$C, cfg$R, cfg$acc, cell$mech, seed)
  truth <- true_kappa(cfg$C, cfg$acc)

  out <- vector("list", length(FLATTEN))
  for (fi in seq_along(FLATTEN)) {
    fl <- FLATTEN[fi]
    fit <- tryCatch(
      suppressWarnings(misskappa::kappa(
        x, estimator = "cat_fiml", weight = "nominal",
        em_options = list(flatten = fl)
      )),
      error = function(e) conditionMessage(e)
    )
    if (is.character(fit)) {
      out[[fi]] <- data.frame(
        dgp = cell$dgp, mech = cell$mech, n = cell$n, rep = cell$rep,
        flatten = fl, coef = NA_character_, est = NA_real_, se = NA_real_,
        null_frac = NA_real_, truth = truth, error = substr(fit, 1, 60),
        stringsAsFactors = FALSE
      )
    } else {
      nf <- fit$null_frac
      if (is.null(nf) || length(nf) != 3L) nf <- rep(NA_real_, 3L)
      out[[fi]] <- data.frame(
        dgp = cell$dgp, mech = cell$mech, n = cell$n, rep = cell$rep,
        flatten = fl, coef = names(fit$estimates),
        est = unname(fit$estimates), se = sqrt(diag(fit$vcov)),
        null_frac = unname(nf), truth = truth, error = NA_character_,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, out)
}

message(sprintf("study 31: %d datasets x %d flatten values on %d cores",
                nrow(cells), length(FLATTEN), CORES))
t0 <- Sys.time()
res <- parallel::mclapply(seq_len(nrow(cells)), one_cell,
                          mc.cores = CORES, mc.preschedule = TRUE)
bad <- vapply(res, function(r) inherits(r, "try-error") || is.null(r), logical(1))
if (any(bad)) message(sprintf("dropped %d worker failures", sum(bad)))
res <- do.call(rbind, res[!bad])
message(sprintf("done in %.1f min", as.numeric(Sys.time() - t0, units = "mins")))

dir.create("results", showWarnings = FALSE)
write.csv(res, "results/calibration-raw.csv", row.names = FALSE)
source("summarize.R")
summ <- summarize_calibration(res)
write.csv(summ, "results/calibration-summary.csv", row.names = FALSE)
print(summ[summ$coef %in% c(NA, "Conger"), ], digits = 3, row.names = FALSE)
