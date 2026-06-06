#!/usr/bin/env Rscript
#
# Experiment 23: paper simulation for the quadratic (Conger / Fleiss) kappas
# under missing data. One question: how do listwise deletion, the
# pairwise-available moment estimator, and robust normal-theory FIML compare
# across data types, missingness mechanisms, sample sizes, and rater counts?
#
# Estimators (all run on the numeric / scored rating matrix):
#   listwise  -- complete-case rows only, then the available-case estimator
#   pairwise  -- pairwise-available moment estimator (kappa_continuous available)
#   nt_fiml   -- robust normal-theory FIML (saturated EM + sandwich SE)
#
# Data-generating processes:
#   normal      -- multivariate normal ratings
#   skewed      -- standardized chi-square noise through a Cholesky factor; the
#                  conditional means E[X_j | X_1] are still LINEAR, so FIML stays
#                  consistent under MAR (it is the missingness model that holds)
#   categorical -- judge-skill model: each rater reports the true category with a
#                  skill probability or guesses from a rater-specific
#                  distribution. The implied E[X_j | X_1] is NONLINEAR, so this is
#                  the honest stress for FIML: its working normal conditional mean
#                  is misspecified and the MAR consistency can break.
#
# Note on "nonlinear MAR": FIML is consistent under *any* MAR mechanism when its
# model is correct, so making the missingness probability a nonlinear function of
# the observed data does not stress it. What stresses FIML is a misspecified
# conditional mean, which the categorical DGP supplies; hence the stress lives in
# the DGP, not in a separate nonlinear-MAR mechanism.

suppressPackageStartupMessages({
  library(misskappa)
})

usage <- function(status = 0L) {
  cat(paste0(
    "Usage: Rscript run_experiment.R [options]\n\n",
    "Options:\n",
    "  --help, -h         Show this help and exit.\n",
    "  --smoke            Cheap check: reps=40.\n",
    "  --reps N           Monte Carlo replications per cell. Default: 1000.\n",
    "  --n N              Reference subjects per replication. Default: 600.\n",
    "  --rate R           Per-entry MCAR / target MAR deletion rate. Default: 0.25.\n",
    "  --seed-base N      Base seed. Default: 23000.\n",
    "  --out-dir PATH     Output directory. Default: script-local results/.\n"
  ))
  quit(save = "no", status = status)
}

script_dir <- local({
  cmd <- commandArgs(FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) == 0L) getwd()
  else dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE))
})

opts <- list(reps = 1000L, n = 600L, rate = 0.25, seed_base = 23000L,
             out_dir = file.path(script_dir, "results"))
args <- commandArgs(TRUE)
i <- 1L
while (i <= length(args)) {
  a <- args[[i]]
  if (a %in% c("--help", "-h")) usage(0L)
  else if (a == "--smoke") { opts$reps <- 40L }
  else if (a == "--reps") { opts$reps <- as.integer(args[[i + 1L]]); i <- i + 1L }
  else if (a == "--n") { opts$n <- as.integer(args[[i + 1L]]); i <- i + 1L }
  else if (a == "--rate") { opts$rate <- as.numeric(args[[i + 1L]]); i <- i + 1L }
  else if (a == "--seed-base") { opts$seed_base <- as.integer(args[[i + 1L]]); i <- i + 1L }
  else if (a == "--out-dir") { opts$out_dir <- args[[i + 1L]]; i <- i + 1L }
  else { cat("Unknown argument:", a, "\n"); usage(1L) }
  i <- i + 1L
}

# ---- data-generating processes ---------------------------------------------
make_Sigma <- function(R, rho = 0.6) {
  S <- matrix(rho, R, R); diag(S) <- 1
  d <- sqrt(seq(0.8, 1.4, length.out = R))
  outer(d, d) * S
}
make_mu <- function(R) seq(-0.45, 0.45, length.out = R)

gen_continuous <- function(n, R, dist) {
  Sigma <- make_Sigma(R); mu <- make_mu(R)
  Z <- if (dist == "normal") matrix(rnorm(n * R), n, R)
       else matrix((rchisq(n * R, df = 3) - 3) / sqrt(6), n, R) # skewed, var 1
  sweep(Z %*% chol(Sigma), 2L, mu, "+")
}

# Judge-skill categorical model: rater j reports the true category with skill
# probability s_j, else guesses from a rater-specific distribution centered at
# a different category (so the rater means differ and kappa_C != kappa_F).
gen_categorical <- function(n, R) {
  C <- 5L; vals <- seq_len(C)
  s <- seq(0.65, 0.78, length.out = R)
  centers <- seq(2, C - 1, length.out = R)
  guess <- vapply(centers, function(ctr) {
    p <- exp(-((vals - ctr)^2) / (2 * 1.2^2)); p / sum(p)
  }, numeric(C))
  Tcat <- sample(vals, n, replace = TRUE)
  X <- matrix(0, n, R)
  for (j in seq_len(R)) {
    skilled <- runif(n) < s[j]
    g <- sample(vals, n, replace = TRUE, prob = guess[, j])
    X[, j] <- ifelse(skilled, Tcat, g)
  }
  X
}

# Two-component Gaussian mixture in which a latent class sets BOTH a mean shift
# and the within-rater correlation, so the observed first rater reveals the class
# and hence the conditional covariance of the remaining raters. The saturated
# normal model assumes a single homoscedastic covariance, so this is the
# continuous stress for FIML: its MAR consistency breaks because the conditional
# second moments depend on the observed anchor. (A pure location mixture, with
# only a common mean shift, does NOT stress FIML -- it stays consistent for
# (mu, Sigma) there despite the nonlinear conditional mean.)
gen_mixture <- function(n, R, rho0 = 0.1, rho1 = 0.85, shift = 1.0) {
  mu <- make_mu(R); d <- sqrt(seq(0.8, 1.4, length.out = R))
  z <- rbinom(n, 1L, 0.5)
  X <- matrix(0, n, R)
  for (cl in 0:1) {
    idx <- which(z == cl)
    if (!length(idx)) next
    rho <- if (cl == 0L) rho0 else rho1
    sh <- if (cl == 0L) -shift else shift
    S <- matrix(rho, R, R); diag(S) <- 1; S <- outer(d, d) * S
    X[idx, ] <- sweep(matrix(rnorm(length(idx) * R), length(idx), R) %*% chol(S),
                      2L, mu + sh, "+")
  }
  X
}

gen_cell <- function(n, R, dist) {
  if (dist == "categorical") gen_categorical(n, R)
  else if (dist == "mixture") gen_mixture(n, R)
  else gen_continuous(n, R, dist)
}

amputate <- function(X, mech, rate) {
  n <- nrow(X); R <- ncol(X)
  if (mech == "mcar") {
    X[matrix(runif(n * R) < rate, n, R)] <- NA
  } else { # observed-anchor MAR: raters 2..R drop with prob rising in rater 1
    a <- log(rate / (1 - rate))
    pr <- plogis(a + 1.1 * scale(X[, 1])[, 1])
    for (j in 2:R) X[runif(n) < pr, j] <- NA
  }
  X
}

# Population Conger / Fleiss from a mean vector and covariance.
truth_kappa <- function(mu, S) {
  R <- length(mu)
  t1 <- sum(S); t2 <- sum(diag(S)); t3 <- sum((mu - mean(mu))^2)
  c(Conger = (t1 - t2) / ((R - 1) * t2 + R * t3),
    Fleiss = (t1 - t2 - t3) / ((R - 1) * (t2 + t3)))
}
coefs <- c("Conger", "Fleiss")

# One replication: long data frame over method x coef.
one_rep <- function(dist, mech, n, R, rate) {
  X <- amputate(gen_cell(n, R, dist), mech, rate)
  X <- X[rowSums(is.finite(X)) > 0, , drop = FALSE]
  colnames(X) <- paste0("r", seq_len(ncol(X)))

  fits <- list(
    listwise = local({
      Xc <- X[stats::complete.cases(X), , drop = FALSE]
      if (nrow(Xc) >= 3L) {
        tryCatch(kappa_continuous(Xc, method = "available", weight = "quadratic"),
                 error = function(e) NULL)
      } else NULL
    }),
    pairwise = tryCatch(kappa_continuous(X, method = "available",
                                         weight = "quadratic"),
                        error = function(e) NULL),
    nt_fiml = tryCatch(kappa_continuous(X, method = "fiml",
                                        weight = "quadratic",
                                        se_type = "sandwich"),
                       error = function(e) NULL)
  )
  do.call(rbind, lapply(names(fits), function(m) {
    f <- fits[[m]]
    est <- se <- c(Conger = NA_real_, Fleiss = NA_real_)
    if (!is.null(f)) { est <- coef(f)[coefs]; se <- sqrt(diag(vcov(f)))[coefs] }
    data.frame(method = m, coef = coefs, estimate = unname(est[coefs]),
               se = unname(se[coefs]), stringsAsFactors = FALSE)
  }))
}

# ---- grid -------------------------------------------------------------------
ref_n <- opts$n; ref_R <- 4L
grid <- rbind(
  # reference: every DGP x mechanism at the reference n and R
  expand.grid(dist = c("normal", "skewed", "categorical", "mixture"),
              mech = c("mcar", "mar"), n = ref_n, R = ref_R,
              stringsAsFactors = FALSE),
  # sample-size sensitivity: normal DGP, smaller n, reference R
  expand.grid(dist = "normal", mech = c("mcar", "mar"),
              n = c(50L, 200L), R = ref_R, stringsAsFactors = FALSE),
  # rater-count sensitivity: normal DGP, two raters, reference n
  expand.grid(dist = "normal", mech = c("mcar", "mar"),
              n = ref_n, R = 2L, stringsAsFactors = FALSE)
)

# Population truth per (dist, R), from a large complete pilot sample.
truth_key <- unique(grid[, c("dist", "R")])
truth_tab <- do.call(rbind, lapply(seq_len(nrow(truth_key)), function(k) {
  set.seed(opts$seed_base + 99000L + k)
  Xp <- gen_cell(2e5L, truth_key$R[k], truth_key$dist[k])
  np <- nrow(Xp)
  tk <- truth_kappa(colMeans(Xp), stats::cov(Xp) * (np - 1) / np)
  data.frame(dist = truth_key$dist[k], R = truth_key$R[k],
             coef = coefs, truth = unname(tk[coefs]), stringsAsFactors = FALSE)
}))

# ---- drive ------------------------------------------------------------------
z975 <- qnorm(0.975)
cell_summaries <- vector("list", nrow(grid))
for (g in seq_len(nrow(grid))) {
  dist <- grid$dist[g]; mech <- grid$mech[g]; n <- grid$n[g]; R <- grid$R[g]
  reps <- vector("list", opts$reps)
  for (b in seq_len(opts$reps)) {
    set.seed(opts$seed_base + 1000L * g + b)
    reps[[b]] <- one_rep(dist, mech, n, R, opts$rate)
  }
  rep_df <- do.call(rbind, reps)
  tr <- truth_tab[truth_tab$dist == dist & truth_tab$R == R, ]
  rep_df$truth <- tr$truth[match(rep_df$coef, tr$coef)]

  key <- interaction(rep_df$method, rep_df$coef, drop = TRUE)
  parts <- lapply(split(rep_df, key), function(d) {
    ok <- is.finite(d$estimate)
    se_ok <- is.finite(d$se) & d$se >= 0
    err <- d$estimate - d$truth
    cover <- ok & se_ok
    data.frame(
      dist = dist, mech = mech, n = n, R = R,
      method = d$method[[1L]], coef = d$coef[[1L]],
      reps = opts$reps, truth = d$truth[[1L]],
      bias = mean(err[ok]),
      rmse = sqrt(mean(err[ok]^2)),
      mean_se = if (any(se_ok)) mean(d$se[se_ok]) else NA_real_,
      coverage95 = if (any(cover)) mean(abs(err[cover]) <= z975 * d$se[cover]) else NA_real_,
      mean_ci_length = if (any(se_ok)) mean(2 * z975 * d$se[se_ok]) else NA_real_,
      valid_reps = sum(ok), stringsAsFactors = FALSE
    )
  })
  cell_summaries[[g]] <- do.call(rbind, parts)
  cat(sprintf("  done %-11s %-4s n=%-3d R=%d\n", dist, mech, n, R))
}
summary_out <- do.call(rbind, cell_summaries)

metadata_out <- data.frame(
  key = c("reps", "n", "rate", "rho", "seed_base", "R_version", "misskappa"),
  value = c(opts$reps, opts$n, opts$rate, 0.6, opts$seed_base,
            R.version.string, as.character(utils::packageVersion("misskappa"))),
  stringsAsFactors = FALSE
)

dir.create(opts$out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(summary_out, file.path(opts$out_dir, "summary.csv"), row.names = FALSE)
write.csv(truth_tab, file.path(opts$out_dir, "truth.csv"), row.names = FALSE)
write.csv(metadata_out, file.path(opts$out_dir, "metadata.csv"), row.names = FALSE)

cat("\n== Conger coverage at reference (n=", ref_n, ", R=", ref_R, ") ==\n", sep = "")
ref <- summary_out[summary_out$coef == "Conger" & summary_out$n == ref_n &
                   summary_out$R == ref_R, c("dist", "mech", "method", "bias",
                                             "coverage95")]
print(format(ref, digits = 3), row.names = FALSE)
cat("\nWrote:\n  ", file.path(opts$out_dir, "summary.csv"),
    "\n  ", file.path(opts$out_dir, "truth.csv"),
    "\n  ", file.path(opts$out_dir, "metadata.csv"), "\n", sep = "")
