#!/usr/bin/env Rscript
#
# Experiment 22: validation of the normal-FIML quadratic (Conger / Fleiss)
# kappa (misskappa::kappa_quadratic_fiml, surfaced as
# kappa_continuous(method = "fiml")). One question: does the saturated-EM
# covariance and the delta-method standard error of the quadratic kappas agree
# with the established references and the bootstrap truth?
#
# This is the kappa counterpart of experiment 17 (alpha). The saturated-EM
# moments are checked against magmaan's estimate_saturated_em_moments(); the
# point and SE are checked against the pairwise-available estimator
# kappa_continuous(method = "available") and a nonparametric case bootstrap (the
# ground truth). The DGP carries rater-specific means so t3 = sum (mu_j - mubar)^2
# is non-trivial -- this exercises the mean block of the kappa gradient, the one
# piece that has no analogue in alpha. The nonnormal cell is the discriminator:
# there the normal-theory delta SE should undershoot while the sandwich SE tracks
# the bootstrap.

suppressPackageStartupMessages({
  library(misskappa)
})

usage <- function(status = 0L) {
  cat(paste0(
    "Usage: Rscript run_experiment.R [options]\n\n",
    "Options:\n",
    "  --help, -h         Show this help and exit.\n",
    "  --smoke            Cheap check: n=400, boot=100.\n",
    "  --n N              Sample size per cell. Default: 800.\n",
    "  --R R              Number of raters. Default: 4.\n",
    "  --boot N           Bootstrap replicates for the truth SE. Default: 600.\n",
    "  --seed-base N      Base seed. Default: 22000.\n",
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

opts <- list(n = 800L, R = 4L, boot = 600L, seed_base = 22000L,
             out_dir = file.path(script_dir, "results"))
args <- commandArgs(TRUE)
i <- 1L
while (i <= length(args)) {
  a <- args[[i]]
  if (a %in% c("--help", "-h")) usage(0L)
  else if (a == "--smoke") { opts$n <- 400L; opts$boot <- 100L }
  else if (a == "--n") { opts$n <- as.integer(args[[i + 1L]]); i <- i + 1L }
  else if (a == "--R") { opts$R <- as.integer(args[[i + 1L]]); i <- i + 1L }
  else if (a == "--boot") { opts$boot <- as.integer(args[[i + 1L]]); i <- i + 1L }
  else if (a == "--seed-base") { opts$seed_base <- as.integer(args[[i + 1L]]); i <- i + 1L }
  else if (a == "--out-dir") { opts$out_dir <- args[[i + 1L]]; i <- i + 1L }
  else { cat("Unknown argument:", a, "\n"); usage(1L) }
  i <- i + 1L
}

has <- function(pkg) requireNamespace(pkg, quietly = TRUE)

# ---- DGP --------------------------------------------------------------------
# Positively correlated raters (agreement) with heterogeneous variances and
# rater-specific means (so the Conger and Fleiss kappas differ and t3 > 0).
make_Sigma <- function(R, rho = 0.6) {
  S <- matrix(rho, R, R); diag(S) <- 1
  d <- sqrt(seq(0.8, 1.4, length.out = R))
  outer(d, d) * S
}
make_mu <- function(R) seq(-0.45, 0.45, length.out = R)

gen <- function(n, Sigma, mu, dist) {
  R <- ncol(Sigma)
  Z <- if (dist == "normal") matrix(rnorm(n * R), n, R)
       else matrix((rchisq(n * R, df = 3) - 3) / sqrt(6), n, R) # skewed, var 1
  sweep(Z %*% chol(Sigma), 2L, mu, "+")
}
amputate <- function(X, mech, rate = 0.18) {
  n <- nrow(X); R <- ncol(X)
  if (mech == "complete") return(X)
  if (mech == "mcar") {
    X[matrix(runif(n * R) < rate, n, R)] <- NA
  } else { # MAR: raters 2..R drop with prob increasing in observed rater 1
    pr <- plogis(scale(X[, 1]) * 0.9 - 0.9)
    for (j in 2:R) X[runif(n) < pr, j] <- NA
  }
  X
}

run_cell <- function(dist, mech, Sigma, mu, n, boot, seed) {
  set.seed(seed)
  X <- amputate(gen(n, Sigma, mu, dist), mech)
  X <- X[rowSums(is.finite(X)) > 0, , drop = FALSE]
  colnames(X) <- paste0("r", seq_len(ncol(X)))
  nb <- nrow(X)

  fit_s <- kappa_continuous(X, method = "fiml", weight = "quadratic", se_type = "sandwich")
  fit_n <- kappa_continuous(X, method = "fiml", weight = "quadratic", se_type = "normal")
  Sig <- fit_s$moments$Sigma

  # Pairwise-available reference (point + SE), continuous quadratic kernel.
  ref <- tryCatch(kappa_continuous(X, method = "available", weight = "quadratic"),
                  error = function(e) NULL)
  ref_est <- if (is.null(ref)) c(Conger = NA, Fleiss = NA) else coef(ref)
  ref_se  <- if (is.null(ref)) c(Conger = NA, Fleiss = NA) else sqrt(diag(vcov(ref)))

  # magmaan oracle: covariance, and the saturated-moment sandwich SE through
  # acov = H^-1 J H^-1, an independent check of the influence-function vcov.
  d_mag <- NA_real_
  se_mag <- c(Conger = NA_real_, Fleiss = NA_real_)
  if (has("magmaan")) {
    mask <- is.finite(X); storage.mode(mask) <- "logical"
    mag <- magmaan::magmaan_core$estimate_saturated_em_moments(
      list(X = X, mask = mask))
    d_mag <- max(abs(Sig - mag$cov[[1L]]))
    G <- misskappa:::.kqf_grad(fit_s$moments$mu, fit_s$moments$Sigma)$G
    se_mag <- sqrt(diag(G %*% as.matrix(mag$acov) %*% t(G)))
  }

  # Bootstrap truth SE (resample subjects, refit the point).
  bvals <- vapply(seq_len(boot), function(b) {
    Xb <- X[sample.int(nb, nb, replace = TRUE), , drop = FALSE]
    fb <- tryCatch(kappa_continuous(Xb, method = "fiml", weight = "quadratic", se_type = "normal"),
                   error = function(e) NULL)
    if (is.null(fb)) c(NA_real_, NA_real_) else coef(fb)[c("Conger", "Fleiss")]
  }, numeric(2))
  se_boot <- apply(bvals, 1L, sd, na.rm = TRUE)

  est <- coef(fit_s)[c("Conger", "Fleiss")]
  se_s <- sqrt(diag(vcov(fit_s)))[c("Conger", "Fleiss")]
  se_n <- sqrt(diag(vcov(fit_n)))[c("Conger", "Fleiss")]

  data.frame(
    dist = dist, mech = mech, n = nb,
    coef = c("Conger", "Fleiss"),
    est = unname(est),
    est_ref = unname(ref_est[c("Conger", "Fleiss")]),
    se_sandwich = unname(se_s),
    se_normal = unname(se_n),
    se_available = unname(ref_se[c("Conger", "Fleiss")]),
    se_magmaan = unname(se_mag[c("Conger", "Fleiss")]),
    se_boot = unname(se_boot),
    max_cov_err_magmaan = d_mag,
    em_iter = fit_s$moments$iterations,
    stringsAsFactors = FALSE
  )
}

# ---- drive ------------------------------------------------------------------
Sigma <- make_Sigma(opts$R)
mu <- make_mu(opts$R)
grid <- rbind(
  data.frame(dist = "normal",    mech = "complete"),
  data.frame(dist = "normal",    mech = "mcar"),
  data.frame(dist = "normal",    mech = "mar"),
  data.frame(dist = "nonnormal", mech = "mcar")
)
rows <- lapply(seq_len(nrow(grid)), function(k) {
  run_cell(grid$dist[k], grid$mech[k], Sigma, mu, opts$n, opts$boot,
           opts$seed_base + k)
})
validation <- do.call(rbind, rows)

metadata <- data.frame(
  key = c("n", "R", "boot", "seed_base", "rho", "R_version", "misskappa",
          "magmaan"),
  value = c(opts$n, opts$R, opts$boot, opts$seed_base, 0.6,
            R.version.string,
            as.character(utils::packageVersion("misskappa")),
            if (has("magmaan")) as.character(utils::packageVersion("magmaan")) else "NA"),
  stringsAsFactors = FALSE
)

dir.create(opts$out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(validation, file.path(opts$out_dir, "validation.csv"), row.names = FALSE)
write.csv(metadata, file.path(opts$out_dir, "metadata.csv"), row.names = FALSE)

cat("\n== quadratic kappa SE: sandwich / normal / available / magmaan / boot ==\n")
print(format(validation[, c("dist", "mech", "coef", "est", "se_sandwich",
                            "se_normal", "se_available", "se_magmaan",
                            "se_boot")],
             digits = 4), row.names = FALSE)
cat("\n== max covariance error vs magmaan oracle ==\n")
print(format(validation[validation$coef == "Conger",
                        c("dist", "mech", "max_cov_err_magmaan")], digits = 3),
      row.names = FALSE)
cat("\nWrote:\n  ", file.path(opts$out_dir, "validation.csv"),
    "\n  ", file.path(opts$out_dir, "metadata.csv"), "\n", sep = "")
