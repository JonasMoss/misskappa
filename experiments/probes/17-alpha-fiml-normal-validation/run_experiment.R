#!/usr/bin/env Rscript
#
# Experiment 17: validation of the normal-FIML continuous-item coefficient
# alpha (misskappa::alpha_continuous). One question: does the saturated-EM
# covariance and its delta-method standard error agree with the established
# implementations? Moments are checked against lavaan's saturated h1 estimator
# and magmaan's estimate_saturated_em_moments(); the alpha point and the
# sandwich SE are checked against coefficientalpha (varphi = 0); and the SE is
# checked against a nonparametric case bootstrap (the ground truth). The
# nonnormal cell is the discriminator: there the normal-theory delta SE should
# undershoot while the sandwich SE tracks the bootstrap.

suppressPackageStartupMessages({
  library(misskappa)
})

usage <- function(status = 0L) {
  cat(paste0(
    "Usage: Rscript run_experiment.R [options]\n\n",
    "Options:\n",
    "  --help, -h         Show this help and exit.\n",
    "  --smoke            Cheap check: n=300, boot=100.\n",
    "  --n N              Sample size per cell. Default: 600.\n",
    "  --p P              Number of items. Default: 6.\n",
    "  --boot N           Bootstrap replicates for the truth SE. Default: 600.\n",
    "  --seed-base N      Base seed. Default: 17000.\n",
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

opts <- list(n = 600L, p = 6L, boot = 600L, seed_base = 17000L,
             out_dir = file.path(script_dir, "results"))
args <- commandArgs(TRUE)
i <- 1L
while (i <= length(args)) {
  a <- args[[i]]
  if (a %in% c("--help", "-h")) usage(0L)
  else if (a == "--smoke") { opts$n <- 300L; opts$boot <- 100L }
  else if (a == "--n") { opts$n <- as.integer(args[[i + 1L]]); i <- i + 1L }
  else if (a == "--p") { opts$p <- as.integer(args[[i + 1L]]); i <- i + 1L }
  else if (a == "--boot") { opts$boot <- as.integer(args[[i + 1L]]); i <- i + 1L }
  else if (a == "--seed-base") { opts$seed_base <- as.integer(args[[i + 1L]]); i <- i + 1L }
  else if (a == "--out-dir") { opts$out_dir <- args[[i + 1L]]; i <- i + 1L }
  else { cat("Unknown argument:", a, "\n"); usage(1L) }
  i <- i + 1L
}

has <- function(pkg) requireNamespace(pkg, quietly = TRUE)

# ---- DGPs -------------------------------------------------------------------
make_Sigma <- function(p, rho = 0.4) {
  S <- matrix(rho, p, p); diag(S) <- 1
  d <- sqrt(seq(0.7, 1.6, length.out = p))
  outer(d, d) * S
}
gen <- function(n, Sigma, dist) {
  p <- ncol(Sigma)
  Z <- if (dist == "normal") matrix(rnorm(n * p), n, p)
       else matrix((rchisq(n * p, df = 3) - 3) / sqrt(6), n, p) # skewed, var 1
  Z %*% chol(Sigma)
}
amputate <- function(X, mech, rate = 0.18) {
  n <- nrow(X); p <- ncol(X)
  if (mech == "complete") return(X)
  if (mech == "mcar") {
    X[matrix(runif(n * p) < rate, n, p)] <- NA
  } else { # MAR: items 2..p drop with prob increasing in observed item 1
    pr <- plogis(scale(X[, 1]) * 0.9 - 0.9)
    for (j in 2:p) X[runif(n) < pr, j] <- NA
  }
  X
}

alpha_from_S <- function(S) {
  p <- ncol(S); (p / (p - 1)) * (1 - sum(diag(S)) / sum(S))
}

run_cell <- function(dist, mech, Sigma, n, boot, seed) {
  set.seed(seed)
  X <- amputate(gen(n, Sigma, dist), mech)
  X <- X[rowSums(is.finite(X)) > 0, , drop = FALSE]
  colnames(X) <- paste0("y", seq_len(ncol(X)))

  fit_s <- alpha_continuous(X, se_type = "sandwich")
  fit_n <- alpha_continuous(X, se_type = "normal")
  Sig <- fit_s$moments$Sigma

  d_lav <- NA_real_
  if (has("lavaan")) {
    Mp <- lavaan:::lav_data_mi_patterns(X)
    lav <- lavaan:::lav_mvn_mi_h1_est_moments(
      y = X, mp = Mp, tol = 1e-12, max_iter = 10000L)
    d_lav <- max(abs(Sig - lav$Sigma))
  }
  d_mag <- NA_real_
  if (has("magmaan")) {
    mask <- is.finite(X); storage.mode(mask) <- "logical"
    mag <- magmaan::magmaan_core$estimate_saturated_em_moments(
      list(X = X, mask = mask))
    d_mag <- max(abs(Sig - mag$cov[[1L]]))
  }
  ca_alpha <- ca_se <- NA_real_
  if (has("coefficientalpha")) {
    ca <- tryCatch({
      invisible(utils::capture.output(
        o <- coefficientalpha::alpha(X, varphi = 0, se = TRUE, test = FALSE,
                                     silent = TRUE)))
      o
    }, error = function(e) NULL)
    if (!is.null(ca)) { ca_alpha <- as.numeric(ca$alpha); ca_se <- as.numeric(ca$se) }
  }

  # bootstrap truth SE
  nb <- nrow(X)
  vals <- vapply(seq_len(boot), function(b) {
    Xb <- X[sample.int(nb, nb, replace = TRUE), , drop = FALSE]
    fb <- tryCatch(alpha_continuous(Xb, se_type = "normal"),
                   error = function(e) NULL)
    if (is.null(fb)) NA_real_ else unname(coef(fb)["alpha"])
  }, numeric(1))

  data.frame(
    dist = dist, mech = mech, n = nb,
    alpha = unname(coef(fit_s)["alpha"]), alpha_ca = ca_alpha,
    se_sandwich = sqrt(vcov(fit_s)[1, 1]),
    se_normal = sqrt(vcov(fit_n)[1, 1]),
    se_ca = ca_se,
    se_boot = sd(vals, na.rm = TRUE),
    max_cov_err_lavaan = d_lav,
    max_cov_err_magmaan = d_mag,
    em_iter = fit_s$moments$iterations,
    stringsAsFactors = FALSE
  )
}

# ---- drive ------------------------------------------------------------------
Sigma <- make_Sigma(opts$p)
grid <- rbind(
  data.frame(dist = "normal",    mech = "complete"),
  data.frame(dist = "normal",    mech = "mcar"),
  data.frame(dist = "normal",    mech = "mar"),
  data.frame(dist = "nonnormal", mech = "mcar")
)
rows <- lapply(seq_len(nrow(grid)), function(k) {
  run_cell(grid$dist[k], grid$mech[k], Sigma, opts$n, opts$boot,
           opts$seed_base + k)
})
validation <- do.call(rbind, rows)

metadata <- data.frame(
  key = c("n", "p", "boot", "seed_base", "R", "misskappa", "lavaan",
          "magmaan", "coefficientalpha"),
  value = c(opts$n, opts$p, opts$boot, opts$seed_base,
            R.version.string,
            as.character(utils::packageVersion("misskappa")),
            if (has("lavaan")) as.character(utils::packageVersion("lavaan")) else "NA",
            if (has("magmaan")) as.character(utils::packageVersion("magmaan")) else "NA",
            if (has("coefficientalpha")) as.character(utils::packageVersion("coefficientalpha")) else "NA"),
  stringsAsFactors = FALSE
)

dir.create(opts$out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(validation, file.path(opts$out_dir, "validation.csv"), row.names = FALSE)
write.csv(metadata, file.path(opts$out_dir, "metadata.csv"), row.names = FALSE)

cat("\n== alpha SE: sandwich / normal-theory / coefficientalpha / bootstrap ==\n")
print(format(validation[, c("dist", "mech", "alpha", "se_sandwich",
                            "se_normal", "se_ca", "se_boot")], digits = 4),
      row.names = FALSE)
cat("\n== max covariance error vs oracles ==\n")
print(format(validation[, c("dist", "mech", "max_cov_err_lavaan",
                            "max_cov_err_magmaan")], digits = 3),
      row.names = FALSE)
cat("\nWrote:\n  ", file.path(opts$out_dir, "validation.csv"),
    "\n  ", file.path(opts$out_dir, "metadata.csv"), "\n", sep = "")
