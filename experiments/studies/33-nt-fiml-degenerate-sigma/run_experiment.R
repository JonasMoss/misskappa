#!/usr/bin/env Rscript

usage <- function() {
  cat(
"Usage: Rscript run_experiment.R [--smoke|--screen|--full] [options]

Study 33: nt-FIML degenerate Sigma_oo rescue timing and calibration.

Modes:
  --smoke              Tiny timing smoke. Defaults: reps=3, n=10,20, miss=0.2.
  --screen             Medium local screen. Defaults: reps=30, n=10,20,40, miss=0,0.2.
  --full               Original study scale. Defaults: reps=1200, n=10,20,40, miss=0,0.2.

Options:
  --reps N             Replicates per n/miss cell.
  --ns LIST            Comma-separated sample sizes, e.g. 10,20,40.
  --misses LIST        Comma-separated missing fractions, e.g. 0,0.2.
  --truth-n N          Population draw size used for the truth approximation.
  --out-dir DIR        Output directory. Defaults to results/<mode>.
  --seed-base N        Base RNG seed. Default: 20260611.
  --soo-rcond X        Rescue rcond passed as em_options$soo_rcond. Default: 1e-10.
  --deg-rcond X        rcond(Sigma) threshold for degenerate stratum. Default: 1e-8.
  --load-all           Load local r-package with devtools::load_all().
  --no-progress        Suppress per-cell progress messages.
  -h, --help           Show this help.
",
    sep = ""
  )
}

this_file <- function() {
  args <- commandArgs(FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0L) return(normalizePath("run_experiment.R", mustWork = FALSE))
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
}

parse_csv_numeric <- function(x) {
  as.numeric(strsplit(x, ",", fixed = TRUE)[[1L]])
}

parse_args <- function(args) {
  opt <- list(
    mode = "smoke",
    reps = NA_integer_,
    ns = NULL,
    misses = NULL,
    truth_n = NA_integer_,
    out_dir = NULL,
    seed_base = 20260611L,
    soo_rcond = 1e-10,
    deg_rcond = 1e-8,
    load_all = FALSE,
    progress = TRUE
  )

  i <- 1L
  while (i <= length(args)) {
    a <- args[[i]]
    if (a %in% c("-h", "--help")) {
      usage()
      quit(status = 0L)
    } else if (a %in% c("--smoke", "--screen", "--full")) {
      opt$mode <- substring(a, 3L)
    } else if (a == "--load-all") {
      opt$load_all <- TRUE
    } else if (a == "--no-progress") {
      opt$progress <- FALSE
    } else if (grepl("^--[^=]+=", a)) {
      parts <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1L]]
      key <- parts[[1L]]
      val <- paste(parts[-1L], collapse = "=")
      opt <- set_arg(opt, key, val)
    } else if (startsWith(a, "--")) {
      key <- sub("^--", "", a)
      if (i == length(args)) stop("missing value for --", key, call. = FALSE)
      i <- i + 1L
      opt <- set_arg(opt, key, args[[i]])
    } else {
      stop("unknown argument: ", a, call. = FALSE)
    }
    i <- i + 1L
  }

  if (is.na(opt$reps)) {
    opt$reps <- switch(opt$mode, smoke = 3L, screen = 30L, full = 1200L,
                       stop("unknown mode: ", opt$mode, call. = FALSE))
  }
  if (is.null(opt$ns)) {
    opt$ns <- switch(opt$mode, smoke = c(10L, 20L), screen = c(10L, 20L, 40L),
                     full = c(10L, 20L, 40L),
                     stop("unknown mode: ", opt$mode, call. = FALSE))
  }
  if (is.null(opt$misses)) {
    opt$misses <- switch(opt$mode, smoke = 0.2, screen = c(0, 0.2),
                         full = c(0, 0.2),
                         stop("unknown mode: ", opt$mode, call. = FALSE))
  }
  if (is.na(opt$truth_n)) {
    opt$truth_n <- switch(opt$mode, smoke = 50000L, screen = 100000L,
                          full = 400000L,
                          stop("unknown mode: ", opt$mode, call. = FALSE))
  }
  opt
}

set_arg <- function(opt, key, val) {
  key <- gsub("-", "_", key, fixed = TRUE)
  if (key == "reps") opt$reps <- as.integer(val)
  else if (key == "ns") opt$ns <- as.integer(parse_csv_numeric(val))
  else if (key == "misses") opt$misses <- parse_csv_numeric(val)
  else if (key == "truth_n") opt$truth_n <- as.integer(val)
  else if (key == "out_dir") opt$out_dir <- val
  else if (key == "seed_base") opt$seed_base <- as.integer(val)
  else if (key == "soo_rcond") opt$soo_rcond <- as.numeric(val)
  else if (key == "deg_rcond") opt$deg_rcond <- as.numeric(val)
  else stop("unknown option: --", gsub("_", "-", key), call. = FALSE)
  opt
}

safe_system <- function(...) {
  out <- tryCatch(system2(..., stdout = TRUE, stderr = FALSE), error = function(e) NA_character_)
  if (length(out) == 0L) NA_character_ else out[[1L]]
}

write_rect <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE, na = "")
}

fit_timed <- function(X, em_options) {
  start <- proc.time()[["elapsed"]]
  err <- NULL
  fit <- tryCatch(
    suppressWarnings(misskappa::kappa(
      X, estimator = "nt_fiml", weight = "quadratic", em_options = em_options
    )),
    error = function(e) {
      err <<- conditionMessage(e)
      NULL
    }
  )
  list(
    fit = fit,
    ok = !is.null(fit),
    error = if (is.null(err)) NA_character_ else err,
    elapsed_s = proc.time()[["elapsed"]] - start
  )
}

value_named <- function(x, name) {
  if (is.null(x) || length(x) == 0L) return(NA_real_)
  nms <- names(x)
  if (!is.null(nms) && name %in% nms) return(as.numeric(x[[name]]))
  NA_real_
}

rcond_sym <- function(S) {
  if (is.null(S) || any(!is.finite(S))) return(NA_real_)
  vals <- tryCatch(eigen(S, symmetric = TRUE, only.values = TRUE)$values,
                   error = function(e) NA_real_)
  if (any(!is.finite(vals))) return(NA_real_)
  scale <- max(abs(vals))
  if (scale == 0) return(0)
  min(vals) / scale
}

make_calibration <- function(raw, truth, ns, misses, deg_rcond) {
  z <- qnorm(0.975)
  rows <- list()
  for (coef in c("Conger", "Fleiss")) {
    for (n in ns) {
      for (miss in misses) {
        for (st in c("clean", "degenerate")) {
          g <- raw[raw$coef == coef & raw$on_ok & raw$n == n & raw$miss == miss, ]
          if (nrow(g) == 0L) next
          stratum <- ifelse(!is.finite(g$rc) | g$rc < deg_rcond, "degenerate", "clean")
          g <- g[stratum == st, ]
          if (nrow(g) == 0L) next
          sd_est <- if (nrow(g) > 1L) stats::sd(g$est, na.rm = TRUE) else NA_real_
          mean_se <- mean(g$se, na.rm = TRUE)
          rows[[length(rows) + 1L]] <- data.frame(
            coef = coef,
            n = n,
            miss = miss,
            stratum = st,
            nrep = nrow(g),
            bias = mean(g$est, na.rm = TRUE) - truth[[coef]],
            mc_sd = sd_est,
            mean_se = mean_se,
            se_sd = if (is.finite(sd_est) && sd_est > 0) mean_se / sd_est else NA_real_,
            cov95 = mean(abs(g$est - truth[[coef]]) <= z * g$se, na.rm = TRUE),
            med_null_frac = stats::median(g$nf, na.rm = TRUE)
          )
        }
      }
    }
  }
  if (length(rows) == 0L) return(data.frame())
  do.call(rbind, rows)
}

make_recovery <- function(raw) {
  cell <- raw[raw$coef == "Conger", ]
  if (nrow(cell) == 0L) return(data.frame())
  cell$off_fail <- !cell$off_ok
  cell$on_fail <- !cell$on_ok
  stats::aggregate(
    cbind(off_fail, on_fail, off_elapsed_s, on_elapsed_s) ~ n + miss,
    data = cell,
    FUN = mean,
    na.rm = TRUE
  )
}

make_timing <- function(raw) {
  cell <- raw[raw$coef == "Conger", ]
  if (nrow(cell) == 0L) return(data.frame())
  stats::aggregate(
    cbind(off_elapsed_s, on_elapsed_s) ~ n + miss,
    data = cell,
    FUN = function(x) stats::median(x, na.rm = TRUE)
  )
}

script <- this_file()
exp_dir <- dirname(script)
repo_root <- normalizePath(file.path(exp_dir, "..", "..", ".."), mustWork = TRUE)
opt <- parse_args(commandArgs(TRUE))
if (is.null(opt$out_dir)) opt$out_dir <- file.path(exp_dir, "results", opt$mode)
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

if (opt$load_all) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("--load-all requires the devtools package.", call. = FALSE)
  }
  suppressMessages(devtools::load_all(file.path(repo_root, "r-package"), quiet = TRUE))
} else {
  suppressPackageStartupMessages(library(misskappa))
}

R <- 5L
C <- 5L
rho <- 0.6
bias <- seq(-0.5, 0.5, length.out = R)

gen <- function(n) {
  a <- stats::rnorm(n)
  Z <- sqrt(rho) * matrix(a, n, R) + sqrt(1 - rho) * matrix(stats::rnorm(n * R), n, R)
  Z <- sweep(Z, 2L, bias, "+")
  br <- c(-Inf, stats::qnorm(seq_len(C - 1L) / C), Inf)
  matrix(as.integer(cut(Z, breaks = br, labels = FALSE)), n, R)
}

kfun <- function(mu, S) {
  t1 <- sum(S)
  t2 <- sum(diag(S))
  t3 <- sum((mu - mean(mu))^2)
  c(
    Conger = (t1 - t2) / ((R - 1) * t2 + R * t3),
    Fleiss = (t1 - t2 - t3) / ((R - 1) * (t2 + t3))
  )
}

set.seed(opt$seed_base)
Xb <- gen(opt$truth_n)
mub <- colMeans(Xb)
Sb <- stats::cov(Xb) * (nrow(Xb) - 1) / nrow(Xb)
truth <- kfun(mub, Sb)

cat(sprintf("truth: Conger=%.4f Fleiss=%.4f\n", truth[["Conger"]], truth[["Fleiss"]]))
cat(sprintf(
  "mode=%s reps=%d ns=%s misses=%s truth_n=%d out_dir=%s\n",
  opt$mode, opt$reps, paste(opt$ns, collapse = ","),
  paste(opt$misses, collapse = ","), opt$truth_n, opt$out_dir
))

rows <- list()
for (n in opt$ns) {
  for (miss in opt$misses) {
    if (opt$progress) {
      cat(sprintf("cell n=%d miss=%.3g reps=%d\n", n, miss, opt$reps))
      flush.console()
    }
    for (rep_id in seq_len(opt$reps)) {
      set.seed(opt$seed_base + n * 100000L + as.integer(round(miss * 1000)) * 1000L + rep_id)
      X <- gen(n)
      if (miss > 0) X[matrix(stats::runif(n * R) < miss, n, R)] <- NA

      off <- fit_timed(X, list())
      on <- fit_timed(X, list(soo_rcond = opt$soo_rcond))
      rc <- if (on$ok) rcond_sym(on$fit$moments$Sigma) else NA_real_
      se <- if (on$ok) sqrt(diag(on$fit$vcov)) else c(Conger = NA_real_, Fleiss = NA_real_)
      if (is.null(names(se)) && length(se) == 2L) names(se) <- c("Conger", "Fleiss")
      if (!is.null(names(se))) {
        se <- se[c("Conger", "Fleiss")]
      }

      for (coef in c("Conger", "Fleiss")) {
        rows[[length(rows) + 1L]] <- data.frame(
          rep_id = rep_id,
          n = n,
          miss = miss,
          coef = coef,
          off_ok = off$ok,
          on_ok = on$ok,
          off_error = off$error,
          on_error = on$error,
          off_elapsed_s = off$elapsed_s,
          on_elapsed_s = on$elapsed_s,
          est = if (on$ok) value_named(on$fit$estimates, coef) else NA_real_,
          se = if (coef %in% names(se)) as.numeric(se[[coef]]) else NA_real_,
          nf = if (on$ok) value_named(on$fit$null_frac, coef) else NA_real_,
          rc = rc,
          truth = truth[[coef]]
        )
      }
    }
  }
}

raw <- do.call(rbind, rows)
recovery <- make_recovery(raw)
timing <- make_timing(raw)
calibration <- make_calibration(raw, truth, opt$ns, opt$misses, opt$deg_rcond)
truth_df <- data.frame(coef = names(truth), truth = as.numeric(truth))
metadata <- data.frame(
  key = c(
    "timestamp", "mode", "reps", "ns", "misses", "truth_n", "seed_base",
    "soo_rcond", "deg_rcond", "load_all", "git_branch", "git_sha",
    "r_version", "misskappa_version"
  ),
  value = c(
    format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    opt$mode,
    as.character(opt$reps),
    paste(opt$ns, collapse = ","),
    paste(opt$misses, collapse = ","),
    as.character(opt$truth_n),
    as.character(opt$seed_base),
    format(opt$soo_rcond, scientific = TRUE),
    format(opt$deg_rcond, scientific = TRUE),
    as.character(opt$load_all),
    safe_system("git", c("-C", repo_root, "branch", "--show-current")),
    safe_system("git", c("-C", repo_root, "rev-parse", "--short", "HEAD")),
    R.version.string,
    as.character(utils::packageVersion("misskappa"))
  )
)

paths <- c(
  raw = file.path(opt$out_dir, "raw.csv"),
  recovery = file.path(opt$out_dir, "recovery.csv"),
  timing = file.path(opt$out_dir, "timing.csv"),
  calibration = file.path(opt$out_dir, "calibration.csv"),
  truth = file.path(opt$out_dir, "truth.csv"),
  metadata = file.path(opt$out_dir, "metadata.csv")
)
write_rect(raw, paths[["raw"]])
write_rect(recovery, paths[["recovery"]])
write_rect(timing, paths[["timing"]])
write_rect(calibration, paths[["calibration"]])
write_rect(truth_df, paths[["truth"]])
write_rect(metadata, paths[["metadata"]])

cat("\n=== recovery ===\n")
print(recovery, row.names = FALSE)
cat("\n=== median timing seconds per fit ===\n")
print(timing, row.names = FALSE)
cat("\n=== calibration ===\n")
print(calibration, row.names = FALSE)
cat("\nwrote:\n")
cat(paste0("  ", unname(paths), collapse = "\n"), "\n")
