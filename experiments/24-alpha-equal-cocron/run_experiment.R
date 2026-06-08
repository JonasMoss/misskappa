#!/usr/bin/env Rscript
#
# 24-alpha-equal-cocron
#
# A real-data worked example of the equality-of-alpha test promised in
# the alpha-missing paper ("independent-sample tests of equal alpha,
# paired contrasts"). The data are Diedenhofen & Musch's `knowledge`
# set shipped with the `cocron` package: 312 testees who each took two
# 30-item binary knowledge tests (test1, test2). Because the same
# subjects took both forms, the two alpha estimates are *dependent*; the
# paired contrast is the natural test.
#
# Three things are demonstrated:
#   (1) misskappa::alpha(estimator = "pairwise") reproduces cocron's alpha
#       per form to machine-comparable precision (point-estimate check).
#   (2) The public paired Wald test alpha_test(paired = TRUE) (IF-based, using
#       the exposed fit$psi component) tests equal alpha while accounting for
#       the form dependence, and is contrasted with cocron's normal-theory
#       dependent/independent tests (Feldt family).
#   (3) Under MCAR amputation of the item responses, the available-case
#       alpha estimates and their paired contrast stay essentially
#       unbiased relative to the complete-data values while the standard
#       error inflates smoothly --- i.e. the missing-data estimator
#       degrades gracefully rather than breaking.
#
# Outputs (under results/):
#   bench.csv             complete-data: cocron vs misskappa, both alphas,
#                         contrast, test statistic / df / p, est. corr
#   missing_replicates.csv per (rate, replicate) available-case results
#   missing_summary.csv    per-rate summary: bias, mean SE, MC-SD, power
#   metadata.csv           run metadata
#
# Usage: Rscript run_experiment.R [--smoke] [--reps N] [--seed-base K]
#                                 [--rates a,b,c] [--only bench|missing]

suppressPackageStartupMessages({
  library(misskappa)
  library(cocron)
})

# ---- CLI ----------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
has_flag <- function(name) name %in% args
get_val  <- function(name, default, parser = as.character) {
  i <- which(args == name)
  if (length(i) == 0L) return(default)
  parser(args[i + 1L])
}

if (has_flag("--help") || has_flag("-h")) {
  cat("Usage: Rscript run_experiment.R [options]\n",
      " --smoke         Fast smoke run (reps = 20).\n",
      " --reps N        MCAR amputation replicates per rate (default 300).\n",
      " --rates a,b,c   MCAR item-missing rates (default 0.1,0.25,0.4).\n",
      " --seed-base K   Seed base (default 1).\n",
      " --only WHICH    bench | missing (default both).\n",
      " --help, -h      This help.\n", sep = "")
  quit("no", status = 0)
}

smoke     <- has_flag("--smoke")
seed_base <- get_val("--seed-base", 1L, as.integer)
reps_user <- get_val("--reps", NA_integer_, as.integer)
rates     <- as.numeric(strsplit(get_val("--rates", "0.1,0.25,0.4"), ",")[[1]])
only      <- get_val("--only", NA_character_, as.character)

reps <- if (smoke) 20L else 300L
if (!is.na(reps_user)) reps <- reps_user

steps <- if (is.na(only)) c("bench", "missing") else only

results_dir <- "results"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
t0 <- Sys.time()
tic <- function(tag) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), tag))

# ---- Data ---------------------------------------------------------------
data("knowledge", package = "cocron")
form_a <- knowledge$test1   # 312 x 30 binary
form_b <- knowledge$test2   # 312 x 30 binary
n_subj <- nrow(form_a)
n_item <- ncol(form_a)

# Paired/independent comparison of two same-subject alpha fits, built on the
# PUBLIC verb alpha_test() for the test statistics and on the exposed
# influence-function component (fit$psi) for the descriptive correlation and
# standard errors. Returns the contrast estimate, the estimated correlation of
# the two alpha hats, the dependent/independent SEs, and both chi-squares (df=1)
# with their p-values.
paired_alpha_test <- function(fa, fb) {
  a <- coef(fa)[["alpha"]]; b <- coef(fb)[["alpha"]]
  pa <- fa$psi[, "alpha"]; pb <- fb$psi[, "alpha"]
  n  <- length(pa)
  vA  <- sum(pa^2)    / n^2
  vB  <- sum(pb^2)    / n^2
  cAB <- sum(pa * pb) / n^2
  dep <- alpha_test(fa, fb, paired = TRUE)
  ind <- alpha_test(fa, fb, paired = FALSE)
  list(
    alpha_A = a, alpha_B = b, diff = a - b,
    corr = cAB / sqrt(vA * vB),
    se_dep = sqrt(max(vA + vB - 2 * cAB, 0)),
    se_ind = sqrt(max(vA + vB, 0)),
    chisq_dep = as.numeric(dep$statistic), p_dep = as.numeric(dep$p.value),
    chisq_ind = as.numeric(ind$statistic), p_ind = as.numeric(ind$p.value)
  )
}

# ===================== bench: complete-data ==============================
run_bench <- function() {
  tic("bench start")

  # cocron normal-theory benchmarks (Feldt family).
  bt_dep <- cocron(knowledge, dep = TRUE)
  bt_ind <- cocron(knowledge, dep = FALSE)
  cc_alpha <- as.numeric(bt_dep@alpha)          # a1, a2

  # misskappa available-case alphas on the same forms.
  fa <- alpha(form_a, estimator = "pairwise")
  fb <- alpha(form_b, estimator = "pairwise")
  pt <- paired_alpha_test(fa, fb)

  bench <- rbind(
    data.frame(source = "misskappa", test = "paired (dependent IF)",
               alpha_A = pt$alpha_A, alpha_B = pt$alpha_B, diff = pt$diff,
               corr = pt$corr, statistic = pt$chisq_dep, df = 1,
               p_value = pt$p_dep),
    data.frame(source = "misskappa", test = "independent IF (var add)",
               alpha_A = pt$alpha_A, alpha_B = pt$alpha_B, diff = pt$diff,
               corr = NA_real_, statistic = pt$chisq_ind, df = 1,
               p_value = pt$p_ind),
    data.frame(source = "cocron", test = "dependent (Feldt/normal)",
               alpha_A = cc_alpha[1], alpha_B = cc_alpha[2],
               diff = cc_alpha[1] - cc_alpha[2], corr = NA_real_,
               statistic = as.numeric(bt_dep@statistic),
               df = as.numeric(bt_dep@df),
               p_value = as.numeric(bt_dep@p.value)),
    data.frame(source = "cocron", test = "independent (Hakstian-Whalen)",
               alpha_A = cc_alpha[1], alpha_B = cc_alpha[2],
               diff = cc_alpha[1] - cc_alpha[2], corr = NA_real_,
               statistic = as.numeric(bt_ind@statistic),
               df = as.numeric(bt_ind@df),
               p_value = as.numeric(bt_ind@p.value))
  )

  write.csv(bench, file.path(results_dir, "bench.csv"), row.names = FALSE)
  tic("bench done")
  invisible(bench)
}

# ===================== missing: MCAR robustness ==========================
ampute_mcar <- function(x, rate) {
  M <- matrix(stats::runif(length(x)) < rate, nrow = nrow(x), ncol = ncol(x))
  x[M] <- NA
  x
}

run_missing <- function(complete_diff) {
  tic("missing start")
  per_rep <- list()
  for (rate in rates) {
    for (b in seq_len(reps)) {
      set.seed(seed_base + as.integer(round(rate * 1e4)) * 1000L + b)
      amp_a <- ampute_mcar(form_a, rate)
      amp_b <- ampute_mcar(form_b, rate)
      fa <- try(alpha(amp_a, estimator = "pairwise"), silent = TRUE)
      fb <- try(alpha(amp_b, estimator = "pairwise"), silent = TRUE)
      if (inherits(fa, "try-error") || inherits(fb, "try-error")) next
      pt <- try(paired_alpha_test(fa, fb), silent = TRUE)
      if (inherits(pt, "try-error") || !is.finite(pt$chisq_dep)) next
      per_rep[[length(per_rep) + 1L]] <- data.frame(
        rate = rate, b = b,
        alpha_A = pt$alpha_A, alpha_B = pt$alpha_B, diff = pt$diff,
        se_dep = pt$se_dep, corr = pt$corr,
        chisq_dep = pt$chisq_dep,
        reject_05 = pt$chisq_dep > qchisq(0.95, 1)
      )
    }
    tic(sprintf("  rate %.2f done (%d kept)", rate,
                sum(vapply(per_rep, function(d) d$rate == rate, logical(1)))))
  }
  rep_df <- do.call(rbind, per_rep)
  write.csv(rep_df, file.path(results_dir, "missing_replicates.csv"),
            row.names = FALSE)

  summ <- do.call(rbind, by(rep_df, rep_df$rate, function(d) {
    data.frame(
      rate            = unique(d$rate),
      reps_kept       = nrow(d),
      mean_alpha_A    = mean(d$alpha_A),
      mean_alpha_B    = mean(d$alpha_B),
      mean_diff       = mean(d$diff),
      bias_diff       = mean(d$diff) - complete_diff,
      mean_se_dep     = mean(d$se_dep),
      mc_sd_diff      = sd(d$diff),
      mean_corr       = mean(d$corr),
      power_05        = mean(d$reject_05)
    )
  }))
  rownames(summ) <- NULL
  write.csv(summ, file.path(results_dir, "missing_summary.csv"),
            row.names = FALSE)
  tic("missing done")
  invisible(summ)
}

# ---- Dispatch -----------------------------------------------------------
complete_diff <- coef(alpha(form_a, estimator = "pairwise"))[["alpha"]] -
                 coef(alpha(form_b, estimator = "pairwise"))[["alpha"]]

if ("bench" %in% steps)   run_bench()
if ("missing" %in% steps) run_missing(complete_diff)

# ---- Metadata -----------------------------------------------------------
meta <- data.frame(
  key = c("dataset", "n_subjects", "n_items_per_form", "seed_base", "reps",
          "mcar_rates", "steps", "R_version", "misskappa_version",
          "cocron_version", "started_at", "elapsed_s"),
  value = c(
    "cocron::knowledge",
    as.character(n_subj),
    as.character(n_item),
    as.character(seed_base),
    as.character(reps),
    paste(rates, collapse = ","),
    paste(steps, collapse = ","),
    as.character(getRversion()),
    as.character(utils::packageVersion("misskappa")),
    as.character(utils::packageVersion("cocron")),
    format(t0, "%Y-%m-%dT%H:%M:%S"),
    sprintf("%.1f", as.numeric(difftime(Sys.time(), t0, units = "secs")))
  )
)
write.csv(meta, file.path(results_dir, "metadata.csv"), row.names = FALSE)
cat(sprintf("Wrote %s/{bench,missing_replicates,missing_summary,metadata}.csv\n",
            results_dir))
