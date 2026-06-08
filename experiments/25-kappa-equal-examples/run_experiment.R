#!/usr/bin/env Rscript
#
# 25-kappa-equal-examples
#
# Real-data worked examples of the equality tests promised in the
# quadratic-kappa paper:
#
#   "The same delta-method covariance does yield a Wald test on
#    kappa^(1) - kappa^(2) for two independent samples ... A
#    straightforward modification handles joint tests on multiple
#    dependent kappa estimates from a single rating matrix ... testing
#    whether the two-rater quadratic kappa is constant across all rater
#    pairs."
#
# We realise both structures with quadratically weighted Conger's kappa
# on published ordinal rating data, plus an optional external benchmark.
#
#   Case A  Within-data homogeneity. irr::anxiety (3 raters x 20 subjects,
#           6-point ordinal). Fit the two-rater quadratic Conger kappa on
#           each of the 3 rater pairs via the paper's NT-FIML estimator
#           (kappa(estimator = "nt_fiml"), which exposes per-subject
#           influence functions) and test that all pairwise kappas are equal
#           with kappa_test(..., paired = TRUE) (df = 2).
#
#   Case B  Independent two-sample. The Westlund & Kurland multiple
#           sclerosis data as tabulated by Landis & Koch (1977): two
#           neurologists classify patients into 4 ordinal diagnostic
#           classes, separately in Winnipeg (n = 149) and New Orleans
#           (n = 69). Test equal quadratic Conger kappa across the two
#           independent samples (variances add; df = 1).
#
#   Case C  External benchmark (OPTIONAL, needs the 'multiagree' package).
#           Vanbelle (2017) FEES weighted-kappa comparison. Skipped with a
#           message unless multiagree is installed; see the report.
#
# Outputs (under results/):
#   caseA_pairs.csv     per-pair quadratic Conger kappa + SE (anxiety)
#   caseA_test.csv      joint homogeneity Wald test across rater pairs
#   caseB_kappa.csv     per-city Conger kappa + SE, each weight scheme
#   caseB_test.csv      two-sample equal-kappa tests, each weight scheme
#   caseC_status.csv    whether the multiagree benchmark ran
#   metadata.csv        run metadata
#
# Usage: Rscript run_experiment.R [--only A,B,C] [--help]

suppressPackageStartupMessages({
  library(misskappa)
  library(irr)
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
      " --only A,B,C,D Run a subset of cases (default all).\n",
      "                A within-data homogeneity (anxiety);\n",
      "                B independent two-sample (MS data);\n",
      "                C FEES benchmark (needs multiagree);\n",
      "                D dependent two-construct on real missing data (mcduff).\n",
      " --smoke        Accepted for interface parity; the cases are\n",
      "                deterministic single-dataset analyses.\n",
      " --help, -h     This help.\n", sep = "")
  quit("no", status = 0)
}
cases <- toupper(strsplit(get_val("--only", "A,B,C,D"), ",")[[1]])

results_dir <- "results"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
t0 <- Sys.time()
tic <- function(tag) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), tag))

# Expand a square cross-classification count table into a 2-column raw
# matrix: one row per subject, columns = the two raters' category codes.
expand_table <- function(M) {
  idx <- which(M > 0, arr.ind = TRUE)
  do.call(rbind, lapply(seq_len(nrow(idx)), function(k) {
    matrix(c(rep(idx[k, 1], M[idx[k, 1], idx[k, 2]]),
             rep(idx[k, 2], M[idx[k, 1], idx[k, 2]])), ncol = 2)
  }))
}

# Joint Wald test that `coefficient` is equal across same-subject fits (df =
# G - 1), via the public kappa_test() verb; also returns the mean pairwise
# correlation of the coefficient estimates (from the exposed fit$psi component).
homogeneity_wald <- function(fits, coefficient = "Conger") {
  wt   <- do.call(kappa_test, c(fits, list(coef = coefficient, paired = TRUE)))
  psi  <- do.call(cbind, lapply(fits, function(f) f$psi[, coefficient]))
  corm <- cov2cor(crossprod(psi))
  list(wt = wt, mean_corr = mean(corm[upper.tri(corm)]))
}

# ===================== Case A: within-data homogeneity ===================
run_caseA <- function() {
  tic("Case A start (anxiety, constant kappa across rater pairs)")
  data("anxiety", package = "irr")
  x <- as.matrix(anxiety)
  R <- ncol(x)
  pairs <- combn(R, 2)
  npair <- ncol(pairs)
  pair_lab <- apply(pairs, 2, function(rs) sprintf("p%d%d", rs[1], rs[2]))

  fits <- vector("list", npair)
  rows <- vector("list", npair)
  for (j in seq_len(npair)) {
    r <- pairs[1, j]; s <- pairs[2, j]
    f <- kappa(x[, c(r, s)], estimator = "nt_fiml")
    fits[[j]] <- f
    rows[[j]] <- data.frame(
      pair  = pair_lab[j], rater_1 = r, rater_2 = s,
      Conger = coef(f)[["Conger"]],
      se     = sqrt(vcov(f)["Conger", "Conger"])
    )
  }
  names(fits) <- pair_lab
  pair_df <- do.call(rbind, rows)
  write.csv(pair_df, file.path(results_dir, "caseA_pairs.csv"),
            row.names = FALSE)

  # Joint Wald test that all pairwise Conger kappas are equal (df = npair-1).
  hw <- homogeneity_wald(fits, "Conger")

  test_df <- data.frame(
    n_subjects = nrow(x), n_raters = R, n_pairs = npair,
    min_kappa  = min(pair_df$Conger), max_kappa = max(pair_df$Conger),
    range      = max(pair_df$Conger) - min(pair_df$Conger),
    mean_pairwise_corr = hw$mean_corr,
    wald_chisq = as.numeric(hw$wt$statistic),
    df         = as.numeric(hw$wt$parameter),
    p_value    = as.numeric(hw$wt$p.value)
  )
  write.csv(test_df, file.path(results_dir, "caseA_test.csv"),
            row.names = FALSE)
  tic("Case A done")
}

# ===================== Case B: independent two-sample ====================
run_caseB <- function() {
  tic("Case B start (multiple sclerosis, two independent samples)")
  # Westlund & Kurland MS data, as tabulated by Landis & Koch (1977).
  # Rows / columns = two neurologists' diagnostic class:
  # 1 = certain, 2 = probable, 3 = possible, 4 = doubtful/unlikely/not MS.
  winnipeg <- matrix(c(38, 5, 0, 1,
                       33, 11, 3, 0,
                       10, 14, 5, 6,
                        3,  7, 3, 10), nrow = 4, byrow = TRUE)
  neworleans <- matrix(c(5, 3, 0, 0,
                         3, 11, 4, 0,
                         2, 13, 3, 4,
                         1,  2, 4, 14), nrow = 4, byrow = TRUE)
  stopifnot(sum(winnipeg) == 149, sum(neworleans) == 69)

  xw <- expand_table(winnipeg)
  xn <- expand_table(neworleans)

  weights <- c("nominal", "linear", "quadratic")
  kappa_rows <- list()
  test_rows  <- list()
  for (w in weights) {
    est <- if (w == "quadratic") "pairwise" else "ipw"
    kw <- kappa(xw, estimator = est, weight = w)
    kn <- kappa(xn, estimator = est, weight = w)
    cw <- coef(kw)[["Conger"]]; sw <- sqrt(vcov(kw)["Conger", "Conger"])
    cn <- coef(kn)[["Conger"]]; sn <- sqrt(vcov(kn)["Conger", "Conger"])
    kappa_rows[[length(kappa_rows) + 1L]] <- data.frame(
      weight = w, city = "Winnipeg",   n = nrow(xw), Conger = cw, se = sw)
    kappa_rows[[length(kappa_rows) + 1L]] <- data.frame(
      weight = w, city = "New Orleans", n = nrow(xn), Conger = cn, se = sn)
    d <- cw - cn
    se <- sqrt(sw^2 + sn^2)                # independent samples: variances add
    z  <- d / se
    test_rows[[length(test_rows) + 1L]] <- data.frame(
      weight = w, kappa_winnipeg = cw, kappa_neworleans = cn,
      diff = d, se_diff = se, z = z,
      p_value = 2 * pnorm(-abs(z)))
  }
  write.csv(do.call(rbind, kappa_rows),
            file.path(results_dir, "caseB_kappa.csv"), row.names = FALSE)
  write.csv(do.call(rbind, test_rows),
            file.path(results_dir, "caseB_test.csv"), row.names = FALSE)
  tic("Case B done")
}

# ===================== Case C: FEES external benchmark ===================
run_caseC <- function() {
  tic("Case C start (FEES via multiagree)")
  have <- requireNamespace("multiagree", quietly = TRUE)
  if (!have) {
    write.csv(data.frame(ran = FALSE,
      note = paste0("multiagree not installed; FEES benchmark skipped. ",
                    "Enable with remotes::install_github('svanbelle/multiagree').")),
      file.path(results_dir, "caseC_status.csv"), row.names = FALSE)
    tic("Case C skipped (no multiagree)"); return(invisible())
  }

  e <- new.env(); utils::data(list = "FEES", package = "multiagree", envir = e)
  FEES <- get("FEES", envir = e)
  pairs <- list(p1 = c("val_TB", "val_TBR"),
                p2 = c("val_MH", "val_MHR"),
                p3 = c("val_CO", "val_COR"))
  allcols <- unlist(pairs, use.names = FALSE)

  # multiagree::delta.pair aligns all pairs to common complete cases (so the
  # dependent kappas share subjects). We reduce to the same aligned rows,
  # which lets kappa_test(paired = TRUE) stack the per-pair influence functions.
  common <- FEES[stats::complete.cases(as.matrix(FEES[, allcols])), ]
  n_aligned <- nrow(common)

  # weight-name map: misskappa <-> multiagree
  wmap <- list(nominal = "unweighted", linear = "equal", quadratic = "squared")

  kappa_rows <- list(); test_rows <- list()
  for (w in names(wmap)) {
    mg <- multiagree::delta.pair(
      data = as.matrix(common[, allcols]), cluster_id = common$subject,
      weight = wmap[[w]], multilevel = TRUE)
    est <- if (w == "quadratic") "pairwise" else "ipw"
    fits <- lapply(pairs, function(cols)
      kappa(as.matrix(common[, cols]), estimator = est, weight = w))
    names(fits) <- names(pairs)
    hw <- homogeneity_wald(fits, "Conger")

    for (i in seq_along(pairs)) {
      kappa_rows[[length(kappa_rows) + 1L]] <- data.frame(
        weight = w, pair = names(pairs)[i],
        our_Conger = coef(fits[[i]])[["Conger"]],
        our_se     = sqrt(vcov(fits[[i]])["Conger", "Conger"]),
        multiagree_kappa = mg$kappa[i, 1], multiagree_se = mg$kappa[i, 2])
    }
    test_rows[[length(test_rows) + 1L]] <- data.frame(
      weight = w, n_aligned = n_aligned,
      our_wald_chisq = as.numeric(hw$wt$statistic),
      our_df         = as.numeric(hw$wt$parameter),
      our_p          = as.numeric(hw$wt$p.value),
      multiagree_T2  = as.numeric(mg$T_test[1]),
      multiagree_p   = as.numeric(mg$T_test[2]))
  }
  write.csv(do.call(rbind, kappa_rows),
            file.path(results_dir, "caseC_kappa.csv"), row.names = FALSE)
  write.csv(do.call(rbind, test_rows),
            file.path(results_dir, "caseC_test.csv"), row.names = FALSE)
  write.csv(data.frame(ran = TRUE,
    note = sprintf(paste0("FEES via multiagree %s: %d subjects aligned across ",
                          "3 rater pairs (common complete cases)."),
                   utils::packageVersion("multiagree"), n_aligned)),
    file.path(results_dir, "caseC_status.csv"), row.names = FALSE)
  tic("Case C done")
}

# ===================== Case D: dependent two-construct test ==============
# dat.mcduff2019: MTurk judges rate smiling images on two attributes
# (smile presence, image positivity). Same judges and items, genuinely
# incomplete -- only ~4% of the item x judge grid is observed. We test
# whether the quadratically weighted Conger kappa is the same for the two
# constructs. Because both kappas are computed on the same items the
# contrast is *dependent*; kappa_test(paired = TRUE) takes the dependence from
# the per-subject influence functions -- the headline being that this runs on
# real missing data.
run_caseD <- function() {
  tic("Case D start (mcduff smile vs positive; real missing data)")
  data("dat.mcduff2019", package = "misskappa")
  d <- dat.mcduff2019
  to_wide <- function(val) {
    m <- tapply(d[[val]], list(d$item, d$judge), function(z) z[1])
    matrix(as.numeric(m), nrow = nrow(m), dimnames = dimnames(m))
  }
  constructs <- list(smile = "rating_smile", positive = "rating_positive")
  mats <- lapply(constructs, to_wide)
  fits <- lapply(mats, function(m)
    kappa(m, estimator = "pairwise"))
  names(fits) <- names(constructs)

  kappa_rows <- do.call(rbind, lapply(names(fits), function(nm) {
    m <- mats[[nm]]; f <- fits[[nm]]
    data.frame(construct = nm, n_items = nrow(m), n_judges = ncol(m),
               observed_frac = mean(!is.na(m)),
               Conger = coef(f)[["Conger"]],
               se = sqrt(vcov(f)["Conger", "Conger"]))
  }))
  write.csv(kappa_rows, file.path(results_dir, "caseD_kappa.csv"),
            row.names = FALSE)

  # Cross-covariance of the two Conger kappas from the exposed influence
  # functions; the public verb supplies the dependent and independent tests.
  ps <- fits$smile$psi[, "Conger"]; pp <- fits$positive$psi[, "Conger"]
  n  <- length(ps)
  vs  <- sum(ps^2)    / n^2
  vp  <- sum(pp^2)    / n^2
  csp <- sum(ps * pp) / n^2
  diff   <- coef(fits$smile)[["Conger"]] - coef(fits$positive)[["Conger"]]
  se_dep <- sqrt(max(vs + vp - 2 * csp, 0))
  se_ind <- sqrt(vs + vp)
  dep <- kappa_test(smile = fits$smile, positive = fits$positive,
                    coef = "Conger", paired = TRUE)
  ind <- kappa_test(smile = fits$smile, positive = fits$positive,
                    coef = "Conger", paired = FALSE)
  test_df <- data.frame(
    diff = diff, corr = csp / sqrt(vs * vp),
    se_dependent = se_dep, se_independent = se_ind,
    wald_chisq = as.numeric(dep$statistic), df = as.numeric(dep$parameter),
    p_dependent = as.numeric(dep$p.value),
    p_independent = as.numeric(ind$p.value))
  write.csv(test_df, file.path(results_dir, "caseD_test.csv"),
            row.names = FALSE)
  tic("Case D done")
}

# ---- Dispatch -----------------------------------------------------------
if ("A" %in% cases) run_caseA()
if ("B" %in% cases) run_caseB()
if ("C" %in% cases) run_caseC()
if ("D" %in% cases) run_caseD()

# ---- Metadata -----------------------------------------------------------
meta <- data.frame(
  key = c("cases", "caseA_data", "caseB_data", "caseC_data", "caseD_data",
          "R_version", "misskappa_version", "irr_version",
          "started_at", "elapsed_s"),
  value = c(
    paste(cases, collapse = ","),
    "irr::anxiety (3 raters x 20, 6-pt ordinal)",
    "Westlund-Kurland MS via Landis & Koch (1977): Winnipeg n=149, New Orleans n=69, 4 ordinal classes",
    "Vanbelle (2017) FEES via multiagree (optional)",
    "dat.mcduff2019 (273 items x 121 judges, ~4% observed, smile & positivity, 1-6 ordinal)",
    as.character(getRversion()),
    as.character(utils::packageVersion("misskappa")),
    as.character(utils::packageVersion("irr")),
    format(t0, "%Y-%m-%dT%H:%M:%S"),
    sprintf("%.1f", as.numeric(difftime(Sys.time(), t0, units = "secs")))
  )
)
write.csv(meta, file.path(results_dir, "metadata.csv"), row.names = FALSE)
cat(sprintf("Wrote %s/{caseA_*,caseB_*,caseC_*,caseD_*,metadata}.csv\n",
            results_dir))
