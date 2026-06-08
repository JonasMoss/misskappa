#!/usr/bin/env Rscript
# Experiment 19 -- alpha-overnight-regimes.
#
# The "bigger and badder" stress sibling of experiment 18: maps where each of
# the three alpha estimators (pairwise + overlap SE, normal-FIML, cat-FIML) and
# the average-n Feldt strawman break -- across extreme measurement models,
# magmaan IG/Pearson nonnormality (robustness-breaking by construction), t3
# (no 4th moment), contamination, small n, planned/zero-overlap/MNAR
# missingness, alpha near the irregular point, and the cat-FIML C^p wall.
#
# Full grid -> online appendix; a few headline cells -> manuscript.
#
# Usage:
#   Rscript run_experiment.R --help
#   Rscript run_experiment.R --smoke                 # reps=3, cheap subset
#   Rscript run_experiment.R --reps 1500 --cores 8 --progress
#   Rscript run_experiment.R --substudy s1,s3 --reps 1500
#
# Compute is here; report.qmd only reads results/.

suppressWarnings(suppressMessages({
  library(misskappa)
  has_magmaan <- requireNamespace("magmaan", quietly = TRUE)
}))

# ---- locate script dir + source modules -------------------------------------
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg))
  dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)) else getwd()
for (f in c("pairwise_se.R", "dgps.R", "generators.R", "mechanisms.R",
            "estimators.R", "metrics.R", "design.R"))
  source(file.path(script_dir, "R", f))

# ---- argument parsing -------------------------------------------------------
opts <- list(smoke = FALSE, reps = 1000L, seed_base = 20260605L,
             substudy = "all", cores = 1L, boot_B = 200L, progress = FALSE,
             out_dir = file.path(script_dir, "results"))

print_help <- function() {
  cat(
    "Experiment 19: alpha-overnight-regimes\n\n",
    "Options:\n",
    "  --smoke            Cheap wiring check: reps=3, n<=150, 2 cells/substudy, cores=1.\n",
    "  --reps N           Replicates per cell. Default 1000.\n",
    "  --substudy CSV     Subset of {s1,s2,s3,s4,s5,torture,all}. Default all.\n",
    "  --seed-base N      Base seed. Default 20260605.\n",
    "  --cores N          Parallel workers over cells (fork). Default 1.\n",
    "  --boot-B N         Case-bootstrap resamples for flagged cells. Default 200.\n",
    "  --out-dir DIR      Output directory. Default <script>/results.\n",
    "  --progress         One line per finished cell.\n",
    "  --help             This message.\n", sep = "")
}

av <- commandArgs(trailingOnly = TRUE)
i <- 1L
while (i <= length(av)) {
  a <- av[[i]]
  takes_val <- a %in% c("--reps", "--substudy", "--seed-base", "--cores", "--boot-B", "--out-dir")
  val <- if (takes_val) { i <- i + 1L; av[[i]] } else NULL
  switch(a,
    "--help"      = { print_help(); quit(status = 0) },
    "--smoke"     = opts$smoke <- TRUE,
    "--progress"  = opts$progress <- TRUE,
    "--reps"      = opts$reps <- as.integer(val),
    "--substudy"  = opts$substudy <- strsplit(val, ",")[[1]],
    "--seed-base" = opts$seed_base <- as.integer(val),
    "--cores"     = opts$cores <- as.integer(val),
    "--boot-B"    = opts$boot_B <- as.integer(val),
    "--out-dir"   = opts$out_dir <- val,
    stop("unknown argument: ", a))
  i <- i + 1L
}
if (!has_magmaan)
  message("NOTE: magmaan not installed -- IG/Pearson cells (ig_*) will be skipped.")

# ---- build (and, for smoke, subset) the design ------------------------------
design <- build_design(opts$substudy)

cell_cat_cost <- function(cell) if (grepl("cat_fiml", cell$estimators)) cell$n_cat^cell$p else 0
if (opts$smoke) {
  opts$reps <- 3L; opts$cores <- 1L; opts$boot_B <- 25L
  design$n <- pmin(design$n, 150L)
  cheap <- vapply(seq_len(nrow(design)), function(k) cell_cat_cost(design[k, ]) <= 2e4, logical(1))
  design <- design[cheap, , drop = FALSE]
  # breadth: <=2 cells per substudy, PLUS guarantee every generator/estimator
  # path (cheapest discrete cat-FIML cell, one IG cell) is exercised at least once.
  keep <- unlist(lapply(split(seq_len(nrow(design)), design$substudy), function(ix) head(ix, 2L)))
  disc <- which(design$dist == "discrete")
  if (length(disc))
    keep <- c(keep, disc[which.min(vapply(disc, function(k) cell_cat_cost(design[k, ]), numeric(1)))])
  ig <- which(design$dist %in% c("ig_skew", "ig_heavy"))
  if (length(ig)) keep <- c(keep, ig[1])
  design <- design[sort(unique(keep)), , drop = FALSE]
}
# IG cells need magmaan; drop them if unavailable.
if (!has_magmaan) design <- design[!design$dist %in% c("ig_skew", "ig_heavy"), , drop = FALSE]

# ---- per-cell driver --------------------------------------------------------
run_cell <- function(cell) {
  Sigma <- cell_sigma(cell$model, cell$p, cell$alpha_target)
  truth <- cell_truth_alpha(cell, Sigma)
  cal <- if (cell$dist %in% c("ig_skew", "ig_heavy")) {
    tg <- ig_targets(cell$dist); ig_calibrate(Sigma, tg$skew, tg$ek)
  } else NULL
  ests <- strsplit(cell$estimators, ",")[[1]]

  reps <- lapply(seq_len(opts$reps), function(rep) {
    rep_seed <- opts$seed_base + cell$cell_id * 100000L + rep
    X <- tryCatch(generate_complete(cell, Sigma, rep_seed, cal), error = function(e) NULL)
    if (is.null(X)) return(NULL)
    set.seed(rep_seed + 777777L)
    X <- apply_mechanism(X, cell$mechanism)
    do.call(rbind, lapply(ests, function(en)
      run_estimator(en, X, do_boot = (en == "pairwise" && cell$bootstrap), boot_B = opts$boot_B)))
  })
  reps_df <- do.call(rbind, reps)

  summ <- summarise_cell(reps_df, truth)
  design_cols <- cell[rep(1L, nrow(summ)),
    c("cell_id", "substudy", "label", "model", "p", "n", "alpha_target",
      "dist", "n_cat", "mechanism")]
  summ <- cbind(design_cols, summ[setdiff(names(summ), "truth_alpha")], truth_alpha = summ$truth_alpha,
                row.names = NULL)
  if (opts$progress) cat(sprintf("[%s] cell %d (%s/%s) done\n",
    format(Sys.time(), "%H:%M:%S"), cell$cell_id, cell$substudy, cell$label))
  list(summary = summ,
       truth = data.frame(cell[1, c("cell_id", "substudy", "label", "model", "p", "n",
         "alpha_target", "dist", "n_cat", "mechanism")], truth_alpha = truth,
         reps = opts$reps, row.names = NULL))
}

cells <- lapply(seq_len(nrow(design)), function(k) design[k, , drop = FALSE])
res <- if (opts$cores > 1L) {
  parallel::mclapply(cells, run_cell, mc.cores = opts$cores, mc.preschedule = FALSE)
} else {
  lapply(cells, run_cell)
}

# ---- collect + write --------------------------------------------------------
ok <- !vapply(res, function(r) inherits(r, "try-error") || is.null(r), logical(1))
summary_df <- do.call(rbind, lapply(res[ok], `[[`, "summary"))
truth_df   <- do.call(rbind, lapply(res[ok], `[[`, "truth"))

dir.create(opts$out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(summary_df, file.path(opts$out_dir, "summary.csv"), row.names = FALSE)
write.csv(truth_df,   file.path(opts$out_dir, "truth.csv"),   row.names = FALSE)

metadata <- data.frame(
  key = c("mode", "reps", "substudy", "seed_base", "cores", "boot_B", "n_cells",
          "magmaan", "misskappa_version", "magmaan_version", "R_version", "run_time"),
  value = c(if (opts$smoke) "smoke" else "full", opts$reps,
            paste(opts$substudy, collapse = ","), opts$seed_base, opts$cores, opts$boot_B,
            nrow(design), has_magmaan,
            as.character(utils::packageVersion("misskappa")),
            if (has_magmaan) as.character(utils::packageVersion("magmaan")) else NA,
            as.character(getRversion()), format(Sys.time())),
  stringsAsFactors = FALSE)
write.csv(metadata, file.path(opts$out_dir, "metadata.csv"), row.names = FALSE)

cat(sprintf("Wrote %d summary rows over %d cells to %s\n",
            nrow(summary_df), sum(ok), opts$out_dir))
