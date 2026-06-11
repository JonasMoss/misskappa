#!/usr/bin/env Rscript
#
# Merge per-cell replicate checkpoints across replicate shards into the final
# replicates.csv + summary.csv. run_experiment.R shards by REPLICATE: each shard
# computes every design cell but only its replicate subset, writing
# <out-dir>/checkpoints/cell-XXX-<dgp>-<mechanism>-nN-shardMofN.csv. A single
# shard's end-of-run summary therefore covers only its own replicates; this
# script reduces over EVERY shard's checkpoints at once so the merged summary is
# over the full replicate set. Use after a sharded (e.g. Modal) run.
#
# Usage: Rscript combine.R --out-dir PATH [--checkpoint-dir PATH]

script_dir <- local({
  cmd <- commandArgs(FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) == 0L) getwd()
  else dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE))
})
source(file.path(script_dir, "summarize.R"))

argv <- commandArgs(trailingOnly = TRUE)
get_opt <- function(flag, default = NULL) {
  i <- match(flag, argv)
  if (is.na(i)) return(default)
  if (i == length(argv)) stop(flag, " needs a value.", call. = FALSE)
  argv[[i + 1L]]
}

out_dir <- get_opt("--out-dir")
if (is.null(out_dir)) stop("combine.R requires --out-dir PATH", call. = FALSE)
checkpoint_dir <- get_opt("--checkpoint-dir", file.path(out_dir, "checkpoints"))

files <- list.files(checkpoint_dir, pattern = "^cell-.*\\.csv$", full.names = TRUE)
if (!length(files)) {
  stop("No cell-*.csv checkpoints found in ", checkpoint_dir, call. = FALSE)
}
message(sprintf("Merging %d checkpoint file(s) from %s", length(files), checkpoint_dir))

replicates <- do.call(rbind, lapply(files, function(p) {
  df <- read.csv(p, stringsAsFactors = FALSE)
  if (!nrow(df)) NULL else df
}))
if (is.null(replicates) || !nrow(replicates)) {
  stop("All checkpoint files were empty.", call. = FALSE)
}

## Replicate rows are deterministic in (cell, rep); distinct shards hold disjoint
## rep ids, so rows are unique. Drop any byte-identical duplicates defensively
## (e.g. a non-sharded shard01of01 run later re-sharded into the same dir).
n_before <- nrow(replicates)
replicates <- unique(replicates)
if (nrow(replicates) < n_before) {
  message(sprintf("Dropped %d duplicate replicate row(s).", n_before - nrow(replicates)))
}

summary <- summarize_replicates(replicates)
write_csv_atomic(replicates, file.path(out_dir, "replicates.csv"))
write_csv_atomic(summary, file.path(out_dir, "summary.csv"))
message(sprintf(
  "Wrote %s (%d rows) and %s (%d rows).",
  file.path(out_dir, "replicates.csv"), nrow(replicates),
  file.path(out_dir, "summary.csv"), nrow(summary)
))
