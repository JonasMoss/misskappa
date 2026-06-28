#!/usr/bin/env Rscript

# Reduce a cell-fan-out run into final CSVs. Each Modal cell wrote its own
# <out-dir>/cells/c<C>_<dgp>_<mech>_n<N>/{summary,replicates}.csv over disjoint
# design cells (the categories column distinguishes cats=5 tau6 from cats=7
# tau6), so the reduction is a plain row-bind. Robust to a partial run: it
# combines whatever cells are present and reports how many.
#
# Usage: Rscript combine.R --out-dir PATH

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[[i + 1L]]
}
out_dir <- get_arg("--out-dir", ".")
cells_dir <- file.path(out_dir, "cells")
if (!dir.exists(cells_dir)) stop("no cells/ dir under ", out_dir, call. = FALSE)

bind_all <- function(name) {
  files <- list.files(cells_dir, pattern = paste0("^", name, "$"),
                      recursive = TRUE, full.names = TRUE)
  if (!length(files)) {
    warning("no ", name, " found under ", cells_dir, call. = FALSE)
    return(NULL)
  }
  parts <- lapply(files, function(f) {
    tryCatch(read.csv(f, stringsAsFactors = FALSE),
             error = function(e) { warning("unreadable: ", f); NULL })
  })
  parts <- parts[!vapply(parts, is.null, logical(1))]
  list(df = do.call(rbind, parts), n_cells = length(parts))
}

s <- bind_all("summary.csv")
if (!is.null(s)) {
  write.csv(s$df, file.path(out_dir, "summary.csv"), row.names = FALSE)
  cat(sprintf("summary.csv: %d rows from %d cells\n", nrow(s$df), s$n_cells))
}
r <- bind_all("replicates.csv")
if (!is.null(r)) {
  write.csv(r$df, file.path(out_dir, "replicates.csv"), row.names = FALSE)
  cat(sprintf("replicates.csv: %d rows from %d cells\n", nrow(r$df), r$n_cells))
}

# Quick coverage check: how many of the 90 expected cells landed, and any kills.
if (!is.null(s)) {
  cells <- unique(s$df[, c("categories", "dgp", "mechanism", "n")])
  cat(sprintf("distinct design cells present: %d (expected 90)\n", nrow(cells)))
  if (all(c("n_killed", "n_killed_mem") %in% names(s$df))) {
    cat(sprintf("total fits killed: %d on time, %d on memory\n",
                sum(s$df$n_killed), sum(s$df$n_killed_mem)))
  }
}
