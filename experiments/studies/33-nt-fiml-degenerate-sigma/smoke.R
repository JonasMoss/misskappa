#!/usr/bin/env Rscript

args0 <- commandArgs(FALSE)
file_arg <- grep("^--file=", args0, value = TRUE)
self <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath("experiments/studies/33-nt-fiml-degenerate-sigma/smoke.R",
                mustWork = TRUE)
}
script <- normalizePath(file.path(dirname(self), "run_experiment.R"), mustWork = TRUE)
args <- c("--smoke", "--load-all", commandArgs(TRUE))
sys <- system2("Rscript", c(script, args))
quit(status = sys)
