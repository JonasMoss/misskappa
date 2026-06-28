#!/usr/bin/env Rscript

# Parallel launcher for the three-sim alpha design (probe 46 / two-sim plan):
#   Sim A: 5-category x 4-item  (5^4 = 625)     all four estimators
#   Sim B: 5-category x 6-item  (5^6 = 15625)   all four + stall-killer
#   Sim C: 7-category x 6-item  (7^6 = 117649)  all four + stall-killer (the
#          realistic 7-point "challenge": cat-FIML kills become visible at large n)
#
# Each (sim, dgp, mechanism) chunk is a separate `run_experiment.R` process, run
# up to `--cores` at a time via mclapply with dynamic scheduling. The per-fit
# stall-killer lives INSIDE each runner (its own mcparallel), so there is no
# nested forking in one R session -- chunks are independent OS processes.
#
# Usage: Rscript run_3sims.R [B=3000] [out_base] [cores=4] [ngrid] [timeout=30]

suppressMessages(library(parallel))

argv <- commandArgs(trailingOnly = TRUE)
arg <- function(i, d) if (length(argv) >= i && nzchar(argv[i])) argv[i] else d

script_args <- commandArgs(FALSE)
this_file <- sub("^--file=", "", script_args[grep("^--file=", script_args)][1L])
script_dir <- dirname(normalizePath(this_file))
runner <- file.path(script_dir, "run_experiment.R")

B        <- as.integer(arg(1, "3000"))
out_base <- arg(2, file.path(script_dir, "results", "3sims"))
cores    <- as.integer(arg(3, "4"))
ngrid    <- arg(4, "100,250,500,1000,2000")
timeout  <- arg(5, "30")
mem_mb   <- arg(6, "2500")   # per-fit RSS budget; kills runaway cat-FIML fits
methods  <- "pairwise,cat_em,normal_fiml,listwise"
mechs    <- c("mcar30", "mar_anchor30")

sims <- list(
  list(name = "simA", cats = 5L,
       dgps = c("ordinal_tau4", "ordinal_congeneric4", "ordinal_twofactor4")),
  list(name = "simB", cats = 5L,
       dgps = c("ordinal_tau6", "ordinal_congeneric6", "ordinal_twofactor6")),
  list(name = "simC", cats = 7L,
       dgps = c("ordinal_tau6", "ordinal_congeneric6", "ordinal_twofactor6"))
)

# Resumable: a chunk that already wrote summary.csv is skipped, so a relaunch
# continues where an interrupted run left off. ci is the stable chunk index
# (drives the per-chunk seed), incremented even for skipped chunks so seeds
# never shift on resume.
chunks <- list()
ci <- 0L
skipped <- 0L
for (s in sims) for (dgp in s$dgps) for (mech in mechs) {
  ci <- ci + 1L
  od <- file.path(out_base, s$name, paste0(dgp, "__", mech))
  if (file.exists(file.path(od, "summary.csv"))) {
    skipped <- skipped + 1L
    next
  }
  dir.create(od, recursive = TRUE, showWarnings = FALSE)
  chunks[[length(chunks) + 1L]] <- list(
    tag = paste(s$name, dgp, mech),
    log = paste0(od, ".log"),
    args = c(runner,
             "--reps", B, "--n-grid", ngrid, "--categories", s$cats,
             "--dgps", dgp, "--mechanisms", mech, "--methods", methods,
             "--fit-timeout", timeout, "--fit-memory-mb", mem_mb,
             "--seed-base", ci * 10000000L, "--out-dir", od))
}
if (skipped > 0L) cat(sprintf("Resuming: %d chunk(s) already complete, skipping.\n", skipped))

cat(sprintf("Launching %d chunks, B=%d, %d cores, n-grid=%s, timeout=%ss, mem=%sMB/fit\n",
            length(chunks), B, cores, ngrid, timeout, mem_mb))
t0 <- Sys.time()
res <- mclapply(chunks, function(ch) {
  st <- Sys.time()
  status <- system2("Rscript", as.character(ch$args),
                    stdout = ch$log, stderr = ch$log)
  data.frame(tag = ch$tag, exit = status,
             mins = round(as.numeric(difftime(Sys.time(), st, units = "mins")), 2),
             stringsAsFactors = FALSE)
}, mc.cores = cores, mc.preschedule = FALSE)

tab <- do.call(rbind, res)
tab <- tab[order(-tab$mins), ]
print(tab, row.names = FALSE)
cat(sprintf("\nTotal wall: %.1f min (%.2f h).  Output under %s\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")),
            as.numeric(difftime(Sys.time(), t0, units = "hours")), out_base))
