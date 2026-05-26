#!/usr/bin/env Rscript
#
# Regenerates the main + appendix simulation tables for kappa-missing.tex.
# Writes:
#   paper/tables/simulations-raw-three-estimators-main.tex     (Conger, identity)
#   paper/supplement/simulations-raw-three-estimators-appendix.tex (Fleiss, BP, linear loss)
#   paper/results/simulations-raw-three-estimators-summary.csv (curated)
#
# Run from the paper/ directory:
#   Rscript scripts/simulations_raw_three_estimators.R
#
# Default mode is "small" (smoke / quick iteration). Set SIM_FULL=1 for the
# Psychometrika-sized run (slow: tens of minutes to hours).

suppressPackageStartupMessages({
  library(misskappa)
})

fmt <- function(x, digits = 4) formatC(x, format = "f", digits = digits)
fmt_bias_sd <- function(b, s, digits = 4) paste0(fmt(b, digits), " (", fmt(s, digits), ")")

set.seed(1)

# --- Global design ---
C <- 5L
R <- 6L

# Marginal distribution p: empirical marginals from Fleiss (1971).
p <- colSums(as.matrix(dat.fleiss1971)) / sum(as.matrix(dat.fleiss1971))

# Sample sizes / replication counts.
sim_full <- Sys.getenv("SIM_FULL", "0") == "1"
if (sim_full) {
  n_truth <- 200000L
  n_mod   <- 4000L
  n_big   <- 40000L
  B_mod   <- 500L
  B_big   <- 200L
} else {
  n_truth <- 5000L
  n_mod   <- 200L
  n_big   <- 1000L
  B_mod   <- 20L
  B_big   <- 8L
}

methods <- c("available", "ipw", "gwet")
weight_main <- "identity"
weight_appendix <- "linear"

guess_mat <- matrix(rep(p, each = R), nrow = R, byrow = TRUE)

# --- Simulators ---
simulate_ratings <- function(n, rho_base, rho_truth_mult, guess) {
  truth <- sample.int(C, n, replace = TRUE, prob = p)
  x_star <- matrix(NA_integer_, nrow = n, ncol = R)
  for (j in seq_len(R)) {
    rho_i <- pmin(pmax(rho_base[j] * rho_truth_mult[truth], 0), 1)
    correct <- stats::rbinom(n, 1, prob = rho_i) == 1
    guessed <- sample.int(C, n, replace = TRUE, prob = guess[j, ])
    x_star[, j] <- ifelse(correct, truth, guessed)
  }
  list(x_star = x_star, truth = truth)
}

apply_missing_mcar <- function(x_star, pi_rater) {
  n <- nrow(x_star)
  M <- matrix(FALSE, nrow = n, ncol = R)
  for (j in seq_len(R)) M[, j] <- stats::runif(n) < pi_rater[j]
  x <- x_star
  x[!M] <- NA_integer_
  x
}

apply_missing_mar_truth <- function(x_star, truth, pi_truth) {
  n <- nrow(x_star)
  p_i <- pi_truth[truth]
  M <- matrix(stats::runif(n * R) < rep(p_i, times = R), nrow = n, ncol = R)
  x <- x_star
  x[!M] <- NA_integer_
  x
}

kappa_all <- function(x, method, weight) {
  misskappa::kappa(x, method = method, weight = weight)$estimates
}

estimate_truth <- function(dgp, weight) {
  x_star <- simulate_ratings(
    n_truth,
    rho_base = dgp$rho_base,
    rho_truth_mult = dgp$rho_truth_mult,
    guess = dgp$guess
  )$x_star
  misskappa::kappa(x_star, method = "available", weight = weight)$estimates
}

run_once <- function(n, dgp, weight) {
  dat <- simulate_ratings(
    n,
    rho_base = dgp$rho_base,
    rho_truth_mult = dgp$rho_truth_mult,
    guess = dgp$guess
  )
  x <- switch(
    dgp$missing,
    "mcar"      = apply_missing_mcar(dat$x_star, dgp$pi_rater),
    "mar_truth" = apply_missing_mar_truth(dat$x_star, dat$truth, dgp$pi_truth),
    stop("Unknown missingness mechanism: ", dgp$missing)
  )
  out <- lapply(methods, function(m) kappa_all(x, method = m, weight = weight))
  names(out) <- methods
  out
}

summarize_kappas <- function(res_list, kappa_name, truth_val) {
  est <- vapply(res_list, function(x) as.numeric(x[[kappa_name]]), numeric(1))
  c(mean = mean(est, na.rm = TRUE),
    sd   = stats::sd(est, na.rm = TRUE),
    bias = mean(est, na.rm = TRUE) - truth_val)
}

extract_method <- function(rep_list, method) {
  lapply(rep_list, function(r) r[[method]])
}

make_summary_rows <- function(n, B, dgp, truth_estimates, weight, kappas) {
  reps <- replicate(B, run_once(n, dgp, weight = weight), simplify = FALSE)
  rows <- lapply(methods, function(m) {
    per_rep <- extract_method(reps, m)
    out <- list(method = m)
    if ("Conger" %in% kappas)
      out$conger <- summarize_kappas(per_rep, "Conger", truth_estimates[["Conger"]])
    if ("Fleiss" %in% kappas)
      out$fleiss <- summarize_kappas(per_rep, "Fleiss", truth_estimates[["Fleiss"]])
    if ("Brennan-Prediger" %in% kappas)
      out$bp <- summarize_kappas(per_rep, "Brennan-Prediger",
                                 truth_estimates[["Brennan-Prediger"]])
    out
  })
  names(rows) <- methods
  rows
}

# --- DGPs ---
dgpA <- list(
  label = "A", name = "Exchangeable + MCAR",
  rho_base = rep(0.92, R), rho_truth_mult = rep(1, C), guess = guess_mat,
  missing = "mcar", pi_rater = c(0.95, 0.90, 0.80, 0.65, 0.50, 0.35))

dgpB <- list(
  label = "B", name = "Non-exchangeable + MCAR",
  rho_base = c(0.98, 0.96, 0.92, 0.85, 0.72, 0.60),
  rho_truth_mult = rep(1, C), guess = guess_mat,
  missing = "mcar", pi_rater = c(0.95, 0.85, 0.70, 0.45, 0.25, 0.15))

dgpC <- list(
  label = "C", name = "Difficulty + MAR",
  rho_base = rep(0.92, R), rho_truth_mult = c(1.05, 1.00, 0.95, 0.85, 0.70),
  guess = guess_mat,
  missing = "mar_truth", pi_truth = c(0.95, 0.90, 0.80, 0.55, 0.25))

dgps <- list(dgpA, dgpB, dgpC)

# --- Run ---
truth_main <- setNames(lapply(dgps, estimate_truth, weight = weight_main),
                       vapply(dgps, `[[`, character(1), "label"))
truth_app  <- setNames(lapply(dgps, estimate_truth, weight = weight_appendix),
                       vapply(dgps, `[[`, character(1), "label"))

summaries_main_mod <- lapply(dgps, function(d)
  make_summary_rows(n_mod, B_mod, d, truth_main[[d$label]], weight_main,
                    c("Conger", "Fleiss", "Brennan-Prediger")))
summaries_main_big <- lapply(dgps, function(d)
  make_summary_rows(n_big, B_big, d, truth_main[[d$label]], weight_main,
                    c("Conger", "Fleiss", "Brennan-Prediger")))
summaries_app_mod  <- lapply(dgps, function(d)
  make_summary_rows(n_mod, B_mod, d, truth_app[[d$label]], weight_appendix,
                    c("Conger")))
summaries_app_big  <- lapply(dgps, function(d)
  make_summary_rows(n_big, B_big, d, truth_app[[d$label]], weight_appendix,
                    c("Conger")))

# --- Output paths (script lives in paper/scripts/, runs with cwd = paper/) ---
tables_dir     <- "tables"
supplement_dir <- "supplement"
results_dir    <- "results"
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supplement_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

method_labels <- c(available = "AC", ipw = "IPW", gwet = "Gwet")

latex_table_conger <- function(summaries_mod, summaries_big, truth_lookup) {
  hdr <- c("\\begin{tabular}{llrrrr}", "\\toprule",
           paste0("DGP & Method & $\\kappa_C$ (truth) & Bias (SD), $n=", n_mod,
                  "$ & Bias (SD), $n=", n_big, "$\\\\"),
           "\\midrule")
  lns <- character()
  for (idx in seq_along(dgps)) {
    d <- dgps[[idx]]
    tC <- as.numeric(truth_lookup[[d$label]][["Conger"]])
    for (m in methods) {
      smod <- summaries_mod[[idx]][[m]]$conger
      sbig <- summaries_big[[idx]][[m]]$conger
      lns <- c(lns, paste0(d$label, " & ", method_labels[[m]], " & ",
                           fmt(tC, 4), " & ",
                           fmt_bias_sd(smod[["bias"]], smod[["sd"]], 4), " & ",
                           fmt_bias_sd(sbig[["bias"]], sbig[["sd"]], 4), "\\\\"))
    }
    if (idx != length(dgps)) lns <- c(lns, "\\addlinespace")
  }
  c(hdr, lns, "\\bottomrule", "\\end{tabular}")
}

latex_table_appendix <- function(kappa_name, key) {
  hdr <- c("\\begin{tabular}{llrrrr}", "\\toprule",
           paste0("DGP & Method & truth & Bias (SD), $n=", n_mod,
                  "$ & Bias (SD), $n=", n_big, "$\\\\"),
           "\\midrule")
  lns <- character()
  for (idx in seq_along(dgps)) {
    d <- dgps[[idx]]
    tK <- as.numeric(truth_main[[d$label]][[kappa_name]])
    for (m in methods) {
      smod <- summaries_main_mod[[idx]][[m]][[key]]
      sbig <- summaries_main_big[[idx]][[m]][[key]]
      lns <- c(lns, paste0(d$label, " & ", method_labels[[m]], " & ",
                           fmt(tK, 4), " & ",
                           fmt_bias_sd(smod[["bias"]], smod[["sd"]], 4), " & ",
                           fmt_bias_sd(sbig[["bias"]], sbig[["sd"]], 4), "\\\\"))
    }
    if (idx != length(dgps)) lns <- c(lns, "\\addlinespace")
  }
  c(hdr, lns, "\\bottomrule", "\\end{tabular}")
}

note_dgps <- function() {
  paste(
    "\\textit{Note.} Rating model: $T_i\\sim\\mathrm{Categorical}(p)$; conditional on $T_i$, rater $j$ reports $T_i$ with probability $\\rho_{ij}$ and otherwise guesses from $p$. Missingness indicators $M_{ij}\\in\\{0,1\\}$ are independent across $(i,j)$ given the specified mechanism. DGP A: exchangeable raters, MCAR rater-specific. DGP B: non-exchangeable, MCAR rater-specific. DGP C: difficulty-dependent skill, MAR via the truth $T_i$. AC = available-case, IPW = inverse-probability-weighted, Gwet = Gwet's reweighting.",
    collapse = " ")
}

snippet_main <- c(
  "% Auto-generated by paper/scripts/simulations_raw_three_estimators.R",
  "\\begin{table}[!ht]\\centering\\small\\setlength{\\tabcolsep}{4pt}",
  "\\caption{Simulation bias and Monte Carlo SD for Conger's $\\kappa_C$ (nominal/identity loss).}\\label{tab:sim-conger}",
  latex_table_conger(summaries_main_mod, summaries_main_big, truth_main),
  "\\vspace{3pt}\\begin{minipage}{0.95\\linewidth}\\footnotesize", note_dgps(),
  "\\end{minipage}\\end{table}")
writeLines(snippet_main, file.path(tables_dir, "simulations-raw-three-estimators-main.tex"))

snippet_appendix <- c(
  "% Auto-generated by paper/scripts/simulations_raw_three_estimators.R",
  "\\subsubsection*{Conger's kappa (linear/absolute loss)}",
  "\\begin{table}[!ht]\\centering\\small\\setlength{\\tabcolsep}{4pt}",
  "\\caption{Simulation bias and Monte Carlo SD for Conger's $\\kappa_C$ (linear/absolute loss).}\\label{tab:sim-conger-linear}",
  latex_table_conger(summaries_app_mod, summaries_app_big, truth_app),
  "\\end{table}",
  "",
  "\\subsubsection*{Fleiss' kappa (identity loss)}",
  "\\begin{center}\\small\\setlength{\\tabcolsep}{4pt}",
  latex_table_appendix("Fleiss", "fleiss"),
  "\\end{center}",
  "",
  "\\subsubsection*{Brennan--Prediger (identity loss)}",
  "\\begin{center}\\small\\setlength{\\tabcolsep}{4pt}",
  latex_table_appendix("Brennan-Prediger", "bp"),
  "\\end{center}")
writeLines(snippet_appendix,
           file.path(supplement_dir, "simulations-raw-three-estimators-appendix.tex"))

# Curated summary CSV.
rows <- list()
for (idx in seq_along(dgps)) {
  d <- dgps[[idx]]
  for (m in methods) {
    smod <- summaries_main_mod[[idx]][[m]]
    sbig <- summaries_main_big[[idx]][[m]]
    for (kn in c("Conger", "Fleiss", "Brennan-Prediger")) {
      key <- c(Conger = "conger", Fleiss = "fleiss", `Brennan-Prediger` = "bp")[kn]
      rows[[length(rows) + 1]] <- data.frame(
        dgp = d$label, method = m, kappa = kn, weight = weight_main,
        n_mod = n_mod, bias_mod = smod[[key]][["bias"]], sd_mod = smod[[key]][["sd"]],
        n_big = n_big, bias_big = sbig[[key]][["bias"]], sd_big = sbig[[key]][["sd"]],
        truth = as.numeric(truth_main[[d$label]][[kn]]))
    }
  }
}
write.csv(do.call(rbind, rows),
          file.path(results_dir, "simulations-raw-three-estimators-summary.csv"),
          row.names = FALSE)

cat(sprintf("Wrote %s/simulations-raw-three-estimators-main.tex\n", tables_dir))
cat(sprintf("Wrote %s/simulations-raw-three-estimators-appendix.tex\n", supplement_dir))
cat(sprintf("Wrote %s/simulations-raw-three-estimators-summary.csv\n", results_dir))
