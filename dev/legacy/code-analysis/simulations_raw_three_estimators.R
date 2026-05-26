#!/usr/bin/env Rscript

suppressWarnings({
  .libPaths(c(file.path("code", "analysis", "Rlib"), .libPaths()))
})

suppressPackageStartupMessages({
  library(misskappa)
})

fmt <- function(x, digits = 4) formatC(x, format = "f", digits = digits)
fmt_bias_sd <- function(bias, sd, digits = 4) paste0(fmt(bias, digits), " (", fmt(sd, digits), ")")

set.seed(1)

# Require the local misskappa API (available/ipw/gwet). If a different misskappa
# is loaded from the global library path, instruct how to install the repo copy.
method_choices <- eval(formals(misskappa::kappa_raw)$method)
required_methods <- c("available", "ipw", "gwet")
if (!all(required_methods %in% method_choices)) {
  stop(
    "This simulation script requires misskappa::kappa_raw(method=available/ipw/gwet).\n",
    "Install the repo version into the local analysis library with:\n",
    "  R CMD INSTALL -l code/analysis/Rlib code/misskappa\n"
  )
}

# --- Global design ---
C <- 5L
R <- 6L

# Baseline marginal distribution p: use the empirical marginal distribution from
# Fleiss (1971) count-form data as a convenient, reproducible default.
data("dat.fleiss1971", package = "misskappa", envir = environment())
p <- colSums(as.matrix(dat.fleiss1971)) / sum(as.matrix(dat.fleiss1971))

# Sample sizes and replication counts (kept modest; increase later if needed).
if (Sys.getenv("SIM_FULL", "0") == "1") {
  n_truth <- 200000L
  n_mod <- 4000L
  n_big <- 40000L
  B_mod <- 500L
  B_big <- 200L
} else {
  n_truth <- 30000L
  n_mod <- 600L
  n_big <- 6000L
  B_mod <- 40L
  B_big <- 15L
}

methods <- c("available", "ipw", "gwet")
weight_main <- "identity"
weight_appendix <- "linear"

guess_mat <- matrix(rep(p, each = R), nrow = R, byrow = TRUE)

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
  for (j in seq_len(R)) {
    M[, j] <- stats::runif(n) < pi_rater[j]
  }
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

kappa_all <- function(x_incomplete, method, weight) {
  misskappa::kappa_raw(x_incomplete, method = method, weight = weight)$estimates
}

estimate_truth <- function(dgp, weight) {
  x_star <- simulate_ratings(
    n_truth,
    rho_base = dgp$rho_base,
    rho_truth_mult = dgp$rho_truth_mult,
    guess = dgp$guess
  )$x_star
  misskappa::kappa_raw(x_star, method = "available", weight = weight)$estimates
}

run_once <- function(n, dgp, weight) {
  dat <- simulate_ratings(
    n,
    rho_base = dgp$rho_base,
    rho_truth_mult = dgp$rho_truth_mult,
    guess = dgp$guess
  )
  x_star <- dat$x_star
  truth <- dat$truth

  x <- switch(
    dgp$missing,
    "mcar" = apply_missing_mcar(x_star, dgp$pi_rater),
    "mar_truth" = apply_missing_mar_truth(x_star, truth, dgp$pi_truth),
    stop("Unknown missingness mechanism: ", dgp$missing)
  )

  out <- lapply(methods, function(m) kappa_all(x, method = m, weight = weight))
  names(out) <- methods
  out
}

summarize_kappas <- function(res_list, kappa_name, truth_val) {
  est <- vapply(res_list, function(x) as.numeric(x[[kappa_name]]), numeric(1))
  c(
    mean = mean(est, na.rm = TRUE),
    sd = stats::sd(est, na.rm = TRUE),
    bias = mean(est, na.rm = TRUE) - truth_val
  )
}

extract_method <- function(rep_list, method) lapply(rep_list, function(one_rep) one_rep[[method]])

make_summary_rows <- function(n, B, dgp, truth_estimates, weight, kappas) {
  reps <- replicate(B, run_once(n, dgp, weight = weight), simplify = FALSE)
  rows <- lapply(methods, function(m) {
    per_rep <- extract_method(reps, m)
    out <- list(method = m)
    if ("Conger" %in% kappas) out$conger <- summarize_kappas(per_rep, "Conger", truth_estimates[["Conger"]])
    if ("Fleiss" %in% kappas) out$fleiss <- summarize_kappas(per_rep, "Fleiss", truth_estimates[["Fleiss"]])
    if ("Brennan-Prediger" %in% kappas) out$bp <- summarize_kappas(per_rep, "Brennan-Prediger", truth_estimates[["Brennan-Prediger"]])
    out
  })
  names(rows) <- methods
  rows
}

# --- DGP definitions (ratings + missingness) ---
# Rating model: T_i ~ Categorical(p); conditional on T_i, rater j reports T_i with
# prob rho_{ij} and otherwise guesses from p.
#
# Missingness indicators M_{ij} are applied to the observed matrix X_{ij}:
# - MCAR: P(M_{ij}=1)=pi_j
# - MAR(T): P(M_{ij}=1 | T_i=c)=pi_c
dgpA <- list(
  label = "A",
  name = "Exchangeable + MCAR",
  rho_base = rep(0.92, R),
  rho_truth_mult = rep(1, C),
  guess = guess_mat,
  missing = "mcar",
  pi_rater = c(0.95, 0.90, 0.80, 0.65, 0.50, 0.35)
)

dgpB <- list(
  label = "B",
  name = "Non-exchangeable + MCAR",
  rho_base = c(0.98, 0.96, 0.92, 0.85, 0.72, 0.60),
  rho_truth_mult = rep(1, C),
  guess = guess_mat,
  missing = "mcar",
  pi_rater = c(0.95, 0.85, 0.70, 0.45, 0.25, 0.15)
)

dgpC <- list(
  label = "C",
  name = "Difficulty + MAR",
  rho_base = rep(0.92, R),
  rho_truth_mult = c(1.05, 1.00, 0.95, 0.85, 0.70),
  guess = guess_mat,
  missing = "mar_truth",
  pi_truth = c(0.95, 0.90, 0.80, 0.55, 0.25)
)

dgps <- list(dgpA, dgpB, dgpC)

# --- Run study ---
truth_main <- lapply(dgps, function(dgp) estimate_truth(dgp, weight = weight_main))
truth_app <- lapply(dgps, function(dgp) estimate_truth(dgp, weight = weight_appendix))
names(truth_main) <- vapply(dgps, `[[`, character(1), "label")
names(truth_app) <- vapply(dgps, `[[`, character(1), "label")

summaries_main_mod <- lapply(dgps, function(dgp) {
  make_summary_rows(
    n = n_mod,
    B = B_mod,
    dgp = dgp,
    truth_estimates = truth_main[[dgp$label]],
    weight = weight_main,
    kappas = c("Conger", "Fleiss", "Brennan-Prediger")
  )
})
summaries_main_big <- lapply(dgps, function(dgp) {
  make_summary_rows(
    n = n_big,
    B = B_big,
    dgp = dgp,
    truth_estimates = truth_main[[dgp$label]],
    weight = weight_main,
    kappas = c("Conger", "Fleiss", "Brennan-Prediger")
  )
})

summaries_app_mod <- lapply(dgps, function(dgp) {
  make_summary_rows(
    n = n_mod,
    B = B_mod,
    dgp = dgp,
    truth_estimates = truth_app[[dgp$label]],
    weight = weight_appendix,
    kappas = c("Conger")
  )
})
summaries_app_big <- lapply(dgps, function(dgp) {
  make_summary_rows(
    n = n_big,
    B = B_big,
    dgp = dgp,
    truth_estimates = truth_app[[dgp$label]],
    weight = weight_appendix,
    kappas = c("Conger")
  )
})

# --- Write a compact LaTeX note ---
out_tex <- file.path("notes", "simulations-raw-three-estimators.tex")
out_snippet_main <- file.path("notes", "simulations-raw-three-estimators-main.tex")
out_snippet_appendix <- file.path("notes", "simulations-raw-three-estimators-appendix.tex")
dir.create(dirname(out_tex), showWarnings = FALSE, recursive = TRUE)

method_labels <- c(available = "AC", ipw = "IPW", gwet = "Gwet")

latex_table_conger <- function() {
  header <- c(
    "\\begin{tabular}{llrrrr}",
    "\\toprule",
    paste0(
      "DGP & Method & $\\kappa_C$ (truth) & Bias (SD), $n=", n_mod, "$ & Bias (SD), $n=", n_big, "$\\\\"
    ),
    "\\midrule"
  )

  lines <- character()
  for (idx in seq_along(dgps)) {
    dgp <- dgps[[idx]]
    tC <- as.numeric(truth_main[[dgp$label]][["Conger"]])
    for (m in methods) {
      s_mod <- summaries_main_mod[[idx]][[m]]$conger
      s_big <- summaries_main_big[[idx]][[m]]$conger
      lines <- c(
        lines,
        paste0(
          dgp$label, " & ",
          method_labels[[m]], " & ",
          fmt(tC, 4), " & ",
          fmt_bias_sd(s_mod[["bias"]], s_mod[["sd"]], 4), " & ",
          fmt_bias_sd(s_big[["bias"]], s_big[["sd"]], 4), "\\\\"
        )
      )
    }
    if (idx != length(dgps)) lines <- c(lines, "\\addlinespace")
  }

  footer <- c("\\bottomrule", "\\end{tabular}")
  c(header, lines, footer)
}

latex_table_conger_linear <- function() {
  header <- c(
    "\\begin{tabular}{llrrrr}",
    "\\toprule",
    paste0(
      "DGP & Method & $\\kappa_C$ (truth) & Bias (SD), $n=", n_mod, "$ & Bias (SD), $n=", n_big, "$\\\\"
    ),
    "\\midrule"
  )

  lines <- character()
  for (idx in seq_along(dgps)) {
    dgp <- dgps[[idx]]
    tC <- as.numeric(truth_app[[dgp$label]][["Conger"]])
    for (m in methods) {
      s_mod <- summaries_app_mod[[idx]][[m]]$conger
      s_big <- summaries_app_big[[idx]][[m]]$conger
      lines <- c(
        lines,
        paste0(
          dgp$label, " & ",
          method_labels[[m]], " & ",
          fmt(tC, 4), " & ",
          fmt_bias_sd(s_mod[["bias"]], s_mod[["sd"]], 4), " & ",
          fmt_bias_sd(s_big[["bias"]], s_big[["sd"]], 4), "\\\\"
        )
      )
    }
    if (idx != length(dgps)) lines <- c(lines, "\\addlinespace")
  }

  footer <- c("\\bottomrule", "\\end{tabular}")
  c(header, lines, footer)
}

latex_table_appendix_kappa <- function(kappa_name, key) {
  header <- c(
    "\\begin{tabular}{llrrrr}",
    "\\toprule",
    paste0(
      "DGP & Method & truth & Bias (SD), $n=", n_mod, "$ & Bias (SD), $n=", n_big, "$\\\\"
    ),
    "\\midrule"
  )

  lines <- character()
  for (idx in seq_along(dgps)) {
    dgp <- dgps[[idx]]
    tK <- as.numeric(truth_main[[dgp$label]][[kappa_name]])
    for (m in methods) {
      s_mod <- summaries_main_mod[[idx]][[m]][[key]]
      s_big <- summaries_main_big[[idx]][[m]][[key]]
      lines <- c(
        lines,
        paste0(
          dgp$label, " & ",
          method_labels[[m]], " & ",
          fmt(tK, 4), " & ",
          fmt_bias_sd(s_mod[["bias"]], s_mod[["sd"]], 4), " & ",
          fmt_bias_sd(s_big[["bias"]], s_big[["sd"]], 4), "\\\\"
        )
      )
    }
    if (idx != length(dgps)) lines <- c(lines, "\\addlinespace")
  }

  footer <- c("\\bottomrule", "\\end{tabular}")
  c(header, lines, footer)
}

latex_note_dgps <- function() {
  lines <- c(
    "\\textit{Note.} Rating model: $T_i\\sim\\mathrm{Categorical}(p)$; conditional on $T_i$, rater $j$ reports $T_i$ with probability $\\rho_{ij}$ and otherwise guesses from $p$.",
    "Missingness indicators $M_{ij}\\in\\{0,1\\}$ are independent across $(i,j)$ given the specified mechanism.",
    paste0(
      "DGP~A: exchangeable raters ($\\rho_{ij}=0.92$) and MCAR with rater-specific $\\pi_j=\\Pr(M_{ij}=1)$ given by $(",
      paste(fmt(dgpA$pi_rater, 2), collapse = ", "),
      ")$."
    ),
    paste0(
      "DGP~B: non-exchangeable raters with $\\rho_{ij}=\\rho_j$ where $\\rho_j=(",
      paste(fmt(dgpB$rho_base, 2), collapse = ", "),
      ")$ and MCAR with $\\pi_j=(",
      paste(fmt(dgpB$pi_rater, 2), collapse = ", "),
      ")$ (higher-skill raters observed more often)."
    ),
    paste0(
      "DGP~C: category-dependent difficulty with $\\rho_{ij}=0.92\\,m(T_i)$ where $m=(",
      paste(fmt(dgpC$rho_truth_mult, 2), collapse = ", "),
      ")$, and MAR with $\\pi_c=\\Pr(M_{ij}=1\\mid T_i=c)=(",
      paste(fmt(dgpC$pi_truth, 2), collapse = ", "),
      ")$."
    ),
    "Main table uses nominal (identity) loss; Appendix table uses linear loss $\\ell(a,b)=|a-b|/(C-1)$ (treating categories as ordered)."
  )
  paste(lines, collapse = " ")
}

lines <- c(
  "% Auto-generated by code/analysis/simulations_raw_three_estimators.R",
  "\\documentclass[11pt]{article}",
  "\\usepackage{amsmath,amssymb,booktabs}",
  "\\usepackage[margin=1in]{geometry}",
  "\\begin{document}",
  "\\section*{Simulations (draft)}",
  "",
  "We simulate raw categorical ratings with $R=6$ raters and $C=5$ categories,",
  "using $T_i\\sim\\mathrm{Categorical}(p)$ with empirical $p$ taken from \\texttt{misskappa::dat.fleiss1971}.",
  paste0("Complete-data truth is approximated by a large Monte Carlo run ($n=", n_truth, "$) under each DGP."),
  "",
  paste0("Moderate sample size: $n=", n_mod, "$ with $B=", B_mod, "$ replications."),
  paste0("Large sample size: $n_{big}=", n_big, "$ with $B_{big}=", B_big, "$ replications."),
  "",
  "\\subsection*{Main: Conger's kappa (nominal/identity loss)}",
  "\\begin{table}[!ht]",
  "\\centering",
  "\\small",
  "\\setlength{\\tabcolsep}{4pt}",
  "\\caption{Simulation bias and Monte Carlo SD for Conger's $\\kappa_C$ (nominal/identity loss).}",
  latex_table_conger(),
  "\\vspace{3pt}",
  "\\begin{minipage}{0.95\\linewidth}",
  "\\footnotesize",
  latex_note_dgps(),
  "\\end{minipage}",
  "\\end{table}",
  "",
  "\\subsection*{Appendix (draft)}",
  "\\subsubsection*{Conger's kappa (linear/absolute loss)}",
  "\\begin{table}[!ht]",
  "\\centering",
  "\\small",
  "\\setlength{\\tabcolsep}{4pt}",
  "\\caption{Simulation bias and Monte Carlo SD for Conger's $\\kappa_C$ (linear/absolute loss).}",
  latex_table_conger_linear(),
  "\\end{table}",
  "",
  "\\subsubsection*{Fleiss' kappa (identity loss)}",
  "\\begin{center}\\small\\setlength{\\tabcolsep}{4pt}",
  latex_table_appendix_kappa("Fleiss", "fleiss"),
  "\\end{center}",
  "",
  "\\subsubsection*{Brennan--Prediger (identity loss)}",
  "\\begin{center}\\small\\setlength{\\tabcolsep}{4pt}",
  latex_table_appendix_kappa("Brennan-Prediger", "bp"),
  "\\end{center}",
  "",
  "\\end{document}"
)

writeLines(lines, out_tex)
cat("Wrote ", out_tex, "\n", sep = "")

snippet_main <- c(
  "% Auto-generated by code/analysis/simulations_raw_three_estimators.R",
  "% Snippet for inclusion in kappa-missing.lyx via \\input{}",
  "\\begin{table}[!ht]",
  "\\centering",
  "\\small",
  "\\setlength{\\tabcolsep}{4pt}",
  "\\caption{Simulation bias and Monte Carlo SD for Conger's $\\kappa_C$ (nominal/identity loss).}\\label{tab:sim-conger}",
  latex_table_conger(),
  "\\vspace{3pt}",
  "\\begin{minipage}{0.95\\linewidth}",
  "\\footnotesize",
  latex_note_dgps(),
  "\\end{minipage}",
  "\\end{table}"
)
writeLines(snippet_main, out_snippet_main)
cat("Wrote ", out_snippet_main, "\n", sep = "")

snippet_appendix <- c(
  "% Auto-generated by code/analysis/simulations_raw_three_estimators.R",
  "% Snippet for inclusion in an appendix via \\input{}",
  "\\subsubsection*{Conger's kappa (linear/absolute loss)}",
  "\\begin{table}[!ht]",
  "\\centering",
  "\\small",
  "\\setlength{\\tabcolsep}{4pt}",
  "\\caption{Simulation bias and Monte Carlo SD for Conger's $\\kappa_C$ (linear/absolute loss).}\\label{tab:sim-conger-linear}",
  latex_table_conger_linear(),
  "\\end{table}",
  "",
  "\\subsubsection*{Fleiss' kappa (identity loss)}",
  "\\begin{center}\\small\\setlength{\\tabcolsep}{4pt}",
  latex_table_appendix_kappa("Fleiss", "fleiss"),
  "\\end{center}",
  "",
  "\\subsubsection*{Brennan--Prediger (identity loss)}",
  "\\begin{center}\\small\\setlength{\\tabcolsep}{4pt}",
  latex_table_appendix_kappa("Brennan-Prediger", "bp"),
  "\\end{center}"
)
writeLines(snippet_appendix, out_snippet_appendix)
cat("Wrote ", out_snippet_appendix, "\n", sep = "")

if (Sys.getenv("SIM_COMPILE", "1") == "1") {
  cmd <- "pdflatex"
  args <- c(
    "-interaction=nonstopmode",
    "-halt-on-error",
    "-output-directory", "notes",
    out_tex
  )
  status <- system2(cmd, args = args)
  if (!identical(status, 0L)) stop("pdflatex failed (exit code ", status, ")")
  cat("Compiled notes/simulations-raw-three-estimators.pdf\n")
}
