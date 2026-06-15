# Probe 41 — verify the direct-on-counts exchangeable CS quadratic kappa against
# (a) the current table-expansion nt_fiml, (b) the moment count kappas
# (fleiss_cuzick, cat_fiml), and (c) calibrate its robust (sandwich) SEs.
#
# Questions:
#   Q1 EXACT.  Does cs_fiml (direct, no table) reproduce kappa_counts(nt_fiml)
#              (table expansion) on PARTIAL counts, point estimate AND SE?
#   Q2 COMPARE.How do listwise / pairwise / fiml relate to each other and to
#              fleiss_cuzick / cat_fiml on complete and partial counts?
#   Q3 CALIB.  Do the robust SEs of the three direct estimators cover at 95%
#              under complete / MCAR (all consistent) and MAR (only fiml)?

suppressMessages(library(misskappa))
source("experiments/workbench/41-counts-exchangeable-cs-direct/direct_cs.R")
set.seed(41)

VALUES <- 1:5
RT <- 5L

# ---- count-data simulator ---------------------------------------------------
# Latent CS-normal ratings rounded to 1..5, then amputated at the entry level
# (MCAR or value-based MAR), then tallied to a counts row per subject.
sim_raw <- function(n, mu = 3, psi = 0.8, e = 1.2) {
  b <- rnorm(n, 0, sqrt(psi))
  Y <- mu + outer(b, rep(1, RT)) + matrix(rnorm(n * RT, 0, sqrt(e)), n, RT)
  Y <- pmin(pmax(round(Y), 1L), 5L)
  matrix(as.integer(Y), n, RT)
}
amputate <- function(Y, mech, rate = 0.3) {
  n <- nrow(Y)
  if (mech == "complete") return(Y)
  if (mech == "mcar") {
    Y[matrix(runif(n * RT) < rate, n, RT)] <- NA
  } else if (mech == "mar") {
    # rater 1 always observed; a rater goes missing with prob increasing in the
    # observed value of rater 1 -> ignorable (depends on observed data only).
    pr <- plogis(-0.5 + 0.8 * scale(Y[, 1])[, 1])
    for (j in 2:RT) Y[runif(n) < pr, j] <- NA
  }
  Y
}
to_counts <- function(Y) {
  t(apply(Y, 1, function(r) tabulate(r[is.finite(r)], nbins = 5L)))
}

# ---- helper: pull (est, se) from any misskappa_estimate ---------------------
es <- function(f) c(F = unname(f$estimates[1]), seF = sqrt(diag(f$vcov))[1],
                    BP = unname(f$estimates[2]), seBP = sqrt(diag(f$vcov))[2])

fit_all <- function(counts, r_total = RT) {
  list(
    fc        = es(kappa_counts(counts, "fleiss_cuzick", "quadratic", VALUES)),
    cat_fiml  = es(kappa_counts(counts, "cat_fiml", "quadratic", VALUES, r_total = r_total)),
    tbl_fiml  = es(kappa_counts(counts, "nt_fiml", "quadratic", VALUES, r_total = r_total)),
    d_fiml    = es(cs_fiml(counts, VALUES, r_total)),
    d_pair    = es(cs_pairwise(counts, VALUES, r_total)),
    d_list    = es(cs_listwise(counts, VALUES, r_total))
  )
}

cat("=====================================================================\n")
cat("Q1/Q2  Cross-estimator comparison\n")
cat("=====================================================================\n\n")

cat("--- (a) dat.fleiss1971 (complete, 30x5, r_total=6) ---\n")
cmp <- fit_all(dat.fleiss1971, r_total = 6L)
print(round(do.call(rbind, cmp), 6))

cat("\n--- (b) simulated PARTIAL counts (n=400, MCAR 30%) ---\n")
Yp <- amputate(sim_raw(400), "mcar", 0.3)
cp <- to_counts(Yp)
cat("row-sum distribution:", paste0(names(table(rowSums(cp))), ":", table(rowSums(cp)),
    collapse = "  "), "\n")
cmpp <- fit_all(cp)
print(round(do.call(rbind, cmpp), 6))

d_exact <- max(abs(cmpp$d_fiml - cmpp$tbl_fiml))
cat(sprintf("\nQ1  max|direct_fiml - table_fiml| (est+se) = %.3e\n", d_exact))

# ---- larger exact-equivalence sweep ----------------------------------------
cat("\n--- Q1 exact-equivalence sweep over 30 random partial-count datasets ---\n")
worst <- 0
for (rep in 1:30) {
  Y <- amputate(sim_raw(sample(150:500, 1)),
                sample(c("mcar", "mar"), 1), runif(1, 0.15, 0.45))
  cc <- to_counts(Y)
  cc <- cc[rowSums(cc) > 0, , drop = FALSE]
  a <- es(cs_fiml(cc, VALUES, RT))
  b <- es(kappa_counts(cc, "nt_fiml", "quadratic", VALUES, r_total = RT))
  worst <- max(worst, max(abs(a - b)))
}
cat(sprintf("worst |direct - table| (est+se) across 30 datasets = %.3e\n", worst))

# =============================================================================
# Q3  Robust-SE calibration
# =============================================================================
cat("\n=====================================================================\n")
cat("Q3  Robust SE calibration  (B reps)\n")
cat("=====================================================================\n")

B <- 600; n <- 200
# truth = CS-FIML estimand on a huge complete rounded sample
big <- to_counts(sim_raw(200000))
truth <- cs_fiml(big, VALUES, RT)$estimates
cat(sprintf("estimand (huge-n complete): Fleiss=%.4f  BP=%.4f\n\n",
            truth[1], truth[2]))

methods <- c("d_fiml", "d_pair", "d_list", "fc")
fit_one <- function(counts) {
  tryCatch(list(
    d_fiml = es(cs_fiml(counts, VALUES, RT)),
    d_pair = es(cs_pairwise(counts, VALUES, RT)),
    d_list = es(cs_listwise(counts, VALUES, RT)),
    fc     = es(kappa_counts(counts, "fleiss_cuzick", "quadratic", VALUES))
  ), error = function(e) NULL)
}

rows <- list()
for (mech in c("complete", "mcar", "mar")) {
  acc <- lapply(methods, function(m) matrix(NA_real_, B, 4))
  names(acc) <- methods
  for (b in seq_len(B)) {
    cc <- to_counts(amputate(sim_raw(n), mech, 0.3))
    cc <- cc[rowSums(cc) > 0, , drop = FALSE]
    r <- fit_one(cc)
    if (is.null(r)) next
    for (m in methods) acc[[m]][b, ] <- r[[m]]
  }
  for (m in methods) {
    A <- acc[[m]]; ok <- is.finite(A[, 1])
    for (coef_i in 1:2) {
      estcol <- if (coef_i == 1) 1 else 3
      secol  <- if (coef_i == 1) 2 else 4
      est <- A[ok, estcol]; se <- A[ok, secol]; tr <- truth[coef_i]
      cov <- mean(abs(est - tr) <= qnorm(0.975) * se)
      rows[[length(rows) + 1]] <- data.frame(
        mech = mech, method = m,
        coef = c("Fleiss", "BP")[coef_i], truth = tr,
        bias = mean(est) - tr, sd_emp = sd(est),
        mean_se = mean(se), se_ratio = mean(se) / sd(est),
        coverage = cov)
    }
  }
}
out <- do.call(rbind, rows)
op <- out; num <- sapply(op, is.numeric); op[num] <- lapply(op[num], round, 4)
print(op, row.names = FALSE)
saveRDS(out, "experiments/workbench/41-counts-exchangeable-cs-direct/calibration.rds")
cat("\ndone.\n")
