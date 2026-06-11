#!/usr/bin/env Rscript
#
# 28-concordance-landscape
#
# Software + literature comparison: where misskappa sits relative to the
# concordance-correlation-coefficient (CCC) world, benchmarked numerically.
#
# Produces four result tables:
#   equivalence.csv      - misskappa quadratic Conger kappa == Lin's CCC, checked
#                          against several independent CCC implementations.
#   bias_correction.csv  - the "unbiased denominator" question: the moment CCC's
#                          squared mean-difference term is biased; Carrasco's
#                          variance-components estimator corrects it. We reproduce
#                          the correction by hand and show its size vs n / shift.
#   repeated.csv         - repeated-measures CCC on bpres via the three CCC
#                          paradigms (variance components, U-statistics, SimplyAgree),
#                          with misskappa's (in)applicability noted.
#   categorical.csv      - misskappa's categorical / vector-valued path on the
#                          Vanbelle CRACKLES data vs irr / irrCAC.
# plus inventory.csv (the software-landscape table) and metadata.csv.
#
# Every comparator is tryCatch-guarded: a missing or failing package degrades to
# an NA row with a recorded note rather than killing the run.

script_arg <- commandArgs(FALSE)
script_file <- sub("^--file=", "", script_arg[grep("^--file=", script_arg)][1L])
script_dir <- if (length(script_file) && !is.na(script_file)) {
  dirname(normalizePath(script_file))
} else getwd()
results_dir <- file.path(script_dir, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

suppressMessages(library(misskappa))
have <- function(p) requireNamespace(p, quietly = TRUE)

# Load a dataset object from a package that may ship it as LazyData (no data/
# dir), e.g. cccrm::bpres. Tries utils::data(), then the namespace lazy-load db.
get_data <- function(name, pkg) {
  e <- new.env()
  suppressWarnings(try(utils::data(list = name, package = pkg, envir = e),
                       silent = TRUE))
  if (exists(name, envir = e, inherits = FALSE)) return(get(name, envir = e))
  ns <- asNamespace(pkg)
  if (exists(name, envir = ns, inherits = FALSE)) return(get(name, envir = ns))
  # last resort: attach the package so search-path lazy data becomes visible
  suppressMessages(suppressWarnings(require(pkg, character.only = TRUE)))
  if (exists(name)) return(get(name))
  stop("cannot load dataset '", name, "' from package '", pkg, "'")
}
note_log <- character(0)
log_note <- function(...) {
  msg <- paste0(...)
  note_log <<- c(note_log, msg)
  message("[note] ", msg)
}

## ---- CCC helpers (explicit formulas, so the conventions are visible) --------

# Lin's concordance correlation coefficient, plug-in (1/n) moments.
ccc_moment <- function(x, y) {
  n <- length(x); mx <- mean(x); my <- mean(y)
  sxx <- sum((x - mx)^2) / n
  syy <- sum((y - my)^2) / n
  sxy <- sum((x - mx) * (y - my)) / n
  2 * sxy / (sxx + syy + (mx - my)^2)
}

# Location-"corrected" CCC: (xbar - ybar)^2 overestimates (mu_x - mu_y)^2 by
# Var(xbar - ybar) = (s_xx + s_yy - 2 s_xy) / n. Subtract it (truncate at 0).
# NB: this is NOT a proper bias correction -- it patches only the denominator's
# location term and ignores the dominant O(1/n) bias, which is the ratio
# nonlinearity E[N/D] != E[N]/E[D]. Kept to show how little it does vs jackknife.
ccc_loc_corrected <- function(x, y) {
  n <- length(x); mx <- mean(x); my <- mean(y)
  sxx <- sum((x - mx)^2) / n
  syy <- sum((y - my)^2) / n
  sxy <- sum((x - mx) * (y - my)) / n
  var_md <- (sxx + syy - 2 * sxy) / n
  2 * sxy / (sxx + syy + max((mx - my)^2 - var_md, 0))
}

# Jackknife bias-corrected CCC: a *proper* second-order O(1/n) bias remover
# (Quenouille). Estimates the same bias the second-order influence function /
# von Mises expansion targets, with no closed form needed.
ccc_jackknife <- function(x, y) {
  n <- length(x); th <- ccc_moment(x, y)
  thi <- vapply(seq_len(n), function(i) ccc_moment(x[-i], y[-i]), numeric(1))
  n * th - (n - 1) * mean(thi)
}

# misskappa: quadratic-weighted Conger kappa on a subjects-by-2 score matrix.
ccc_misskappa <- function(x, y) {
  X <- cbind(x, y); storage.mode(X) <- "double"
  unname(coef(kappa(X, estimator = "pairwise"))["Conger"])
}

## ============================================================================
## Table 1 + bias correction: a non-repeated continuous paired dataset.
## bpres (Carrasco) restricted to the first replicate gives one systolic-BP
## reading per method per subject -- a clean paired continuous comparison.
## ============================================================================

equiv_rows <- list()
bias_rows  <- list()

bpres <- if (have("cccrm")) get_data("bpres", "cccrm") else NULL

if (have("cccrm")) {
  b1 <- bpres[bpres$NM == 1, c("ID", "METODE", "SIS")]
  w <- reshape(b1, idvar = "ID", timevar = "METODE", direction = "wide")
  w <- w[complete.cases(w), ]
  x <- w$SIS.1; y <- w$SIS.2; n <- nrow(w)
  meandiff <- mean(x) - mean(y)

  add_equiv <- function(impl, fun, detail) {
    est <- tryCatch(fun(), error = function(e) {
      log_note(impl, " failed: ", conditionMessage(e)); NA_real_
    })
    equiv_rows[[length(equiv_rows) + 1L]] <<- data.frame(
      dataset = "bpres (systolic, 1st replicate)", n = n,
      implementation = impl, ccc = est, detail = detail,
      stringsAsFactors = FALSE)
  }

  add_equiv("misskappa (quadratic Conger kappa)", function() ccc_misskappa(x, y),
            "pairwise-available moment estimator")
  add_equiv("moment plug-in (1/n), by hand", function() ccc_moment(x, y),
            "Lin (1989) definition")
  if (have("DescTools"))
    add_equiv("DescTools::CCC", function() DescTools::CCC(x, y)$rho.c$est,
              "Lin's CCC, z-transform CI")
  else log_note("DescTools not available")
  add_equiv("cccrm::cccUst (U-statistic)",
            function() unname(cccrm::cccUst(b1, ry = "SIS", rmet = "METODE")[1]),
            "King-Chinchilli-Carrasco U-statistic")
  add_equiv("cccrm::ccc_vc (variance components)",
            function() unname(cccrm::ccc_vc(b1, ry = "SIS", rind = "ID",
                                            rmet = "METODE")$ccc[1]),
            "Carrasco & Jover REML -- BIAS-CORRECTED")

  # --- bias-correction anchor on bpres ---
  vc_anchor <- tryCatch(unname(cccrm::ccc_vc(b1, ry = "SIS", rind = "ID",
                                             rmet = "METODE")$ccc[1]),
                        error = function(e) NA_real_)
  ust_anchor <- tryCatch(unname(cccrm::cccUst(b1, ry = "SIS", rmet = "METODE")[1]),
                         error = function(e) NA_real_)
  bias_rows[[length(bias_rows) + 1L]] <- data.frame(
    setting = "bpres systolic, 1st replicate", n = n, mean_diff = meandiff,
    ccc_moment = ccc_moment(x, y),
    ccc_loc_corrected = ccc_loc_corrected(x, y),
    ccc_vc = vc_anchor, ccc_uststat = ust_anchor,
    stringsAsFactors = FALSE)
} else {
  log_note("cccrm not available: Table 1 / repeated-measures comparators skipped")
}

# A second equivalence anchor on irr::anxiety (ties to the repo's other examples).
if (have("irr")) {
  data("anxiety", package = "irr")
  ax <- as.matrix(anxiety[, 1:2]); storage.mode(ax) <- "double"
  equiv_rows[[length(equiv_rows) + 1L]] <- data.frame(
    dataset = "irr::anxiety (raters 1-2)", n = nrow(ax),
    implementation = "misskappa (quadratic Conger kappa)",
    ccc = ccc_misskappa(ax[, 1], ax[, 2]),
    detail = "equal rater means -> location term vanishes",
    stringsAsFactors = FALSE)
  equiv_rows[[length(equiv_rows) + 1L]] <- data.frame(
    dataset = "irr::anxiety (raters 1-2)", n = nrow(ax),
    implementation = "moment plug-in (1/n), by hand",
    ccc = ccc_moment(ax[, 1], ax[, 2]),
    detail = "Lin (1989) definition",
    stringsAsFactors = FALSE)
  if (have("DescTools"))
    equiv_rows[[length(equiv_rows) + 1L]] <- data.frame(
      dataset = "irr::anxiety (raters 1-2)", n = nrow(ax),
      implementation = "DescTools::CCC",
      ccc = DescTools::CCC(ax[, 1], ax[, 2])$rho.c$est,
      detail = "Lin's CCC", stringsAsFactors = FALSE)
}

# --- bias-correction as a function of (n, mean-shift): synthetic, averaged ---
# Single draws are seed-noisy, so average the moment and corrected CCCs over many
# replicates. The gap is the mean shift of the location correction, which is
# O(1/n) and largest when the mean difference is marginal relative to the noise.
set.seed(28000)
n_bias_rep <- 400L
for (nn in c(20L, 50L, 200L)) {
  for (shift in c(0.5, 2.0)) {
    mom <- numeric(n_bias_rep); loc <- numeric(n_bias_rep); md <- numeric(n_bias_rep)
    for (r in seq_len(n_bias_rep)) {
      z <- rnorm(nn); x <- 10 + 3 * z + rnorm(nn, 0, 1.2)
      y <- (10 + shift) + 3 * z + rnorm(nn, 0, 1.2)
      mom[r] <- ccc_moment(x, y); loc[r] <- ccc_loc_corrected(x, y)
      md[r] <- mean(x) - mean(y)
    }
    bias_rows[[length(bias_rows) + 1L]] <- data.frame(
      setting = sprintf("synthetic (n=%d, shift=%.1f, %d reps)", nn, shift, n_bias_rep),
      n = nn, mean_diff = mean(md),
      ccc_moment = mean(mom),
      ccc_loc_corrected = mean(loc),
      ccc_vc = NA_real_, ccc_uststat = NA_real_,
      stringsAsFactors = FALSE)
  }
}

equivalence <- do.call(rbind, equiv_rows)
bias_df <- do.call(rbind, bias_rows)
bias_df$gap_vc      <- bias_df$ccc_vc - bias_df$ccc_moment
bias_df$gap_loc     <- bias_df$ccc_loc_corrected - bias_df$ccc_moment
write.csv(equivalence, file.path(results_dir, "equivalence.csv"), row.names = FALSE)
write.csv(bias_df, file.path(results_dir, "bias_correction.csv"), row.names = FALSE)

## ---- proper bias correction vs the denominator patch -----------------------
## Controlled DGP with a KNOWN true CCC, so we can measure how much of the
## moment estimator's bias each correction actually removes. The denominator
## patch removes a sliver; the jackknife (a proper second-order corrector)
## removes nearly all of it -- the point that the "unbiased denominator" is not
## a real bias correction.
## DGP: x = lam*z + e_x, y = bias_sd*totsd + lam*z + e_y; corr(x,y) = rho;
## true CCC = 2*rho / (2 + bias_sd^2).
gen_pair <- function(n, bias_sd, rho) {
  lam <- 1; se <- sqrt(lam^2 * (1 - rho) / rho); totsd <- sqrt(lam^2 + se^2)
  z <- rnorm(n)
  x <- lam * z + rnorm(n, 0, se)
  y <- bias_sd * totsd + lam * z + rnorm(n, 0, se)
  cbind(x, y)
}
true_ccc <- function(bias_sd, rho) 2 * rho / (2 + bias_sd^2)
set.seed(28100)
n_meth_rep <- 3000L
meth_grid <- expand.grid(n = c(8L, 15L, 30L, 100L), bias_sd = c(0, 1), rho = 0.5)
meth_rows <- lapply(seq_len(nrow(meth_grid)), function(i) {
  nn <- meth_grid$n[i]; b <- meth_grid$bias_sd[i]; r <- meth_grid$rho[i]
  tru <- true_ccc(b, r)
  M <- L <- J <- numeric(n_meth_rep)
  for (k in seq_len(n_meth_rep)) {
    X <- gen_pair(nn, b, r)
    M[k] <- ccc_moment(X[, 1], X[, 2])
    L[k] <- ccc_loc_corrected(X[, 1], X[, 2])
    J[k] <- ccc_jackknife(X[, 1], X[, 2])
  }
  bias_m <- mean(M) - tru
  pct <- function(corr) if (abs(bias_m) < 1e-8) NA_real_ else
    100 * (mean(corr) - mean(M)) / (tru - mean(M))
  data.frame(
    n = nn, bias_sd = b, true_ccc = tru,
    bias_moment = bias_m,
    bias_loc_patch = mean(L) - tru,
    bias_jackknife = mean(J) - tru,
    pct_removed_loc_patch = pct(L),
    pct_removed_jackknife = pct(J),
    stringsAsFactors = FALSE)
})
bias_methods <- do.call(rbind, meth_rows)
write.csv(bias_methods, file.path(results_dir, "bias_methods.csv"), row.names = FALSE)

## ============================================================================
## Table 2: repeated-measures CCC on the full bpres (2 methods x 2 replicates).
## The three paradigms the CCC literature uses for repeated measurement.
## ============================================================================

rep_rows <- list()
add_rep <- function(impl, paradigm, fun, applies = TRUE) {
  res <- tryCatch(fun(), error = function(e) {
    log_note(impl, " failed: ", conditionMessage(e)); c(NA, NA, NA, NA)
  })
  rep_rows[[length(rep_rows) + 1L]] <<- data.frame(
    implementation = impl, paradigm = paradigm,
    ccc = res[1], se = res[2], lcl = res[3], ucl = res[4],
    stringsAsFactors = FALSE)
}

if (have("cccrm")) {
  add_rep("cccrm::ccc_vc", "variance components (REML)", function() {
    v <- cccrm::ccc_vc(bpres, ry = "SIS", rind = "ID", rmet = "METODE",
                       rtime = "NM")$ccc
    c(v["CCC"], v["SE CCC"], v["LL CI 95%"], v["UL CI 95%"])
  })
  add_rep("cccrm::cccUst", "U-statistics", function() {
    u <- cccrm::cccUst(bpres, ry = "SIS", rmet = "METODE", rtime = "NM")
    c(u["CCC"], u["SE CCC"], u["LL CI 95%"], u["UL CI 95%"])
  })
}
if (have("SimplyAgree")) {
  add_rep("SimplyAgree::agree_reps (ccc)", "U-statistics (K-C-C)", function() {
    wr <- reshape(bpres[, c("ID", "NM", "METODE", "SIS")],
                  idvar = c("ID", "NM"), timevar = "METODE", direction = "wide")
    names(wr)[names(wr) == "SIS.1"] <- "m1"
    names(wr)[names(wr) == "SIS.2"] <- "m2"
    wr <- wr[complete.cases(wr[, c("m1", "m2")]), ]
    ar <- suppressWarnings(SimplyAgree::agree_reps(
      x = "m1", y = "m2", id = "ID", data = wr, delta = 5, ccc = TRUE))
    cc <- ar$ccc.xy
    c(cc[[1]], NA, cc[[2]], cc[[3]])
  })
} else log_note("SimplyAgree not available")

# misskappa has no mixed-model CCC, but its vector path DOES handle the replicate
# structure: treat the 2 methods as raters and the 2 replicates as vector
# components. With equal (identity) weights this is the U-statistic route, so it
# should track cccUst. Record it as a real value, not "none".
if (have("cccrm")) {
  add_rep("misskappa (vector path)", "U-statistics (replicates as vector)",
          function() {
    ids <- sort(unique(bpres$ID))
    Av <- array(NA_real_, c(length(ids), 2L, 2L))
    for (ii in seq_along(ids)) {
      sub <- bpres[bpres$ID == ids[ii], ]
      for (m in 1:2) for (k in 1:2) {
        v <- sub$SIS[sub$METODE == m & sub$NM == k]
        if (length(v) == 1L) Av[ii, m, k] <- v
      }
    }
    keep <- apply(is.finite(Av), 1L, all)
    Av <- Av[keep, , , drop = FALSE]
    fit <- kappa(Av, estimator = "pairwise", weight = "quadratic")
    se <- sqrt(vcov(fit)["Conger", "Conger"])
    est <- unname(coef(fit)["Conger"])
    c(est, se, est - 1.96 * se, est + 1.96 * se)
  })
}

repeated <- do.call(rbind, rep_rows)
write.csv(repeated, file.path(results_dir, "repeated.csv"), row.names = FALSE)

## ============================================================================
## U-statistic CCC (cccrm::cccUst): the closest prior art to misskappa's
## quadratic / vector path. Two demonstrations written to results/:
##   delta_dial.csv        - one estimator spans kappa (delta=0) to CCC (delta=1)
##   repeated_identity.csv - cccUst(Dmat) == misskappa vector kappa (W = Dmat),
##                           exactly on the diagonal, diverging off-diagonal
##                           (cccUst's absolute vs misskappa's signed cross-terms).
## ============================================================================

dial_rows <- list(); ident_rows <- list()
if (have("cccrm")) {
  ## --- delta dial on a 2-method ordinal data set ---
  set.seed(28200)
  nd <- 150L
  zc <- sample(1:4, nd, TRUE, c(.4, .3, .2, .1))
  jit <- function() pmin(pmax(zc + sample(c(0, 0, 1, -1), nd, TRUE), 1L), 4L)
  d1 <- jit(); d2 <- jit()
  longd <- data.frame(id = rep(seq_len(nd), 2L), met = rep(1:2, each = nd),
                      y = c(d1, d2))
  Wd <- cbind(d1, d2); storage.mode(Wd) <- "double"
  dial_rows[[1]] <- data.frame(
    delta = 0, distance = "0/1 (nominal)", target = "unweighted kappa",
    cccUst = unname(cccrm::cccUst(longd, "y", "met", delta = 0)[1]),
    misskappa = unname(coef(kappa(Wd, estimator = "ipw", weight = "nominal"))["Conger"]),
    cross_check = if (have("irr")) suppressWarnings(irr::kappa2(Wd)$value) else NA_real_,
    cross_name = "irr::kappa2", stringsAsFactors = FALSE)
  dial_rows[[2]] <- data.frame(
    delta = 1, distance = "squared", target = "Lin CCC = quadratic kappa",
    cccUst = unname(cccrm::cccUst(longd, "y", "met", delta = 1)[1]),
    misskappa = unname(coef(kappa(Wd, estimator = "pairwise"))["Conger"]),
    cross_check = if (have("DescTools")) DescTools::CCC(d1, d2)$rho.c$est else NA_real_,
    cross_name = "DescTools::CCC", stringsAsFactors = FALSE)

  ## --- cccUst(Dmat) == misskappa vector (W); diagonal matches, off-diag diverges ---
  ## Components 1 and 2 are given anti-correlated method differences so the off-
  ## diagonal cross-term (signed in misskappa, absolute in cccUst) actually bites.
  set.seed(28300)
  nv <- 300L; Tt <- 3L
  tru <- matrix(rnorm(nv * Tt), nv, Tt) %*% diag(c(2, 1.5, 1.8))
  e <- rnorm(nv)                                  # latent that anti-couples comp 1 vs 2
  Mm1 <- tru + matrix(rnorm(nv * Tt, 0, 0.8), nv, Tt)
  Mm2 <- tru + matrix(rnorm(nv * Tt, 0, 0.8), nv, Tt) +
    cbind(0.9 * e, -0.9 * e, rep(0, nv)) +        # opposite-sign bias on comp 1 vs 2
    matrix(rep(c(0.3, -0.2, 0.1), each = nv), nv, Tt)
  longv <- data.frame(SUBJ = rep(seq_len(nv), 2L * Tt),
                      VNUM = rep(rep(seq_len(Tt), each = nv), 2L),
                      MET = rep(1:2, each = nv * Tt),
                      Y = c(as.vector(Mm1), as.vector(Mm2)))
  Av <- array(NA_real_, c(nv, 2L, Tt)); Av[, 1, ] <- Mm1; Av[, 2, ] <- Mm2
  kvq <- getFromNamespace("kappa_vector_quadratic", "misskappa")
  add_ident <- function(label, Dmat) {
    u <- tryCatch(unname(cccrm::cccUst(longv, "Y", "MET", rtime = "VNUM", Dmat = Dmat)[1]),
                  error = function(e) NA_real_)
    v <- tryCatch(unname(kvq(Av, method = "pairwise", W = Dmat)$estimates["Conger"]),
                  error = function(e) NA_real_)
    ident_rows[[length(ident_rows) + 1L]] <<- data.frame(
      weights = label, cccUst = u, misskappa_vector = v, diff = u - v,
      stringsAsFactors = FALSE)
  }
  wdiag <- c(3, 1, 2)
  add_ident("diagonal weights (W diagonal)", diag(wdiag))
  Woff <- diag(wdiag); Woff[1, 2] <- Woff[2, 1] <- 1.5   # PSD: det of (1,2) block = 0.75
  add_ident("off-diagonal weights (W[1,2]=1.5)", Woff)

  write.csv(do.call(rbind, dial_rows), file.path(results_dir, "delta_dial.csv"),
            row.names = FALSE)
  write.csv(do.call(rbind, ident_rows), file.path(results_dir, "repeated_identity.csv"),
            row.names = FALSE)
} else {
  log_note("cccrm not available: delta-dial and U-statistic identity skipped")
}

## ============================================================================
## Table 3: categorical / vector-valued agreement on Vanbelle CRACKLES.
## Validate misskappa's scalar engine against irr / irrCAC at a single site,
## then show the vector-valued aggregate that only misskappa produces.
## ============================================================================

cat_rows <- list()
add_cat <- function(scope, impl, coef_name, fun) {
  est <- tryCatch(fun(), error = function(e) {
    log_note(impl, " failed: ", conditionMessage(e)); NA_real_
  })
  cat_rows[[length(cat_rows) + 1L]] <<- data.frame(
    scope = scope, implementation = impl, coefficient = coef_name,
    estimate = est, stringsAsFactors = FALSE)
}

V <- dat.vanbelle2019            # 20 x 28 x 6 binary array (lazy-loaded)
exp_cols <- paste0("EXP", 1:4)
site1 <- V[, exp_cols, 1]        # 20 x 4 binary matrix, upper-posterior left
storage.mode(site1) <- "integer"

# --- single-site scalar validation: misskappa vs irr vs irrCAC ---
site1_df <- as.data.frame(site1)
add_cat("EXP group, site U1 (scalar)", "misskappa (ipw, nominal)", "Fleiss",
        function() unname(coef(kappa(site1, estimator = "ipw",
                                     weight = "nominal"))["Fleiss"]))
add_cat("EXP group, site U1 (scalar)", "misskappa (cat_fiml, nominal)", "Fleiss",
        function() unname(coef(kappa(site1, estimator = "cat_fiml",
                                     weight = "nominal"))["Fleiss"]))
if (have("irr"))
  add_cat("EXP group, site U1 (scalar)", "irr::kappam.fleiss", "Fleiss",
          function() irr::kappam.fleiss(site1_df)$value)
if (have("irrCAC"))
  add_cat("EXP group, site U1 (scalar)", "irrCAC::fleiss.kappa.raw", "Fleiss",
          function() irrCAC::fleiss.kappa.raw(site1_df)$est$coeff.val)
# Conger's kappa at the same site: misskappa vs irrCAC.
add_cat("EXP group, site U1 (scalar)", "misskappa (ipw, nominal)", "Conger",
        function() unname(coef(kappa(site1, estimator = "ipw",
                                     weight = "nominal"))["Conger"]))
if (have("irrCAC"))
  add_cat("EXP group, site U1 (scalar)", "irrCAC::conger.kappa.raw", "Conger",
          function() irrCAC::conger.kappa.raw(site1_df)$est$coeff.val)

# --- vector-valued aggregate over all six sites: only misskappa ---
exp_array <- V[, exp_cols, ]     # 20 x 4 x 6
add_cat("EXP group, all 6 sites (vector)", "misskappa vector (Hamming)", "Conger",
        function() unname(coef(kappa(exp_array, estimator = "pairwise"))["Conger"]))
add_cat("EXP group, all 6 sites (vector)", "misskappa vector (Hamming)", "Fleiss",
        function() unname(coef(kappa(exp_array, estimator = "pairwise"))["Fleiss"]))
add_cat("EXP group, all 6 sites (vector)", "irr / irrCAC / DescTools", "Conger",
        function() stop("no vector-valued path"))

categorical <- do.call(rbind, cat_rows)
write.csv(categorical, file.path(results_dir, "categorical.csv"), row.names = FALSE)

## ============================================================================
## Vanbelle's multilevel (cluster-collapse) kappa vs misskappa's vector kappa.
## Both face "sites within patient". Vanbelle pools the sites under a within-
## cluster homogeneity assumption (one scalar); misskappa keeps them as a vector.
## They coincide iff homogeneity holds; otherwise they differ by the between-site
## covariance of the raters' marginals in the chance term (misskappa corrects
## chance per-site, Vanbelle pools the marginals first).
## ============================================================================

vec_conger <- function(A) unname(coef(kappa(A, estimator = "pairwise"))["Conger"])
pooled_conger <- function(A) {            # stack sites as rows = Vanbelle's point estimate
  d <- dim(A); M <- matrix(aperm(A, c(1, 3, 2)), nrow = d[1] * d[3], ncol = d[2])
  storage.mode(M) <- "integer"
  unname(coef(kappa(M, estimator = "ipw", weight = "nominal"))["Conger"])
}

# --- real CRACKLES: 7 observer groups vs Vanbelle (2019) Table 3 "All" Conger ---
vb_groups <- c("EXP", "NOR", "RUS", "WAL", "NLD", "PUL", "STU")
vb_published <- c(EXP = 0.56, NOR = 0.58, RUS = 0.20, WAL = 0.53,
                  NLD = 0.49, PUL = 0.40, STU = 0.37)   # PUL = her "PLN" row
vanbelle_fork <- do.call(rbind, lapply(vb_groups, function(g) {
  A <- V[, paste0(g, 1:4), ]; storage.mode(A) <- "double"
  data.frame(group = g,
             misskappa_vector = vec_conger(A),
             vanbelle_pooled = pooled_conger(A),
             vanbelle_published = unname(vb_published[g]),
             stringsAsFactors = FALSE)
}))
write.csv(vanbelle_fork, file.path(results_dir, "vanbelle_fork.csv"), row.names = FALSE)

# --- homogeneity is the fork: synthetic sites, identical vs differing marginals ---
sim_sites <- function(n, site_p, skill, seed) {
  set.seed(seed); S <- length(site_p); R <- length(skill)
  A <- array(NA_real_, c(n, R, S))
  for (s in seq_len(S)) {
    truth <- stats::rbinom(n, 1, site_p[s])
    for (r in seq_len(R)) A[, r, s] <- abs(truth - stats::rbinom(n, 1, 1 - skill[r]))
  }
  A
}
vb_skill <- c(0.90, 0.85, 0.80, 0.78)
Ah <- sim_sites(120L, rep(0.30, 6), vb_skill, 28401L)                 # homogeneous
Ax <- sim_sites(120L, c(.05, .30, .05, .30, .50, .50), vb_skill, 28402L)  # heterogeneous
vanbelle_homog <- data.frame(
  scenario = c("homogeneous (all 6 sites p=0.30)", "heterogeneous (sites differ)"),
  misskappa_vector = c(vec_conger(Ah), vec_conger(Ax)),
  vanbelle_pooled  = c(pooled_conger(Ah), pooled_conger(Ax)),
  stringsAsFactors = FALSE)
vanbelle_homog$diff <- vanbelle_homog$misskappa_vector - vanbelle_homog$vanbelle_pooled
write.csv(vanbelle_homog, file.path(results_dir, "vanbelle_homog.csv"), row.names = FALSE)

## ============================================================================
## Software-landscape inventory (authored here so the report renders one table).
## ============================================================================

inventory <- read.csv(text = '
package,language,computes,repeated_measures,missing_data,inference,reference
cccrm,R,"Lin CCC (generalized, weighted)",yes (longitudinal + replicates),listwise only,"VC delta + U-stat + bootstrap","Carrasco & Jover 2003; King-Chinchilli-Carrasco 2007"
SimplyAgree,R,"Lin CCC + Bland-Altman LoA + reliability",yes (nested/replicate),listwise only,"U-stat / mixed-model + bootstrap","Caldwell 2022 (JOSS)"
epiR (epi.ccc),R,"Lin CCC (paired)",no,listwise only,"z-transform","Lin 1989"
DescTools (CCC),R,"Lin CCC (paired)",no,na.rm only,"z-transform","Lin 1989"
agRee,R,"agreement indices incl. CCC",partial,limited,"bootstrap / Bayesian","Feng et al."
MethComp,R,"method comparison (Deming, LoA, ICC)",yes,listwise only,"mixed model / MCMC","Carstensen et al."
blandr,R,"Bland-Altman limits of agreement",no,listwise only,"normal-theory","Datta 2017"
multiagree,R,"multilevel multirater weighted kappa",yes (multilevel clustering),listwise/aligned,"multilevel Hotelling T2 / delta","Vanbelle 2017/2019"
irr,R,"kappa / ICC / agreement (categorical + ICC)",via ICC only,listwise only,"asymptotic","Gamer et al."
irrCAC,R,"chance-corrected agreement (AC1/AC2, Fleiss, alpha)",no,partial (weights),"asymptotic","Gwet 2014"
SAS %ccc / PROC,SAS,"Lin CCC (VC + U-stat)",yes,listwise,"VC delta / U-stat","Carrasco et al. 2013"
Stata (concord),Stata,"Lin CCC (paired)",no,listwise,"z-transform","Lin 1989"
pingouin,Python,"Lin CCC (paired)",no,listwise,"z-transform","Lin 1989"
misskappa,R,"Conger/Fleiss/BP kappa + alpha; quadratic kappa == Lin CCC",via vector path (= U-stat CCC),"IPW (MCAR) + FIML (MAR)","influence-function sandwich + Wald tests","van Oest & Moss 2026; Moss 2024"
', stringsAsFactors = FALSE, check.names = FALSE)
write.csv(inventory, file.path(results_dir, "inventory.csv"), row.names = FALSE)

## ---- metadata + run log -----------------------------------------------------

pkg_ver <- function(p) if (have(p)) as.character(packageVersion(p)) else "absent"
metadata <- data.frame(
  key = c("R_version", "misskappa", "cccrm", "SimplyAgree", "DescTools",
          "irrCAC", "irr", "multiagree", "epiR"),
  value = c(paste(R.version$major, R.version$minor, sep = "."),
            pkg_ver("misskappa"), pkg_ver("cccrm"), pkg_ver("SimplyAgree"),
            pkg_ver("DescTools"), pkg_ver("irrCAC"), pkg_ver("irr"),
            pkg_ver("multiagree"), pkg_ver("epiR")),
  stringsAsFactors = FALSE)
write.csv(metadata, file.path(results_dir, "metadata.csv"), row.names = FALSE)
writeLines(if (length(note_log)) note_log else "all comparators ran",
           file.path(results_dir, "run_notes.txt"))

cat("\n==== equivalence ====\n");  print(equivalence, row.names = FALSE)
cat("\n==== bias_correction ====\n"); print(bias_df, row.names = FALSE)
cat("\n==== repeated ====\n");     print(repeated, row.names = FALSE)
cat("\n==== categorical ====\n");  print(categorical, row.names = FALSE)
cat("\nNotes:\n"); cat(paste0(" - ", note_log, collapse = "\n"), "\n")
cat("\nWrote results to ", results_dir, "\n")
