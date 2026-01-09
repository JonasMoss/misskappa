# misskappa package audit

This is a living engineering note for the `misskappa` R package (`code/misskappa/` in this repo). It documents what exists today and what to improve next, without committing to removing major subsystems (EM / “quadratic”) yet.

## What the package currently exposes

Public exported API (`code/misskappa/NAMESPACE:3`):

- `kappa_raw()` for categorical ratings matrices
- `kappa_continuous()` for continuous ratings matrices
- `kappa_counts()` for count tables

All three wrappers dispatch into a single Rcpp layer (`code/misskappa/src/rcpp_interface.cpp:29`) with a `method` argument.

### Supported methods (as implemented)

- Raw categorical (`code/misskappa/R/kappa.R:34` + `code/misskappa/src/rcpp_interface.cpp:137`)
  - `method = "available"` → `misskappa::kappanp::kappa(...)`
  - `method = "ipw"` → `misskappa::kappanp::kappa(..., use_ipw=true)`
  - `method = "gwet"` → `misskappa::kappanp::kappa(..., use_gwet=true)`
  - `method = "ml"` → `emdiscrete::{preprocess_raw, run_em}` + `misskappa::kappaml::kappa(...)`
  - `method = "quadratic"` → `misskappa::kappaqp::kappa(...)` (special path)

- Continuous (`code/misskappa/R/kappa.R:100` + `code/misskappa/src/rcpp_interface.cpp:32`)
  - `method = "available"` / `"ipw"` / `"gwet"` → `misskappa::kappanp::kappa_continuous(..., use_ipw=..., use_gwet=...)`
  - `method = "quadratic"` → `misskappa::kappaqp::kappa(...)` (special path)

- Counts (`code/misskappa/R/kappa.R:152` + `code/misskappa/src/rcpp_interface.cpp:193`)
  - `method = "available"` → `misskappa::kappanp::kappa_counts(...)`
  - `method = "ml"` → `emdiscrete::{preprocess_counts, run_em}` + `misskappa::kappaml::kappa_counts(...)`
  - `method = "quadratic"` → `misskappa::kappaqp::kappa_counts(...)` (special path)

## C++ subsystem map

The C++ code is split by estimator family:

- `kappanp` (nonparametric / “MCAR-like” estimators): `code/misskappa/src/kappanp.cpp:1`
  - Implements available-case (no weights), IPW, and “Gwet” weighting via `calculate_inverse_weights()` (`code/misskappa/src/kappanp.cpp:22`).
  - Raw categorical uses unique-row compression to reduce kernel work (`code/misskappa/src/kappanp.cpp:71`).
  - Continuous uses explicit O(n²) loops (`code/misskappa/src/kappanp.cpp:375`).

- `emdiscrete` (EM engine): `code/misskappa/src/emdiscrete.h:1`, `code/misskappa/src/emdiscrete.cpp:1`
  - Provides `preprocess_raw`, `preprocess_counts`, `run_em` (`code/misskappa/src/emdiscrete.h:56`).
  - Computes an asymptotic covariance for `theta_hat` (Louis-style observed information) inside `run_em` (`code/misskappa/src/emdiscrete.h:45` and `code/misskappa/src/emdiscrete.cpp:312`).

- `kappaml` (kappa from EM output): `code/misskappa/src/kappaml.cpp:42`
  - Builds estimators as smooth functionals of `theta_hat` and uses `em_res.var` as input for delta-method variance (`code/misskappa/src/kappaml.cpp:169`).

- `kappaqp` (“quadratic” special path): `code/misskappa/src/kappaqp.cpp:111`
  - A separate codepath that computes moment-based estimates and a covariance matrix (ported legacy logic; see comments at `code/misskappa/src/kappaqp.cpp:24`).

- Shared weight/loss utilities:
  - Categorical weight matrices + continuous loss factories: `code/misskappa/src/common.cpp:1`

## Audit findings

### Public interface / consistency

- R docs don’t match implementation for method choices:
  - `kappa_raw()` and `kappa_continuous()` R wrappers include `"gwet"` (`code/misskappa/R/kappa.R:35`, `code/misskappa/R/kappa.R:101`), but `.Rd` usage does not list it (`code/misskappa/man/kappa_raw.Rd:6`, `code/misskappa/man/kappa_continuous.Rd:6`).
  - Continuous weights: C++ supports `"identity"` loss for continuous, but R’s `weight` choices exclude it (`code/misskappa/R/kappa.R:102` vs `code/misskappa/src/rcpp_interface.cpp:52`).
  - `kappa_counts()` has a dead code branch checking `if (method == "ipw")` even though `"ipw"` isn’t in the `method` choices (`code/misskappa/R/kappa.R:152` and `code/misskappa/R/kappa.R:181`).

- Naming in returned objects is fairly consistent:
  - Raw returns estimates named `Conger`, `Fleiss`, `Brennan-Prediger` (`code/misskappa/R/kappa.R:66`).
  - Continuous returns `Conger`, `Fleiss` (`code/misskappa/R/kappa.R:127`).
  - Counts returns `Fleiss`, `Brennan-Prediger` (`code/misskappa/R/kappa.R:193`).

### Code quality / maintainability

- Hard coupling: `code/misskappa/src/misskappa.h:4` includes `emdiscrete.h` and aliases `Result` to `emdiscrete::Result`. This makes EM a “core dependency” even for the non-EM paths; it increases the cost of cleanly separating estimator families later.

- Complexity hot-spots:
  - `kappanp::kappa_continuous` is O(n²·R²) style loops (`code/misskappa/src/kappanp.cpp:375`), unlike categorical which compresses unique patterns (`code/misskappa/src/kappanp.cpp:71`).
  - `kappaqp.cpp` is dense, low-level, and separate from the nonparametric logic (`code/misskappa/src/kappaqp.cpp:26`).

- Packaging hygiene:
  - The package directory contains interactive/session artifacts (`code/misskappa/.RData`, `code/misskappa/.Rhistory`).
  - `workspace/` is excluded from build (`code/misskappa/.Rbuildignore:10`) but contains compiled binaries and legacy experiments that can confuse maintenance.

### Tests (what exists vs what’s missing)

What exists now:

- Testthat tests that mostly check *cross-method equivalence* on complete data and *non-equality* on incomplete data:
  - Raw: `code/misskappa/tests/testthat/test-raw.R:23`
  - Continuous: `code/misskappa/tests/testthat/test-continuous.R:9`
  - Counts: `code/misskappa/tests/testthat/test-counts.R:16`

What’s missing (recommended “reasonable tests of the code itself”):

- Deterministic “known value” tests using tiny hand-constructed examples where you can compute kappas by hand (and also verify variances are finite/PSD).
- Stress tests for edge cases:
  - all missing; only one rater; one category; degenerate weights; raters with zero observations (IPW/Gwet should error cleanly) (`code/misskappa/src/kappanp.cpp:35`).
- Property-based checks (small random):
  - invariance to permuting raters for Fleiss/BP in counts/available-case settings.
  - `ipw` matches `available` when missingness is uniform MCAR (weights become constant).
- Performance/regression checks (lightweight):
  - a benchmark-style test (not for CRAN) ensuring `kappa_continuous` doesn’t blow up on moderate `n`.

### Documentation

- Package-level docs mention “bootstrap-based confidence intervals” (`code/misskappa/R/misskappa-package.R:5`), but there is no bootstrap implementation in the package code.
- Consider adding an `articles/` vignette or `README` inside the package explaining:
  - when to use `available` vs `ipw` vs `ml` vs `quadratic`,
  - the missingness assumptions (MCAR/PMCAR vs MAR) at a high level,
  - what the variance estimates represent.

## Suggested next steps (prioritized)

### 1) Public interface cleanup (low risk)

- Decide whether `"gwet"` is officially supported; if yes, expose it consistently and implement the missing continuous path (or remove it from choices).
- Align docs (`man/*.Rd`) with the actual `method` and `weight` options, regenerate with roxygen.
- Remove dead branches in R wrappers (e.g. `kappa_counts()` checking `"ipw"`).

### 2) Code-quality refactors (medium risk, high payoff)

- Decouple `misskappa::Result` and `kNaInteger` from `emdiscrete` by moving shared types to a small `result.h`/`types.h`, so estimator families can be trimmed independently later.
- Reduce duplication in `rcpp_interface.cpp` around:
  - category/value mapping,
  - loss/weight selection,
  - method dispatch.

### 3) “Real” tests / validation suite (medium risk)

- Add a new test file focused on small, exact examples (no randomness).
- Add targeted error/edge-case tests.
- Add at least one simulation-based consistency smoke test (small `n`, fixed seed) validating expected qualitative ordering (e.g., `ipw` vs `available`) under a controlled MCAR mechanism.

### 4) Performance and numerical stability (medium/high risk)

- Profile/optimize `kappa_continuous` (consider unique-row compression for repeated values or precomputations).
- Add guards and clearer error messages for degenerate denominators.
