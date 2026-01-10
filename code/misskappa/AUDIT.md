# misskappa package audit

This is a living engineering note for the `misskappa` R package (`code/misskappa/` in this repo). It documents what exists today and what to improve next, without committing to removing major subsystems (EM / “quadratic”) yet.

## Roadmap (cleanup + portability)

Primary goals:

1. **Simple R package**: easy defaults, clear method choices, stable return shapes, good error messages.
2. **Portable C++ core**: isolate algorithmic code from R/Rcpp so it can be reused for a future Python package (e.g. via pybind11) with minimal changes.

### Milestone 0: hygiene and baseline checks (done)

- [x] Add repo workflows (`Justfile`) and make `R CMD check` / tests runnable in one command.
- [x] Remove tracked session artifacts (`.RData`, `.Rhistory`, `.Rproj.user`) and add ignore rules.
- [x] Align docs with R wrappers and fix continuous `"gwet"` dispatch in C++.

### Milestone 1: R package UX (low risk, high impact)

Deliverables / tasks:

- [x] Decide and document which `method` values are *official* for each function (especially `"gwet"`).
  - Decision: `"gwet"` is official but positioned as compatibility/comparison rather than the primary recommended interface.
- [x] Add lightweight user-facing guidance: “which method should I use?” and “what assumptions does it correspond to?” (README or vignette).
- [x] Add a small S3 class (e.g. `"misskappa_estimate"`) with `print()` and `as.data.frame()` so results are easy to use downstream.
- [x] Correct package-level documentation to match what’s implemented (no phantom “bootstrap CI” claim).
- [x] Consider a single top-level entrypoint `kappa()` that dispatches on input type (raw vs counts vs continuous), while keeping `kappa_*()` functions.
  - Implemented `kappa(x, type = c("auto","raw","continuous","counts"), ...)` with conservative heuristics for `"auto"`.
- [ ] Standardize argument naming/behavior across functions (especially how `gwet` is exposed: method vs option).
  - Current decision: keep `"gwet"` as a `method` value (official, compatibility/comparison), rather than introducing a separate option at this stage.
  - Counts: the counts format aggregates away per-rater missingness patterns, so `"ipw"`/`"gwet"` are not currently supported for `kappa_counts()`.
- [x] Add minimal examples for each `kappa_*()` function.

### Milestone 2: isolate a portable C++ “core” (medium risk, foundational)

Target architecture:

- **Core library** (pure C++): no `Rcpp`, no `R` headers, no `NA_INTEGER`, no `Rcpp::stop`; expose a small API and return errors via a portable `Result<T>` / `StatusOr<T>` type.
- **Bindings layer** (Rcpp): translate R matrices/options ↔ core types, translate errors ↔ `Rcpp::stop`.

Concrete steps:

- [ ] Create `result.h` / `types.h` defining a portable `Result<T>` and common types (so `misskappa.h` no longer aliases `Result` from `emdiscrete`).
- [ ] Make `misskappa.h` self-contained (no EM includes); treat EM as one backend behind a clean API.
- [ ] Define core entrypoints (names illustrative): `misskappa::{raw,counts,continuous}::estimate(...)` returning a single `Estimation` struct (estimates + vcov + metadata).
- [ ] Decide on a portable public data model (avoid exposing Armadillo types if Python reuse is a priority).

### Milestone 3: style, tooling, and CI (low/medium risk)

Goals:

- Consistent formatting and a “boring” codebase: fewer bespoke conventions; easier diffs; easier review.

Suggested tooling:

- [ ] C++: add `.clang-format` and a `just fmt-cpp` recipe; run formatting in CI.
- [ ] R: add `styler` / `lintr` config and `just fmt-r` / `just lint-r` recipes (optional).
- [ ] CI: add a minimal GitHub Actions workflow running `R CMD check`, plus a C++ compile-only job with different compilers/standards.

### Milestone 4: simplify/trim estimator families (higher risk, optional later)

Once Milestone 2 exists (portable core + clean dispatch), it becomes much cheaper to:

- Remove or quarantine the bespoke `"quadratic"` path if it’s not needed.
- Split EM (`ml`) into an optional component (or separate package) if you decide it’s too heavyweight.
- Replace dense legacy code with clearer implementations backed by the test suite.

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

- Recent fixes made to align docs and implementation:
  - `"gwet"` is now documented in `.Rd` for `kappa_raw()` and `kappa_continuous()` and is reachable in the C++ continuous backend.
  - Continuous `"identity"` weight is now exposed in `kappa_continuous()`.
  - Removed dead code in `kappa_counts()` that checked for `"ipw"`.

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
  - The repo previously tracked interactive/session artifacts (`code/misskappa/.RData`, `code/misskappa/.Rhistory`, `code/misskappa/.Rproj.user`); these should stay untracked and ignored.
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

- Decide whether `"gwet"` is officially supported and, if so, whether it should remain a `method` or become an option (e.g. `reweight = "gwet"`).
- Add minimal usage examples and a short “method selection” guide (README/vignette).
- Add a small S3 class + printer for stable downstream consumption (optional, but helps “easy-to-use”).

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
