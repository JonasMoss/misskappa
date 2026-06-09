# misskappa repo TODO

Active repo-level backlog. Paper-specific TODOs live at
`papers/<slug>/dev/todo.md` (one per paper); the split decision and
section mapping live at `papers/split-plan.md`. The port
plan with the eight-step roadmap is at `dev/notes/port-plan.md`.

## Roadmap status

- [x] Step 1: skeleton + legacy move into `dev/legacy/`.
- [x] Step 2: port loss layer (`common.cpp` -> `src/loss.cpp` on Eigen + Result).
- [x] Step 3: available-case estimator + variance.
- [x] Step 4: IPW + Gwet estimators.
- [x] Step 5: FIML / EM.
- [x] Step 6: inference consolidation.
- [x] Step 7: R package wiring (Rcpp glue + S3 + sim list).
- [ ] Step 8: paper conversion (LyX -> .tex) + scripts wiring.

## Active backlog

- [ ] **R package release (GitHub-first): finish the unified-API refactor.**
      Decided 2026-06-06: ship a clean *estimators + CIs* surface, defer
      hypothesis tests, release on GitHub first (CRAN later). The public
      surface is narrowed to one entry point per input shape driven by a shared
      `estimator =` selector; the engine builds + runs (smoke-tested all six
      paths). State is mid-refactor — code works, tests/docs lag.

      *Done (mostly committed; `NAMESPACE` + `R/kappa.R` still dirty):*
      shared `estimator` vocabulary `pairwise | ipw | cat_fiml | nt_fiml` on
      `alpha()`, `kappa()`, `kappa_counts()`; `se_type` dropped (sandwich SE
      only); available-case + Gwet pulled from public `kappa()` (reachable for
      sims via internal `estimate_kappa_raw()`); default weight `identity` ->
      `nominal`; `kappa_continuous()` + `kappa_gwise()` demoted to internal
      (so public estimators are just `alpha`/`kappa`/`kappa_counts`);
      `kappa_quadratic_counts()` removed; `kappa_quadratic()` simplified to
      `(x, values)` empirical-IF only (normal/elliptical modes dropped, C++
      glue regenerated to 2-arg) and now returns IFs; new
      `kappa_quadratic_fiml()` (robust NT-FIML, quadratic kernel) backs
      `estimator = "nt_fiml"`; `confint(transform = "fisher")` delta-method CI
      + `print()` CI columns; `joint_vcov()`/`wald_test()` kept internal.

      *Open, next session, in order:*
      1. **Update the test suite to the new API** — ~31 old-API call sites in
         `test-kappa.R`, `test-influence.R`, `test-parity-irrcacsmoke.R`,
         `test-kappa-quadratic-fiml.R`: `method=` -> `estimator=`; drop
         `se_type` / `vcov=` / `relative_kurtosis`; route dropped public
         `available`/`gwet` through internal `estimate_kappa_raw()`; weight
         default is now `nominal`. R CMD check tests currently fail on these.
      2. **Fix stale roxygen cross-refs**, then re-roxygenise with the
         `pkgload` loader (bare `roxygenise()` chokes on the dataset string
         docs): `alpha_cat_fiml`/`alpha_continuous`/`kappa_quadratic`/
         `kappa_counts` descriptions still say `alpha(method=, type=)` /
         `kappa(method="available")`.
      3. **Sync `_pkgdown.yml`** — Estimators list still names now-internal
         `kappa_continuous` + `kappa_gwise`; reduce to `alpha`/`kappa`/
         `kappa_counts`.
      4. `just r-check` to confirm `R CMD check` is clean; commit
         `NAMESPACE` + `R/kappa.R` + doc fixes.
      5. Version decision: pre-1.0 reset (`0.5.0`) vs keep legacy `2.0.0`.

- [x] **GitHub-installability: vendor the C++ library into the R package.**
      `dev/vendor-cpp.sh` (= `just vendor`) copies the canonical `src/*.cpp`,
      `src/*.hpp`, and `include/misskappa/*.hpp` into `r-package/src/`, each
      stamped with a `@generated` banner (and cleaned up by that banner on the
      next run, so renames/deletes don't leave stale copies). `Makevars` is now
      self-contained — `PKG_CPPFLAGS = -I. -DEIGEN_NO_EXCEPTIONS`, C++23, no
      `PKG_LIBS`, no `../../include` — and R compiles every vendored `.cpp`.
      Eigen comes from `LinkingTo: RcppEigen`. Dropped `-fno-exceptions/-rtti`
      (Rcpp needs exceptions; the library never throws, so it's a no-op). The R
      justfile recipes depend on `vendor`, so editing only the canonical sources
      is enough; `just vendor-check` is the CI/pre-commit drift guard. Verified
      by `R CMD build` + `R CMD INSTALL` of the tarball into a throwaway lib with
      `MISSKAPPA_*` env scrubbed: `alpha()`/`kappa()` run, full testthat suite
      201/201. `opt`/`test-opt` still build the static lib for the C++ ctest.

- [ ] **Hypothesis-testing front door (deferred from the release).**
      `t.test`-style `kappa_test()` / `alpha_test()` returning `htest` over the
      already-built `joint_vcov()` / `wald_test()` engine. Use snake_case, not
      dotted names: `base::kappa` is a generic, so `kappa.test` would
      mis-dispatch as an S3 method. One-sample (`theta = theta0`) and
      independent two-sample (variances add) work for any estimator;
      paired/dependent needs stacked IFs via `joint_vcov()` and works for any
      fit whose `fit$psi` is non-NULL. Same-data method/weight/coefficient
      contrasts stay `wald_test()` territory. Exploratory runners:
      `experiments/studies/24-alpha-equal-cocron/`, `experiments/studies/25-kappa-equal-examples/`.
- [ ] **Measure alpha FIML feasibility before expanding scope.**
      Exact EM over the full fixed-category pattern space is plausible for
      small Likert batteries but scales as `C^R`. Add guardrails/documentation
      around feasible `C, R`, benchmark memory/time on representative `C x R`
      grids, and only then consider approximate or sparse EM variants
      (active-support EM, low-rank/composite-likelihood starts, or other
      approximations) as a separate follow-up.
- [ ] **Add an alpha experiment suite.**
      Simulate simple complete-data, MCAR, and MAR categorical item batteries;
      compare complete-data alpha, pairwise/available covariance alpha, and
      categorical FIML alpha for bias, SE calibration, and coverage. Add an
      applied example on `psych::bfi`, reporting one alpha per Big Five
      personality scale under the available-case and FIML paths, with clear
      notes about item scoring/reversal and the fixed-category support used.
      Starter categorical smoke/calibration runners now live in
      `experiments/probes/15-alpha-categorical-smoke/` and
      `experiments/probes/16-alpha-calibration-sweep/`. Experiment 16 has a
      capped `B = 200` calibration pass: `5^6, n = 4000` was too slow
      for an ordinary run, and `5^6` FIML retains downward bias through
      `n = 1000`. Independent saturated-EM checks agree with the package
      estimate, and the quadrature alpha truth agrees with large complete-data
      Monte Carlo, so treat the residual high-dimensional categorical-FIML
      bias as a finite-sample fact rather than an implementation blocker.
      Manuscript notes can state this matter-of-factly and leave debiasing,
      penalized/smoothed EM, bootstrap/jackknife corrections, or approximate
      sparse variants as optional follow-up. The normal-FIML comparison,
      applied `psych::bfi` example, alpha-specific Louis spectrum diagnostic,
      and explicit high-dimensional FIML guardrails remain open. The proposed
      paper-facing simulation grid is in
      `dev/notes/alpha-missing-simulation-study.md`: use an essential-tau
      setup derived from Zhang-Yuan's six-item calibration, include
      Zhang-Yuan's congeneric cell, add stronger congeneric and two-factor
      non-congeneric cells, borrow Enders/Savalei missingness rates, and use
      clean observed-anchor MAR rather than truth-dependent missingness.
- [ ] **Step 8: paper conversion (LyX -> .tex) + scripts wiring.**
      Continue the manuscript split/wiring work under the paper-local todo
      files once the current manuscript tree noise is settled.

## Deferred

- [x] **Add component-separable vector-valued pairwise agreement.**
      C++ now has internal pairwise/IPW estimators for rectangular
      `subjects x raters x features` data, passed as rater-major flattened
      matrices. The implemented loss class is component-separable with
      diagonal feature weights: Hamming, absolute, squared, and RMS. The R
      wrapper `kappa_vector()` is internal/unexported for experiments. The
      math and FIML boundary are in
      `dev/notes/component-separable-vector-kappa.md`; the CRACKLES pilot uses
      `experiments/workbench/26-crackles-vector-kappa/`.
- [x] **Add full-weight quadratic vector covariance route.**
      The internal R backend `kappa_vector_quadratic()` estimates Conger and
      Fleiss from `(mu, Sigma)` with a full symmetric PSD feature metric `W`,
      using either pairwise-available covariance moments or saturated
      normal-theory FIML. Tensor formulas live in
      `dev/notes/quadratic-vector-tensor-kappa.md`; the CRACKLES experiment
      compares pairwise covariance and NT-FIML under richer component
      missingness.
- [ ] **Categorical FIML for vector profiles.**
      Deferred. Component-wise missingness makes full-profile categorical
      equality and general finite-vector losses latent. A saturated
      multinomial FIML over rater-feature profiles is the natural extension,
      but the support is combinatorial and needs a separate design pass.
- [ ] **POD-pointer overloads — only when a non-Eigen consumer shows up.**
      The Eigen API is already thin for the two consumers we have (in-tree
      C++ tests and R via Rcpp+RcppEigen). If we ever add a Python / Julia /
      standalone-CLI binding, add overloads in
      `include/misskappa/estimate.hpp` that take
      `(const int* ratings, int n, int R, const double* W, int C)` and write
      results into caller buffers, with the Eigen versions becoming
      one-liners over `Eigen::Map`. Until then, not worth the surface
      doubling.
- [ ] **Re-evaluate counts-format and continuous-rating scope after submission.**
      Counts-format and continuous-rating estimators are currently present in
      the library and R package. Revisit their public-surface priority after
      the paper submission path is settled.

## Done

- [x] **Prototype coefficient alpha estimators.**
      `estimate_alpha_available()` computes the scored-categorical
      pairwise-covariance alpha baseline with per-subject influence rows, and
      `estimate_alpha_available_continuous()` exposes the same MCAR sandwich
      path for numeric item scores. `estimate_alpha_fiml()` reuses the raw
      categorical EM / Louis machinery to fit the saturated full-response
      distribution and map theta to alpha. The R package exposes numeric alpha
      through `alpha(method = c("available", "fiml"))` and the finite-category
      multinomial EM path through `alpha_cat_fiml()`. Counts input,
      feasibility benchmarking, approximate EM, and experiments remain separate
      follow-ups.
- [x] **Settle the inference-surface policy.**
      `misskappa_estimate` now has a consistent R-facing inference contract:
      all estimators with subject-level influence rows carry them in `fit$psi`
      for `joint_vcov()` and multi-fit `wald_test()`, and
      `stats::influence(fit)` is registered as the public S3 accessor for the
      same matrix.
- [x] **Expose influence functions for FIML estimators.**
      Raw FIML and counts-format FIML now return per-subject influence
      functions from the same reduced observed-score / Louis-information
      path used for their vcov. C++ and R tests assert that
      `crossprod(psi) / n^2` reconstructs the reported covariance, and
      `joint_vcov()` / multi-fit `wald_test()` can include FIML fits.
- [x] **Expose influence functions for counts and continuous estimators.**
      Counts-format available-case and continuous MCAR estimators now return
      per-subject IF matrices through the shared `misskappa_estimate` surface,
      so `joint_vcov()` and `wald_test()` can combine them with other
      same-subject fits.
- [x] **Add Wald tests for misskappa estimates.**
      `wald_test()` now tests single-fit and cross-fit linear hypotheses
      against `misskappa_estimate` objects. Single-fit tests use `vcov()`;
      multi-fit tests use `joint_vcov()` and therefore require aligned
      per-subject influence functions.
- [x] **Keep ASan opt-in.**
      The default `dev` preset stays a portable debug build without
      sanitizers. The `asan` preset remains available for targeted local runs
      once the required libasan / libubsan runtimes are installed; making it
      the default would make ordinary build/test loops too host-dependent.
- [x] **Add `.clang-format` and a `just fmt-cpp` recipe.**
      The formatter config follows the active C++ style with two-space
      indentation and a 100-column limit. The recipe formats active C++
      surfaces only (`include/`, `src/`, `tests/unit/`, `r-package/src/`),
      leaving legacy and vendored third-party code untouched.
- [x] **Wire CI.**
      `.github/workflows/ci.yml` now runs a C++ dev build/test job and an
      R package check job on Ubuntu 24.04. The R job builds the opt static
      library, installs the local `dev/irrcacsmoke` oracle, and runs
      `R CMD check --no-manual` against `r-package/`.
- [x] **Add closed rectangular g-wise / Frechet agreement.**
      `estimate_gwise` and `kappa_gwise()` now cover complete
      subjects-by-raters designs with Frechet nominal, Frechet absolute,
      Frechet quadratic, and Hubert disagreement kernels. They return
      Cohen-type and Fleiss-type chance-corrected coefficients with
      influence-function vcov. Categorical chance terms use exact finite
      support enumeration; continuous distances remain direct over item
      tuples. Missing-rater support remains intentionally out of scope.
- [x] **Remove the O(n^2) raw-estimator chance kernels.**
      `estimate_available`, `estimate_ipw`, and `estimate_gwet` now compute
      chance-disagreement V-statistic row sums, column sums, and totals from
      weighted category/rater totals in `src/estimate_raw.cpp`, preserving the
      existing influence-function covariance without looping over all subject
      pairs.
