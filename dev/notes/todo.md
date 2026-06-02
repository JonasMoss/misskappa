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

- [ ] **Step 8: paper conversion (LyX -> .tex) + scripts wiring.**
      Continue the manuscript split/wiring work under the paper-local todo
      files once the current manuscript tree noise is settled.

## Deferred

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

- [x] **Settle the inference-surface policy.**
      `misskappa_estimate` now has a consistent R-facing inference contract:
      all estimators with subject-level influence rows expose them through
      `influence()` for `joint_vcov()` and multi-fit `wald_test()`, while
      closed-form quadratic estimators intentionally return `NULL` because
      their covariance is built at the reduced moment-summary level.
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
