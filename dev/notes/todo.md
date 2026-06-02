# misskappa repo TODO

Active repo-level backlog. Paper-specific TODOs live at
`papers/<slug>/dev/todo.md` (one per paper); the split decision and
section mapping live at `papers/split-plan.md`. The port
plan with the eight-step roadmap is at `dev/notes/port-plan.md`.

## Currently in flight

- [x] Step 1: skeleton + legacy move into `dev/legacy/`.
- [x] Step 2: port loss layer (`common.cpp` -> `src/loss.cpp` on Eigen + Result).
- [x] Step 3: available-case estimator + variance.
- [x] Step 4: IPW + Gwet estimators.
- [x] Step 5: FIML / EM.
- [ ] Step 6: inference consolidation.
- [x] Step 7: R package wiring (Rcpp glue + S3 + sim list).
- [ ] Step 8: paper conversion (LyX -> .tex) + scripts wiring.

## Backlog

- [ ] **POD-pointer overloads — only when a non-Eigen consumer shows up.**
      The Eigen API is already thin for the two consumers we have (in-tree
      C++ tests and R via Rcpp+RcppEigen). If we ever add a Python / Julia /
      standalone-CLI binding, add overloads in
      `include/misskappa/estimate.hpp` that take
      `(const int* ratings, int n, int R, const double* W, int C)` and write
      results into caller buffers, with the Eigen versions becoming
      one-liners over `Eigen::Map`. Until then, not worth the surface
      doubling.
- [ ] Decide whether to enable the `asan` preset by default once Fedora's
      libasan / libubsan runtime packages are installed locally.
- [ ] Add `.clang-format` and a `just fmt-cpp` recipe.
- [ ] Re-evaluate whether counts-format input and continuous ratings should
      come back into the library after the paper is submitted.

## Done

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
