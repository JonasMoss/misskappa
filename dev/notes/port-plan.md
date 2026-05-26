# Port plan

Roadmap for restructuring `misskappa` from an R-package-in-paper-repo into
a standalone C++23 library with R bindings and a paper as subprojects.
Approved plan stub at `~/.claude/plans/ok-so-its-been-glowing-duckling.md`.

## Eight steps

1. **Skeleton.** New tree (`include/`, `src/`, `tests/`, `r-package/`,
   `paper/`, `dev/`); legacy moved under `dev/legacy/`; empty
   `libmisskappa.a` builds; sanity doctest passes. **Done.**
2. **Loss layer.** Port `common.cpp` math (eight categorical weight matrices,
   five continuous loss kernels) to `src/loss.cpp` on Eigen + `Result<T>`.
   `tests/unit/loss_test.cpp` exercises symmetry, zero diagonal, and
   closed-form known values per weighting.
3. **Available-case + variance.** Port the categorical path from
   `dev/legacy/misskappa/src/kappanp.cpp`. Add irrCAC oracle fixture in
   `tests/fixtures/` via `tests/tools/regen_oracle.R`.
4. **IPW + Gwet.** Port the IPW + Gwet branches; property test that IPW
   under uniform MCAR collapses to available-case.
5. **FIML / EM.** Port `emdiscrete.{h,cpp}` + `kappaml.cpp` to
   `src/estimate_fiml.cpp` and `include/misskappa/em.hpp`. Bounded
   iteration; `Error::not_converged` on non-convergence. Louis
   observed-information + delta-method variance.
6. **Inference consolidation.** Move shared influence-function / delta-method
   plumbing into `src/inference.cpp` behind `include/misskappa/inference.hpp`.
7. **R package wiring.** `r-package/src/Makevars` -> prebuilt
   `../../build-opt/libmisskappa.a`. Public R: `kappa(x, method, weight, ...)`,
   S3 print/coef/vcov/as.data.frame/confint, `sim` list, datasets.
8. **Paper conversion.** LyX -> .tex; relevant `dev/legacy/notes/*` -> `paper/dev/notes/`;
   `simulations_raw_three_estimators.R` -> `paper/scripts/`; results -> `paper/results/`.

## Math kept from legacy

The new tree is rewritten on Eigen + `Result<T>`; legacy code provides the
math. Source-of-truth files in `dev/legacy/misskappa/src/`:

- `common.cpp` — eight weight matrices + five continuous losses.
- `kappanp.cpp` — `calculate_inverse_weights` (Gwet/IPW), unique-row
  compression for raw, variance formulas.
- `emdiscrete.{h,cpp}` — `preprocess_raw`, `run_em`, Louis observed
  information.
- `kappaml.cpp` — delta-method variance over EM output.
- `R/` — S3 layer, `kappa()` dispatcher, dataset wrappers (mostly portable
  with path rewrites and a Makevars change).

## Dropped from scope (Phase 1)

- `kappaqp.cpp` — legacy quadratic special path. Not on the critical path
  for the paper. Re-evaluate after submission.
- `kappa_counts()` and counts-format input.
- `kappa_continuous()` for continuous ratings.

All three can come back later if a downstream need surfaces.
