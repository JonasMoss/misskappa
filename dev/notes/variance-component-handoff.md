# Handoff — finish the compound-symmetry ("variance-component") FIML covariance

**For:** Codex. **Branch:** `ten-hove`. **Date:** 2026-06-14.

## What this is

Commit `ff96f5c` added a structured within-subject covariance option to the
normal-FIML path: `em_options = list(covariance = "variance_components")` fits
`Sigma = psi*J + e*I` (compound symmetry) instead of the saturated covariance,
for quadratic `kappa()` (Conger/Fleiss) and continuous `alpha()`. It is a
saturated EM with a compound-symmetry projection in the M-step
(`project_variance_components` in `src/estimate_normal_fiml.cpp`). This note is
the finishing work that investigation (probes 39, 40) showed is needed.

## Already verified — do NOT redo

- Point estimate == the fixed-rater MLE: matches `lmer(value ~ 0 + item +
  (1|subject), REML=FALSE)` to **5e-8** on MISSING data (MCAR + MAR, well- and
  mis-specified). See `experiments/probes/40-alpha-variance-component-ci/`.
- Probe 39 (`experiments/probes/39-ten-hove-icc/`): lme4 / magmaan SEM parity,
  delta-SE vs parametric bootstrap.
- C++ unit tests green: `tests/unit/normal_fiml_test.cpp` (12 cases). The IF
  reconstructs the vcov; the identifiability guard relaxes from complete-pair
  (saturated) to *connected* co-observation (compound symmetry).
- Alpha CI behaviour (probe 40): MCAR -> consistent + best-calibrated even under
  misspecification (alpha = f(tr,sum) and the CS projection preserves both); MAR
  + parallel -> consistent; MAR + congeneric -> inconsistent (saturated is the
  fix); pairwise is MAR-invalid. This is the documented behaviour, not a bug.

## SCOPE — read before touching anything

This pack is **Conger / Fleiss / Brennan-Prediger** (and coefficient alpha) for
agreement among the **fixed set of raters in hand**, under missing data. The two
within-subject covariance models are *saturated* (`nt_fiml`) and *compound
symmetry* (this option). Inference is the i.i.d.-subjects influence-function
sandwich, which is the **correct, honest SE** for these fixed-rater estimands.

**Explicitly OUT of scope — do not implement:**
- Random-rater ICC(A,1)/ICC(C,1), a rater variance component `sigma^2_r`,
  G-theory generalizability, crossed / cross-classified random effects, or any
  "generalize to a rater universe" SE. (The crossed model breaks subject
  independence and would need a multiway subject x rater cluster-robust SE; we
  are not going there.) The compound-symmetry option is the *within-subject*
  covariance only; the rater spread is carried by the free mean vector, which is
  exactly right for fixed-rater Conger/Fleiss/BP.

## Tasks

### 1. Fix the naming / docs (the "ten Hove" conflation) — DO THIS

The option is currently aliased `"ten_hove"` and docstrings call it "the scalar
ten Hove covariance model". That misleads: ten Hove's recommended estimator
(MLE-RE) is the *random-rater* crossed model with three variance components,
which this is NOT. This is the fixed-rater compound-symmetry within-subject
covariance.

- In `r-package/R/kappa_quadratic_fiml.R` `.normal_covariance_model()`: drop the
  `"ten_hove"` alias (branch is unreleased, no deprecation cycle needed). Keep
  `"variance_components"` as the canonical string (it matches the C++ enum
  `NormalCovarianceModel::variance_components` and is accurate for the
  within-subject subject-variance + residual decomposition).
- Replace every "ten Hove covariance model" docstring phrase with
  "compound-symmetry (fixed-rater) covariance model `Sigma = psi J + e I`" in:
  `r-package/R/kappa_quadratic_fiml.R`, `r-package/R/alpha_continuous.R`,
  `r-package/R/kappa.R` (search: `ten Hove`, `ten_hove`). Optionally add one
  sentence: "treats the raters as fixed; it is not ten Hove's random-rater
  MLE-RE / ICC." Re-roxygenise with `just r-docs` so `man/*.Rd` match.

### 2. Brennan-Prediger parity in the FIML path — DO (in scope, confirmed)

BP is included in every estimator in this pack; the FIML path is the gap. The
empirical quadratic path returns **Conger / Fleiss / Brennan-Prediger**
(`r-package/R/kappa.R:274`), but the normal-FIML path
(`scalar_quadratic_gradient` in `src/estimate_normal_fiml.cpp`) returns only
**Conger / Fleiss**. Add BP to the FIML path for both the saturated and the
compound-symmetry covariance models.

Design note (this is why it is not a one-liner): BP's chance term is the
**uniform** expected disagreement, not the observed marginals. For quadratic
weights it depends only on the rating *support*, not on `(mu, Sigma)`:
- continuous scale `[lo, hi]`: `D_e^BP = (hi - lo)^2 / 6`;
- category values `v_1..v_q`: `D_e^BP = (1/q^2) * sum_{a,b} (v_a - v_b)^2`.

The observed disagreement `D_o` IS a smooth function of `(mu, Sigma)` (same
moments as Conger/Fleiss), so `BP = 1 - D_o / D_e^BP` slots straight into the
existing gradient/delta-method machinery — the denominator is a constant, so
`grad(BP) = grad(D_o) / D_e^BP`. The only new requirement is a `values` / scale
argument on the FIML entry point, mirroring the empirical path's `values`. Plumb
it through `kappa_quadratic_fiml()` -> `rcpp_*` -> the gradient functor, return a
3-vector named `c("Conger","Fleiss","Brennan-Prediger")`, extend the IF/vcov to
3 rows, and add a unit test pinning BP against the empirical quadratic on
complete data.

### 3. Folding / public-surface decision — leave as-is unless told otherwise

Per probe 39's decision the compound-symmetry option is positioned as the
quadratic/Conger paper's online-appendix comparator, not a headline recommended
estimator. Keep the current `em_options` switch; do not promote it to a
top-level `estimator=` value or change defaults.

## Validation & tests (rough guide)

Build/run loop (edit canonical `src/`, never the vendored `r-package/src/`):
- `just opt` — C++ opt build + ctest (runs `tests/unit/normal_fiml_test.cpp`).
- `just r-dev` — fast R reinstall (links the prebuilt `.a`).
- `just vendor-check` — drift guard; must pass (vendored mirror == `src/`).
- `just r-docs` — regenerate `man/*.Rd` after docstring edits.
- `just r-check` — `R CMD check`; must stay clean.

**Oracles (all in-repo / installed — no new deps):**
- `lme4` fixed-rater MLE `lmer(value ~ 0 + item + (1|subject), REML=FALSE)` —
  oracle for the compound-symmetry point estimate (already 5e-8, probe 40).
- The **empirical quadratic path** (`kappa(x, weight="quadratic")`, default
  estimator) — the trusted in-repo BP/Conger/Fleiss implementation; it is the
  oracle for task 2 on complete data.
- `magmaan` one-factor SEM — secondary parity check (probe 39).

**Task 1 (naming) — checks:**
- `grep -rn "ten_hove\|ten Hove" r-package/R man` returns only the deliberate
  "not ten Hove's MLE-RE" disclaimer (no alias, no "ten Hove covariance model").
- `covariance="variance_components"` still works; `covariance="ten_hove"` now
  errors. Update/remove any test asserting the old alias
  (`tests/testthat/test-kappa-quadratic-fiml.R`, `test-alpha-continuous.R`).
- `just r-docs` + `just r-check` clean; `just vendor-check` passes.

**Task 2 (BP) — checks:**
- New unit test (C++ in `normal_fiml_test.cpp`, mirrored in testthat): on
  COMPLETE data, FIML BP == empirical-quadratic BP to ~1e-8, for a fixed
  `values`/range, under both saturated and compound-symmetry covariance.
- Extend `check_influence_contract` to **k=3**: `crossprod(psi)/n^2 == vcov`
  for the 3-coefficient (Conger/Fleiss/BP) IF block.
- Missing-data smoke: BP finite + `vcov` finite under missingness (both
  covariance models), connected-pattern case for compound symmetry.
- Regression: re-run `experiments/probes/40-alpha-variance-component-ci/
  parity_and_coverage.R` (alpha must be unchanged); optionally add a BP row to
  the kappa coverage to confirm it tracks Conger/Fleiss.
- Scope guard: confirm no `sigma^2_r` / ICC(A,1) / ICC(C,1) code crept in.

## Pointers

- Estimator core + gradients: `src/estimate_normal_fiml.cpp`
  (`project_variance_components`, `scalar_quadratic_gradient`, `alpha_gradient`,
  `variance_component_jacobian`). Vendored mirror under `r-package/src/` is
  regenerated by `just vendor`; edit the canonical `src/` copy.
- R wrappers: `kappa_quadratic_fiml.R`, `alpha_continuous.R`, dispatch in
  `kappa.R`.
- Identifiability guard: `src/detail_pattern_checks.hpp`
  (`require_connected_observation`).
- Tests: `tests/unit/normal_fiml_test.cpp`; R: `tests/testthat/
  test-kappa-quadratic-fiml.R`, `test-alpha-continuous.R`.
- Evidence: `experiments/probes/39-ten-hove-icc/NOTES.md`,
  `experiments/probes/40-alpha-variance-component-ci/NOTES.md`.
