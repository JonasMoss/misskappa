# Notes

Flat folder of legacy notes. Filenames carry a type prefix:

- `draft-*` — paper-style drafts.
- `note-*` — short focused notes.
- `scratch-*` — unfinished scratchpads.
- `ref-*` — reference / textbook-style theory notes.
- `impl-*` — implementation notes.
- `derive-*` — derivations.
- `verify-*` — verification / sanity checks.
- `example-*` — usage examples.
- `simulations-*` — auto-generated simulation tables.

## Drafts

- `draft-kappas-incomplete-data.lyx`: incomplete-data draft (EM machinery + Louis method).
- `draft-kappas-incomplete-data-v2.lyx`: older/alternate version.
- `draft-pairwise-complete-inference.lyx`: pairwise-complete / PMCAR draft.
- `draft-pairwise-inference-v2.lyx`: older/alternate draft, includes influence-function material.

## Notes

- `note-count-data-exchangeable-missingness.lyx`: count-data + exchangeable-missingness kernel used in EM-style derivations.
- `note-gwet-method-meaning.tex`: what `method = "gwet"` does in the current `misskappa` code.
- `note-horvitz-beta-ipw.lyx`: IPW / Horvitz-style scratch.
- `note-missing-compatible-distributions.lyx`: which joint distributions are compatible with a given missing-data mechanism.

## Scratch

- `scratch-em-ideas-and-derivations.lyx`: long EM/ML ideas-and-derivations scratchpad.
- `scratch-pairwise-estimator-definitions.lyx`: estimator-definition scratchpad for the pairwise drafts.

## References (semiparametric / theory)

- `ref-conditional-expectation-projections.lyx`
- `ref-semiparametric-efficiency.lyx`
- `ref-semiparametric-efficiency-v2.lyx`: older/alternate version.
- `ref-van-der-vaart-asymptotics.lyx`
- `ref-wasserman-all-of-statistics.lyx`

## Implementation

- `impl-misskappa-package-notes.lyx`: R-package implementation notes.
- `impl-rust-ffi-concept.lyx`: Rust/FFI concept notes.

## Derivations / verifications / examples

- `derive-eif-inference-kappa-raw.tex`: EIF-based inference formulas for the `kappa_raw()` backend.
- `verify-count-form-vs-irrcac-fleiss1971.tex`: count-form computation cross-checked against `irrCAC` on `dat.fleiss1971`.
- `example-misskappa-package-overview.tex`: short package-usage example write-up.

## Simulation tables

- `simulations-raw-three-estimators.tex`: combined three-estimators simulation write-up.
- `simulations-raw-three-estimators-main.tex`: auto-generated simulation table snippet (included in the legacy LyX manuscript via `\input{...}`).
- `simulations-raw-three-estimators-appendix.tex`: auto-generated extra simulation tables (appendix / online-supplement candidate).

## Misc

- `misskappa.bib`: BibTeX file used by the LyX drafts.
- `revision-plan.md`: internal working plan for revising the legacy `kappa-missing.lyx` manuscript.
