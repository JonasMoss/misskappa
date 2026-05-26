# TODO for manuscript

For a structured revision plan (targets + hard/soft TODOs), see `notes/revision-plan.md`.

## Content (LyX reality check)

- [ ] Fix the truncated sentence “We w” in the Introduction (`kappa-missing.lyx` around line ~354).
- [ ] Fill in Section 2.3 “Computation” (currently a stub heading).
- [ ] Fill in Section 2.4 “Illustrations” (currently a stub heading).
- [x] Integrate simulations into the manuscript (LyX now inputs `notes/simulations-raw-three-estimators-main.tex` in the Simulations section; regenerate via `Rscript code/analysis/simulations_raw_three_estimators.R`).
- [ ] Fill in “Concluding remarks” (currently a stub section heading; bibliography starts immediately).
- [x] Count-form discussion exists (currently a subsection under “Other estimators”) with a Fleiss-1971 example; decide if it should be promoted to a main section and add any missing estimation/inference/computation detail.
- [x] `misskappa` package section exists with raw + count examples; consider adding `coef()`/`vcov()` output and/or mentioning the `kappa()` dispatcher for the “one interface” story.
- [ ] Gwet/irrCAC section exists; add citations for Gwet/irrCAC and ensure the definition matches what the package actually computes (and what we claim about (in)consistency).

## Cleanup

- [ ] Proofread for obvious rough edges/typos (e.g. abstract grammar: “and hold …”, “well-known form” → “well-known from”, etc.).
- [ ] Expand bibliography (currently very small): add missing-data basics (MCAR/MAR), Hájek/IPW references, and anything needed to justify inference machinery without over-citing.
- [ ] Decide bibliography/citation style for Psychometrika (LyX currently uses `natbib` + `plainnat`; switch if you truly need APA/apacite).
- [ ] Housekeeping: treat `kappa-missing.lyx` as the source of truth; `kappa-missing.pdf` is a compiled artifact and may lag behind.
