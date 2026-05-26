# paper TODO

Paper-specific TODO. Repo-level TODO is in `../../dev/notes/todo.md`.

## Conversion

- [ ] One-shot convert `dev/legacy/kappa-missing.lyx` -> `paper/kappa-missing.tex`
      via `lyx -e latex` plus cleanup pass. Document the conversion in
      `paper/dev/notes/`.
- [ ] Move relevant content from `dev/legacy/notes/{em,mcar-ipw,semiparametrics,gwet-case.tex,inference-eif-raw.tex}`
      into `paper/dev/notes/`.
- [ ] Port `dev/legacy/code-analysis/simulations_raw_three_estimators.R` to
      `paper/scripts/`, adjusted to write raw output to `results/raw/` and
      curated CSVs to `results/`.

## Manuscript

- [ ] Fix the truncated "We w" sentence in the Introduction.
- [ ] Fill in Section 2.3 "Computation".
- [ ] Fill in Section 2.4 "Illustrations".
- [ ] Fill in "Concluding remarks".
- [ ] Tighten the abstract phrasing once the rest of the paper settles.
- [ ] Add a Symbol Glossary section near the front.
- [ ] Expand bibliography: missing-data foundations (MCAR/MAR, IPW/Hajek),
      Tsiatis / van der Vaart / Robins-Rotnitzky-Zhao, U-statistics for
      asymptotic variance.

## Build infrastructure

- [ ] Wire `paper/justfile` recipes: `pdf`, `sims`, `tables`, `figures`,
      `paper`, `archive`, `clean`.
- [ ] Add `paper/scripts/build_tables.R` and `paper/scripts/build_figures.R`
      stubs.
- [ ] Add `paper/tables/<slug>_stats.tex` macro file with `??` fallbacks.
