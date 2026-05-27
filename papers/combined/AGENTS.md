# papers/combined/AGENTS.md

Profile: Math

Target journal: Psychometrika.

## Status (2026-05-27): combined draft being split into three papers

This directory holds the **combined source-of-truth draft** of the
original kappa-missing manuscript. As of 2026-05-27 the work is being
broken into three independent papers, each with its own sibling
directory at the repo root:

- `papers/ipw/` — practical IPW paper, target Psychometrika
- `papers/fiml/` — semiparametric-efficient estimation under MAR
- `papers/quadratic/` — closed-form quadratic / Lin's CCC

These spinoffs are intended for simultaneous submission as a coordinated
paper family. Cross-references among Papers A/B/C are allowed and expected;
use provisional labels or "companion paper" language until final titles and
bibliography keys are settled.

The decision and section mapping live in `papers/combined/dev/split-plan.md`.
New work should go in the appropriate `papers/<slug>/` directory; this
directory remains canonical only until each section has been migrated
and marked `% MIGRATED to papers/<slug>/`. Do not delete content from
`papers/combined/kappa-missing.tex` until the corresponding spinoff compiles
standalone.

If you are arriving cold and would otherwise edit `papers/combined/kappa-missing.tex`,
check `papers/combined/dev/split-plan.md` first to find the spinoff directory
that owns the section you care about.

## Scope

This manuscript presents four estimators for weighted kappa coefficients
(Cohen / Fleiss / Brennan-Prediger) with incomplete categorical ratings,
organised by data form and missingness assumption:

- **Raw (rater-identified) data.** Hájek-type IPW (MCAR) and FIML / EM
  (MAR). IPW is consistent under MCAR with no parametric assumption on the
  rater joint; FIML is consistent under MAR and attains the semiparametric
  bound under the discrete-data full model.
- **Counts data with exchangeable raters.** Available-on-counts (the
  pairwise plug-in, which under exchangeability coincides with IPW) and
  FIML / EM (hypergeometric completion of partial counts).

The available-case (raw), Gwet, and closed-form quadratic estimators
appear as a short comparison section. Available-case is consistent under
narrower conditions and inefficient; Gwet is generally inconsistent; the
quadratic estimator is closed-form efficient in the two-rater / pairwise-
complete case (`E[l(X^*, X^*) | X]` constant in `X`).

The paper develops consistency and inference (influence functions for
IPW, Louis observed information for FIML), shows that Gwet's reweighting
is inconsistent under MAR, and demonstrates the bias / efficiency
ordering in a small simulation study. A worked empirical example
exercises the companion R package.

## Conventions

- Density: dense, methods-developer voice. No tutorials, no padding.
- Citations: Tsiatis (2006), van der Vaart (1998), and Robins / Rotnitzky
  / Zhao for missing-data foundations — three sentences in the intro, no
  more. Plus core kappa definitions (Cohen, Fleiss, Brennan-Prediger) and
  the Gwet / irrCAC stack.
- Length: target ~12 published pages (Math profile lower quartile for
  Psychometrika). Online supplement carries heavier derivations.
- Voice: impersonal ("The estimator solves..."), past tense for results.
- Style: see `papers/combined/STYLE.md`. Booktabs always; Okabe-Ito palette;
  every cited number flows from `tables/<slug>_stats.tex`.

## Build

```sh
just pdf            # build the manuscript PDF
just sim <name>     # run one simulation
just sims           # run all simulations
just tables         # rebuild tables/ from results/
just figures        # rebuild figures/ from results/
just paper          # tables + figures + pdf (no sims)
just archive        # zip the archive-bound subset for OSF
just clean          # remove LaTeX build products
```

Heavy work lives under `scripts/`. Cheap recipes reuse existing
`results/` and complete in seconds.

## Reproducibility

Each simulation runner writes raw per-replicate output to
`results/raw/<run-id>/`; a summarising step writes curated CSVs into
`results/`. Tables and figures read CSVs, never call simulation code.
Curated CSVs are tracked; raw is ignored.
