# paper/AGENTS.md

Profile: Math

Target journal: Psychometrika.

## Scope

This manuscript presents two estimators for weighted kappa coefficients
(Cohen / Fleiss / Brennan-Prediger) with incomplete categorical ratings:

- Available-case (pairwise complete observations).
- IPW / Hajek (inverse-probability-weighted).

Both are MCAR estimators. The paper proves consistency under MCAR, gives
joint inferential theory via influence functions, and shows that Gwet's
earlier inferential method (the comparison baseline) is inconsistent. A
small simulation study under MCAR / MAR illustrates the order-of-magnitude
differences in finite samples; a worked empirical example demonstrates the
companion R package.

FIML / EM is implemented in the `misskappa` library and used in the
simulation, but the efficient-missing-data story (EM efficiency, sandwich
robustness) is deliberately deferred to a follow-up paper.

## Conventions

- Density: dense, methods-developer voice. No tutorials, no padding.
- Citations: Tsiatis (2006), van der Vaart (1998), and Robins / Rotnitzky
  / Zhao for missing-data foundations — three sentences in the intro, no
  more. Plus core kappa definitions (Cohen, Fleiss, Brennan-Prediger) and
  the Gwet / irrCAC stack.
- Length: target ~12 published pages (Math profile lower quartile for
  Psychometrika). Online supplement carries heavier derivations.
- Voice: impersonal ("The estimator solves..."), past tense for results.
- Style: see `paper/STYLE.md`. Booktabs always; Okabe-Ito palette;
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
