# paper TODO

Paper-specific TODO. Repo-level TODO is in `../../dev/notes/todo.md`.

The 2025 rewrite refocused the manuscript on four estimators across two data
forms: (raw, counts) × (FIML, IPW). The available-case, Gwet, and quadratic
estimators are now a short comparison section. See `paper/AGENTS.md` for the
declared scope.

## Open

### Simulations

- [ ] Run the full-grade simulation (`SIM_FULL=1 Rscript
      scripts/simulations_raw_three_estimators.R`). Smoke-run output is fine
      for layout, but the Psychometrika-bound numbers need `B_mod=500`,
      `B_big=200`. Note: includes FIML; expect tens of minutes to hours.
- [ ] Add a counts simulation script (`scripts/simulations_counts.R`)
      mirroring the raw script for the available-on-counts vs counts-FIML
      comparison. Single MAR-counts cell suffices.
- [ ] Decide on a final main-text table shape (4 methods × 3 DGPs gets wide;
      consider splitting into one MCAR cell and one MAR cell).

### Pre-paper exploration (in `experiments/`, not bound to manuscript yet)

Open questions about the current sim design before committing the paper to
specific DGPs / claims. Each is one focused experiment folder per
`experiments/AGENTS.md`. Findings may motivate manuscript changes, may
soften / sharpen specific claims, or may end in the trash — that's the
point of running them outside `paper/`.

- [ ] `01-coverage-iif-louis` — Do IPW influence-function and FIML Louis
      Wald CIs hit nominal coverage at moderate n? Single column added to
      the bias/SD table answers the obvious referee question. Cheap, both
      variance estimators are already implemented.
- [x] `02-rater-model-sensitivity` — Do A/B/C bias and efficiency orderings
      survive a Dawid-Skene rater model (per-rater confusion matrix)? The
      current latent-truth-plus-guess model is easy to describe but very
      constrained. Runner/report landed; the `n = 1000`, `B = 50` run
      showed unstable orderings, so choose the manuscript DGP deliberately.
- [ ] `03-ac-vs-ipw-efficiency` — Where exactly does AC inefficiency
      become visible? Sweep π-variability against rater-exchangeability.
      Smoke run currently shows IPW SD > AC SD under DGP A (MCAR +
      exchangeable + rater-varying π), which is the opposite of the
      paper's framing. If no realistic MCAR cell shows AC > IPW SD, the
      Section 5.1 "inefficient" remark needs softening — AC vs FIML is
      the cleaner contrast.
- [x] `04-counts-sampling-misspec` — Counts-FIML under violations of
      Assumption (S). Rater-specific dropout aggregated into counts is
      realistic and out-of-model. Runner/report landed and showed large
      counts-FIML bias under heterogeneous rater MCAR after aggregation,
      while raw FIML/IPW stayed near truth.
- [ ] `05-fiml-sparsity-scaling` — Does FIML bias under MCAR-exchangeable
      shrink at `1/n` as theory predicts, or does the `C^R` parameter-space
      dimension create a slow regime? Starter runner/report landed with
      a low-rep grid. High-rep run still needed to drive Monte Carlo SE
      below `0.001` and settle the manuscript caveat.

### Manuscript

- [ ] Tighten the inference recipe (Section 3.3): the Louis-information
      paragraph is currently a sketch. Either expand with the explicit form
      or move it to the supplement and cite.
- [ ] Add a worked numeric example for inference (Section 7) — at least one
      coefficient with its Wald CI from `vcov()`.
- [ ] Expand citations: `Rubin1976-rb`, `Robins1994-jp`, `Gwet2014-pj`,
      `Conger1980-co` are now in `kappa-missing.bib`. Add U-statistics /
      EM references if a referee asks for them.
- [ ] Decide whether the appendix proof of Proposition
      \ref{prop:ac-consistency} should also state and prove the
      corresponding Gwet limit explicitly (currently in a `\begin{rem}`).
- [ ] Sentence-level pass against `paper/STYLE.md`: no em-dashes, no
      semicolons, no LLM-isms, no antithesis. The current draft mostly
      complies but a final sweep is worth doing.

### Build infrastructure

- [ ] Wire `paper/justfile` recipes that are still TODO: confirm `archive`,
      `tables`, `figures` are functional. `pdf` and `sim` are working.
- [ ] Add `paper/scripts/build_tables.R` (stub for `just tables`) — currently
      the simulation script writes the tables directly.
- [ ] Add `paper/tables/<slug>_stats.tex` macro file with `??` fallbacks for
      any cited numbers in the prose (none currently, but if the worked
      example adds one, this becomes load-bearing).

## Done (2025 rewrite)

- [x] Rewrite scope in `paper/AGENTS.md` to (raw, counts) × (FIML, IPW).
- [x] Rewrite `paper/kappa-missing.tex` around four headline estimators.
      Adds Section 3.2 (FIML, raw), Section 4 (counts + exchangeable raters),
      Section 5 (other estimators, including quadratic special case),
      symbol glossary. Fills the previously empty "Concluding remarks".
- [x] Add missing bibliography entries: `Rubin1976-rb`, `Robins1994-jp`,
      `Gwet2014-pj`, `Conger1980-co`.
- [x] Extend `scripts/simulations_raw_three_estimators.R` with FIML
      (`em_options = list(max_iter = 50000, tol = 1e-7)` to handle harder
      DGPs).
- [x] Counts-FIML model write-up: `paper/dev/notes/counts-fiml-model.md`
      (carried over from earlier work).
- [x] Fix the truncated "We w" sentence in the Introduction.
- [x] Fill Section 7 (R package) and Section 8 (Concluding remarks).
