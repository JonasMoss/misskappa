# Three-paper split plan

Decision (2026-05-27): the combined `papers/combined/kappa-missing.tex` manuscript
is being broken into three independent papers, each with its own
audience and venue. Single source of truth for the split lives in this
file. See also `/home/jmoss/.claude/plans/i-think-we-might-tender-lerdorf.md`
for the discussion that led here.

The combined draft remains the canonical reference for written content
until each spinoff paper has migrated its sections; do not delete
content from `papers/combined/kappa-missing.tex` until the corresponding
`papers/<slug>/` directory compiles standalone.

## Papers

| Slug | Working title | Target venue | Profile | Pages |
|------|---------------|--------------|---------|-------|
| `papers/ipw` | Inverse-probability weighting for kappa with missing raters | Psychometrika | Math (lower-quartile) | ~10–12 |
| `papers/fiml` | Semiparametric-efficient estimation of kappa under MAR | Psychometrika or Biometrika | Math | ~14–18 |
| `papers/quadratic` | Closed-form quadratic kappa and Lin's concordance with missing pairwise data | Biometrics or similar | Math | ~12–14 |

## Section mapping from the combined draft

Source sections are identified by their slot in the current
`papers/combined/kappa-missing.tex` (rewrite of 2026-05-26).

### Paper A (IPW)

- Sec 1 — Introduction (rewrite, trim to IPW scope)
- Sec 2 — Setup (full notation, MCAR/MAR, glossary; self-contained)
- Sec 3.1 — Raw IPW (Hájek-type, Prop 1)
- Sec 3.3 — Inference for IPW (influence-function / V-statistic delta)
- Sec 4.1–4.2 — Counts model + pairwise plug-in (Prop 3)
- Sec 5.1 — Available-case (Prop 5 + appendix proof)
- Sec 5.2 — Gwet's reweighting
- Sec 6 — Simulations (IPW / AC / Gwet rows only; drop FIML)
- Sec 7 — R package, IPW paths only
- Sec 8 — Concluding remarks (IPW-focused)

Centerpiece empirical result: experiment 03 (AC vs IPW). Reframe AC
as **inconsistent under non-exchangeable raters** rather than
"inefficient".

### Paper B (FIML / efficiency)

- Trimmed Setup (cite Paper A for the bulk of notation)
- Sec 3.2 — Raw FIML (Prop 2, semiparametric efficiency under MAR)
- Sec 3.3 — Louis information for FIML with finite-sample caveat
  (Exp 01) or a finite-sample correction
- Sec 4.3 — Counts FIML (Prop 4)
- Sec 4.4 — Moore–Penrose inference + **Remark 1 (counts
  identifiability of κ despite non-identified θ)**
- Sec 6 — Simulations (FIML rows; possibly new sims)
- Sec 7 — R package, FIML paths
- Sec 8 — Concluding remarks (efficiency-focused)

Key empirical inputs: Exp 01 (Louis SE conservatism — needs
resolution), Exp 04 (counts-FIML misspecification), Exp 05 (FIML
sparsity, in progress).

Open work before submission: resolve the Louis-SE finite-sample issue
(at minimum a documented caveat; ideally a fix).

### Paper C (quadratic / Lin's CCC)

- Setup focused on continuous ratings and quadratic loss
- Lin's CCC framing: quadratic-weight κ_X ≡ Lin's CCC for continuous
  bivariate data; survey CCC-with-missing-data literature (sparse)
- Closed-form moment estimator on pairwise-available data (seeded by
  current Sec 5.3, ~17 lines) — full derivation, asymptotic
  distribution, comparison to IPW
- Extension to ≥ 2 raters (Conger / Fleiss quadratic variants)
- MCAR vs MAR consistency on pairwise-available data — resolve
  during writing; foreground if MAR-consistent, otherwise frame as
  MCAR-efficient closed form
- New simulations comparing closed-form vs raw IPW on continuous /
  quadratic-loss data (currently no quadratic-focused experiment
  exists; a `07-quadratic-*` slug is queued)

## Shared material policy

- **Symbol glossary, notation, MCAR/MAR setup.** Each paper carries
  its own. Paper A's is full; Papers B and C cite Paper A and restate
  only what they need.
- **Bibliography.** Each paper folder owns its own `.bib`, per
  `papers/combined/STYLE.md` ("One .bib per paper folder. Do not share .bib
  files across papers; copy the entries you need."). Initial copy
  comes from `papers/combined/kappa-missing.bib`.
- **R package material.** Each paper points at `misskappa` and
  documents only the entry points it uses.

## Repository layout

All four papers live under `papers/`:

```
papers/combined/        # source-of-truth combined draft; do not delete
papers/ipw/
papers/fiml/
papers/quadratic/
experiments/            # unchanged; reports annotate which paper(s) they feed
```

Each `papers/<slug>/` directory mirrors `papers/combined/` layout per `papers/combined/AGENTS.md`:

```
papers/<slug>/
  AGENTS.md             # cloned from papers/combined/AGENTS.md, paper-specific scope
  STYLE.md              # cloned from papers/combined/STYLE.md (shared style)
  justfile              # cloned, slug substituted
  kappa-missing-?.tex   # paper-specific tex
  kappa-missing-?.bib   # paper-specific bib (subset of combined .bib)
  dev/                  # paper-specific notes + todo
  figures/              # paper-specific figures
  tables/               # paper-specific tables
  scripts/              # paper-specific simulation/build scripts
  results/              # curated CSVs
  supplement/           # online supplement source
```

## Sequencing

1. **Paper A first.** Cleanest cut; lands quickly in Psychometrika;
   builds citation surface for B and C.
2. **Paper B second**, in parallel with A in review. Substantial —
   Louis-SE work is non-trivial.
3. **Paper C last.** Needs new sims and Lin's CCC literature survey.
   Build framing while Paper B is in late drafting.

## Migration policy

- Copy (don't move) sections from `papers/combined/kappa-missing.tex` into the
  spinoff directories. Mark migrated sections in `papers/combined/kappa-missing.tex`
  with a `% MIGRATED to papers/<slug>/` comment so source-of-truth status
  remains explicit.
- Once a paper compiles standalone and has been read end-to-end, retire
  the corresponding sections from `papers/combined/kappa-missing.tex`.
- `papers/combined/dev/todo.md` items get partitioned across the three
  `papers/<slug>/dev/todo.md` files as their target paper becomes
  unambiguous; items that touch all three remain in `papers/combined/dev/`.

## Cross-paper experiments

`experiments/` continues to live at the repository root. Each
`report.qmd` should annotate which paper(s) it feeds. Existing mapping:

- Exp 01 — IPW coverage feeds Paper A; Louis-SE finding feeds Paper B
- Exp 02 — rater-model sensitivity feeds Paper A and Paper B
- Exp 03 — AC vs IPW efficiency feeds Paper A (centerpiece)
- Exp 04 — counts-FIML misspecification feeds Paper B (centerpiece)
- Exp 05 — FIML sparsity scaling feeds Paper B
- Exp 06 — raw estimator scaling feeds engineering / library, not a paper
- Exp 07 — `07-quadratic-*` (to be created) feeds Paper C
