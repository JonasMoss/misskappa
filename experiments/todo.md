# experiments TODO

Open experiments backlog. Each is one focused folder under `experiments/`,
following `experiments/AGENTS.md`. Pick one, do it end-to-end (runner +
report), move the entry to "Done" with a one-line summary of the finding.

Experiments are exploration first. They may motivate paper / library
changes, or they may end in the trash. Do not bind manuscript prose to
results until the relevant experiment has landed.

Multiple agents may work in parallel, one experiment per agent. Pick a
slug from "Open", move it to "In progress" with your owner string, and
proceed. Do not touch entries other agents are working on.

## Open

(none)

## In progress

### `05-fiml-sparsity-scaling` — owner: codex

**Question.** Does the smoke-run FIML bias under MCAR-exchangeable
(`~−0.015` at `n = 1000`, `B = 8`) shrink at the `1/n` rate as theory
predicts, or does the `C^R` parameter-space dimension create a slow
regime?

**Approach.** Fix DGP A (or a cleaner exchangeable + MCAR variant).
Sweep `n ∈ {500, 2000, 8000, 32000}` with `B` large enough to drive
the Monte Carlo SE on bias below `0.001`. Plot bias × `n` vs `n`. Two
parameter-space sizes:

- P1: `C = 5`, `R = 6` → `15625` patterns (current paper).
- P2: `C = 3`, `R = 4` → `81` patterns (low-dim sanity).

Also sweep `prune_tol ∈ {1e-12, 1e-9, 1e-6}` at one `n` to confirm the
manuscript hypothesis that pruning shifts `vcov` but not the point
estimate.

**Why it matters.** Decides whether the residual FIML bias seen in the
smoke run is finite-sample MLE bias (clears at full grade) or a
structural concern that warrants a manuscript caveat.

## Done

### `01-coverage-iif-louis`

Runner + report landed under `experiments/01-coverage-iif-louis/`. At
`n ∈ {500, 2000, 8000}`, `B = 50`: IPW Wald coverage is near nominal
under MCAR (DGPs A and B, `cov_95 ∈ [0.92, 0.98]`, `SD(z) ≈ 1`); FIML
Louis SE is systematically too large under DGP A and **does not shrink
at the parametric rate** (mean SE / MC SD ratio grows from 3.6× at
`n = 500` to 9.3× at `n = 8000`), pointing at pseudo-inverse handling
of unidentified θ-directions in the Louis information. Under DGP C the
FIML estimator is biased, not merely conservative, with `cov_95`
collapsing from 0.94 (n = 500) to 0 (n = 8000) and mean `z` growing as
`~√n`: DGP C is technically MNAR for FIML's `X*`-only model because
missingness depends on the latent truth `T_i`, not on observed entries.
Two manuscript implications: (a) Section 3.3 needs a "Louis SE is
conservative in finite samples" caveat or a finite-sample correction;
(b) Section 6 should replace DGP C with a cleanly MAR mechanism (e.g.
missingness conditional on observed ratings of other raters) or rewrite
the prose to flag the MNAR boundary.

### `02-rater-model-sensitivity`

Runner + report landed under
`experiments/02-rater-model-sensitivity/`; the larger `n = 1000`,
`B = 50` run showed unstable Conger identity-loss orderings between
the truth-plus-guess and Dawid-Skene models, so the paper DGP needs a
deliberate choice rather than assuming model invariance.

### `04-counts-sampling-misspec`

Runner + report landed under `experiments/04-counts-sampling-misspec/`.
Counts-FIML is badly biased when heterogeneous rater MCAR data are
aggregated first (Fleiss bias +0.092 at `n = 500`, `B = 30`), while raw
FIML and raw IPW stay near truth. Value-dependent dropout biases all
methods.

### `03-ac-vs-ipw-efficiency`

Runner + report landed under `experiments/03-ac-vs-ipw-efficiency/`. At
`n = 1000`, `B = 500`, IPW wins MSE for Conger's kappa exactly when
raters are non-exchangeable AND `π_j` varies (AC bias up to `+0.085`,
MSE ratio AC/IPW up to `29×`). Under exchangeable raters AC is unbiased
and slightly tighter than IPW, so the "AC is inefficient" framing is
backwards in finite samples; the cleaner manuscript framing is "AC is
inconsistent under non-exchangeability". Within-subject correlation in
`M` shifts variance but does not flip the ordering.
