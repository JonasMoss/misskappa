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

Runner + report landed under `experiments/03-ac-vs-ipw-efficiency/`.
Conger × identity, `n = 1000`, `B = 500`, factorial over exchangeability
× `var(π_j)` × within-subject `corr(M_{i,j}, M_{i,k})`.

**Exchangeable raters (truth = 0.843).** AC and IPW are both unbiased
to MC noise in every cell. AC SD ≤ IPW SD throughout, with the gap
growing in `var(π_j)`: at `var(π_j) = 0`, SDs coincide (IPW reduces to
AC); at the high `π`-variance level, AC SD `≈ 0.010` vs IPW SD
`≈ 0.012`. AC wins MSE in every exchangeable cell. The "AC is
inefficient" claim does not show up here.

**Non-exchangeable raters (truth = 0.695).** AC inherits a large bias
from rater-specific observation rates: `+0.047` (mid `π`-var),
`+0.085` (high `π`-var); IPW residual bias is within MC noise (`|t| < 2`
at `B = 500`). AC's SD is slightly smaller than IPW's (`0.011` vs
`0.016` at high `π`-var, no `corr_M`), but the squared bias dominates.
MSE ratios AC/IPW: `11×` at mid `π`-var, `29×` at high `π`-var. The
exception is `var(π_j) = 0`: AC's pair weighting becomes uniform, so
both estimators are consistent and indistinguishable in MSE.

**Within-subject correlation in `M`.** Across `corr_M ∈ {0, 0.4, 0.8}`,
SDs move by a few thousandths in both directions but the AC vs IPW
ordering never flips. AC's bias under non-exchangeability is set by
the marginal `π_j` schedule, so it barely moves with `corr_M`.

**Implications for the manuscript.** Section 5.1's "AC is inefficient"
framing is backwards in finite samples: at moderate `n` under MCAR
with rater-varying `π`, AC has the lower variance whenever it is
consistent. The clean contrast is "AC is inconsistent under
non-exchangeable raters", not a variance argument. Suggest rewriting
Section 5.1 around the bias mechanism (with a pointer to this
experiment for the simulation evidence) and dropping the inefficiency
claim, or restricting it to "asymptotic efficiency relative to FIML"
where it is on firmer ground.
