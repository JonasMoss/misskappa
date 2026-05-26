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

### `05-fiml-sparsity-scaling`

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

## In progress

### `03-ac-vs-ipw-efficiency` — owner: claude

**Question.** Under what conditions on `π`-variability and rater-
exchangeability does IPW's MSE actually beat AC's — i.e. when is
the "AC is inefficient" claim visible in finite samples?

**Approach.** Hold the rater model fixed (latent truth + guess). Sweep
over (i) exchangeability vs not, (ii) variance of the `π_j` vector,
(iii) within-pair correlation of `M_{i,j} M_{i,k}` (introduce a per-
subject latent missingness factor). Report bias, SD, and MSE for AC and
IPW at one moderate `n`.

**Why it matters.** The smoke run of the paper sim has IPW SD `> AC SD`
under MCAR + exchangeable + rater-varying `π` (DGP A), the opposite of
the paper's framing. If no realistic MCAR cell shows `SD(IPW) < SD(AC)`,
the manuscript Section 5.1 "inefficient" prose should soften — AC
inconsistent under non-exchangeability is the cleaner contrast.

### `01-coverage-iif-louis` — owner: claude

**Question.** Do IPW influence-function Wald CIs and FIML Louis Wald CIs
for the weighted kappas achieve nominal coverage at moderate `n`?

**Approach.** Reuse DGPs A / B / C from
`paper/scripts/simulations_raw_three_estimators.R`. For each
DGP × estimator × `n`, run `B` replicates; on each, fit
`kappa(x, method = m)`, extract estimate and `sqrt(diag(vcov))`, form
the 95% Wald CI, check inclusion of the truth. Report empirical
coverage (90 / 95 / 99), mean CI width, and the standardised residual
`(κ − truth) / SE` (which should be approximately `N(0, 1)`).

Conger × identity is the headline; extend to Fleiss / Brennan-Prediger
only if the answer is ambiguous. Sample sizes span small / moderate /
large.

**Why it matters.** The first referee question. If FIML undercovers
visibly at moderate `n`, Section 3.3 of the manuscript needs a finite-
sample correction or an explicit caveat.

## Done

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
