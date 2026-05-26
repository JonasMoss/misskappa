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

### `02-rater-model-sensitivity`

**Question.** Do the A / B / C bias and efficiency orderings reported in
`paper/scripts/simulations_raw_three_estimators.R` survive a Dawid-Skene
rater model?

**Approach.** Replicate the design of `simulations_raw_three_estimators.R`,
replacing the latent-truth-plus-guess rating model with per-rater
confusion matrices `θ^j_{c,c'} = P(X_j = c' | T = c)`. Construct three
DGPs in the D-S framework that match the spirit of A / B / C:

- A': all raters share a single confusion matrix (exchangeable). MCAR
  rater-specific dropout, same `π` as paper DGP A.
- B': raters have heterogeneous confusion matrices (non-exchangeable).
  MCAR rater-specific dropout, same `π` as paper DGP B.
- C': raters share a confusion matrix; missingness depends on the
  observed rating of the previous rater (a clean MAR mechanism, not
  truth-dependent).

Fit AC / IPW / FIML / Gwet on each cell. Report bias and SD side by side
with the paper sims. Drives the DGP choice for the final manuscript.

**Why it matters.** Current paper sims use a constrained rater model
(identity confusion plus uniform mixture). If orderings change under
D-S, the paper should run D-S; if they hold, the simpler model stays.

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

### `04-counts-sampling-misspec` — owner: codex

**Question.** How badly does counts-FIML break when Assumption (S)
(uniform hypergeometric subsampling) fails?

**Approach.** Generate raw rater-identified data with rater-specific
dropout, aggregate to counts, run `kappa_counts(..., method = "fiml")`.
Two failure flavours:

- M1: rater-specific MCAR dropout (`π_j` varies); after aggregation the
  observed counts are biased toward the more-observed raters' marginal.
- M2: value-dependent dropout (low-category raters drop more often).
  NMAR for any of the estimators.

Compare counts-FIML bias to (a) rater-identified FIML on the same
underlying data and (b) raw-data IPW. Establishes the scope claim in
Section 4 of the manuscript.

**Why it matters.** Counts-FIML's Assumption (S) is plausible only under
MCAR with exchangeable raters. The paper says this; the simulation
should show how badly it fails when (S) is violated.

## Done

(none yet)
