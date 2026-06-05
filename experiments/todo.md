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

## Paper mapping

The combined `kappa-missing.tex` draft (deleted 2026-06-02, in git
history) was split into three papers (see `papers/split-plan.md`). Each
`report.qmd` should annotate which paper(s) it feeds.

- `01-coverage-iif-louis` — IPW coverage feeds **Paper A**; Louis-SE
  conservatism finding feeds **Paper B** (needs resolution before
  submission).
- `02-rater-model-sensitivity` — DGP-choice question for **Paper A**
  and **Paper B**.
- `03-ac-vs-ipw-efficiency` — **Paper A** centerpiece (AC
  inconsistency under non-exchangeable raters).
- `04-counts-sampling-misspec` — **Paper B** (counts-FIML scope claim
  needs to be sharpened to MCAR + actually exchangeable raters).
- `05-fiml-sparsity-scaling` — **Paper B** (FIML bias scaling /
  manuscript caveat).
- `06-raw-estimator-scaling` — library engineering; feeds no paper
  directly.
- `07-quadratic-edgeworth-coverage` — Edgeworth coverage / quadratic
  variance work feeding **Paper C** (created in parallel by another
  agent; report under `experiments/07-quadratic-edgeworth-coverage/`).
- `08-quadratic-vs-ipw` (open, below) — **Paper C** centerpiece.
- `09-joint-vcov-pilot` — joint inference across non-independent kappa
  estimates on the same data; pilot motivates a per-subject IF /
  paired-spec C++ entry point that would feed **Paper A** (UC4
  Hausman test) and possibly a future standalone note (UC1 pairwise
  homogeneity).
- `10-louis-spectrum` — **Paper B** diagnostic for the rank-truncated
  Louis-information pseudo-inverse.
- `12-clean-mar-dgp` — **Paper B** replacement candidates for the
  inherited fake-MAR DGP C.
- `15-alpha-categorical-smoke` — standalone alpha-missing paper
  smoke/feasibility check for the categorical alpha package path.
- `16-alpha-calibration-sweep` — standalone alpha-missing paper
  starter calibration and feasibility sweep for categorical alpha.

## Open

### `08-quadratic-vs-ipw` — feeds Paper C

**Question.** How does the closed-form quadratic-kappa moment
estimator on pairwise-available data compare to raw IPW
(continuous / quadratic-loss) in finite samples, and is it
MAR-consistent on pairwise-available data?

**Approach.** Continuous bivariate (and ≥ 3-rater) DGPs with known
Lin's CCC. Three missingness mechanisms: (a) MCAR baseline,
(b) MAR conditional on observed ratings of other raters,
(c) MAR conditional on rater identity only. Sweep `n` and missingness
rate. Compare bias, SD, MSE, Wald coverage between
`misskappa::kappa_quadratic()` and `misskappa::kappa(method = "ipw")`
(with quadratic loss). Tie variance estimates back to the U-statistic
asymptotics from Moss (2024).

**Why it matters.** Decides whether Paper C's headline result is
MCAR-only or extends to MAR on pairwise-available data, and settles
the efficiency-gain claim against Paper A's IPW. None of experiments
01–06 are quadratic-focused; without 07 there is no empirical anchor
for Paper C.

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

### `16-alpha-calibration-sweep`

Runner + report under `experiments/16-alpha-calibration-sweep/`.
Computes fixed population alpha truths by one-dimensional quadrature
under latent ordinal item models, then compares alpha-available and
categorical alpha-FIML across complete, MCAR, and anchor-MAR cells.
The committed starter grid spans `C^R ∈ {81, 1024, 15625}`,
`n ∈ {150, 500}`, and `B = 8`; it records bias, MC SD, mean SE,
Wald coverage, timing, observed support, and IF-vcov reconstruction.
Mechanical checks pass, but sparse `5^6` FIML cells show finite-sample
downward bias at low `n`; use the runner for a larger replication and
`n` sweep before manuscript claims.

### `15-alpha-categorical-smoke`

Runner + report under `experiments/15-alpha-categorical-smoke/`. First
mechanical alpha diagnostic for the standalone alpha-missing paper:
complete-data FIML-vs-available equality, small MCAR/MAR categorical
item batteries compared against the paired complete-data alpha, and
starter `C^R` timing rows for saturated categorical EM. This is a smoke
check only; the full bias, SE calibration, coverage, and Louis sweep
remains open.

### `12-clean-mar-dgp`

Runner + report under `experiments/12-clean-mar-dgp/`. Replaces the
old latent-truth DGP C with genuinely MAR mechanisms. The baseline
anchor-MAR design always observes rater 1 and lets raters 2--4 be
observed according to rater 1's observed category; at
`n ∈ {500, 2000, 8000}`, `B = 120`, AC/IPW/Gwet retain positive
Conger bias around `+0.03` to `+0.04`, while FIML is near the full-data
truth. A small candidate search shows more dramatic clean-MAR examples:
sequential four-category MAR gives largest non-FIML bias `+0.144`
versus FIML bias `-0.002` at `n = 4000`, `B = 80`. Recommendation:
use the anchor-MAR design as the main-text replacement and keep the
severe sequential designs as supplement/stress-test evidence.

### `09-joint-vcov-pilot`

Runner + report under `experiments/09-joint-vcov-pilot/`. Joint
inference across non-independent kappa estimates on the same data,
using the analytical IF-based `misskappa::joint_vcov()` helper that
landed in the same commit cycle. `n = 500`, `reps = 200`. Whole sweep
runs in ~2 s.

Four sub-questions:

- **UC8** (Fleiss-vs-BP same-fit contrast from the existing 3×3
  vcov): reject 0/200 under uniform marginals (population contrast = 0),
  200/200 under skewed. Methodological remark only.
- **UC4** (Hausman test, AC vs IPW on same incomplete data): main
  positive result. Reject 0.045 under exchangeable + MCAR null
  (~ nominal), 1.0 under the experiment-03 non-exchangeable cell.
  Analytical SE matches MC SD of the contrast to ~2%. Per-dataset
  diagnostic for AC bias under non-exchangeability; feeds **Paper A**.
- **UC1** (all 15 pairwise Cohen kappas + joint Wald homogeneity
  test, chi-sq df=14): power 200/200 against one miscalibrated rater;
  size 0.16 at `n = 500` (chi-sq approximation residual; drops to
  ~0.075 at `n = 2000`).
- **UC2** (three weighted kappas on same data): identity vs
  quadratic statistically distinct in 87% of replicates at `n = 500`
  despite joint correlation ≈ 0.94 between linear and quadratic.
  Useful weight-sensitivity demo.

Pilot motivated the influence-function exposure that now ships in
`misskappa::influence()` and `misskappa::joint_vcov()` for the
categorical raw estimators.

### `10-louis-spectrum`

Runner under `experiments/10-louis-spectrum/`. One DGP-A draw at
`n = 8000` decomposes Conger's FIML Louis variance by eigendirection of
the reduced observed-information matrix. With the untruncated inverse,
effectively all variance comes from relative eigenvalues below `1e-3`;
with `info_rcond = 5e-5`, the diagnostic variance matches the package
`vcov()` value (`1.26e-5`) and removes the near-kernel blow-up.

### `01-coverage-iif-louis`

Runner + report landed under `experiments/01-coverage-iif-louis/`. At
`n ∈ {500, 2000, 8000}`, `B = 50`: IPW Wald coverage is near nominal
under MCAR (DGPs A and B, `cov_95 ∈ [0.92, 0.98]`, `SD(z) ≈ 1`). The
FIML Louis SE blow-up was resolved by the rank-truncated inverse
(`info_rcond = 5e-5`): DGP-A mean SE / MC SD ratios are now
`1.22, 1.09, 1.12` rather than the old `3.6, 7.9, 9.3`, and DGP-B
ratios are `1.00, 1.03, 1.10`. Remaining subnominal coverage in some
A/B cells is driven by finite-sample FIML point-estimate bias. Under
DGP C the FIML estimator is biased, not merely conservative, with
`cov_95` collapsing from 0.94 (n = 500) to 0 (n = 8000) and mean `z`
growing as `~√n`: DGP C is technically MNAR for FIML's `X*`-only model
because missingness depends on the latent truth `T_i`, not on observed
entries. Section 6 should replace DGP C with a cleanly MAR mechanism
or rewrite the prose to flag the MNAR boundary.

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
