# Study 33 — nt-FIML degenerate-Σ rescue: does it buy trustworthy inference?

**Date:** 2026-06-11. **Verdict: no — do not enable the rescue. nt-FIML's n≥40
limit stands, with the mechanism corrected.**

## Background / what this corrects

The small-n nt-FIML "Numerical error" failures in study 29 were attributed to the
sandwich solve `solve(H)` on the observed information. That was the *old pure-R*
path. In the live C++ backend (`ef961ac`) the failures are upstream: the EM
E-step inverts each pattern's observed-block covariance `Σ_oo`, and at small n
with missingness the fitted saturated covariance goes **singular**. Instrumented
sweep (p∈{4,5,6} × n∈{10,20,40}, 5 cats, 30% MCAR, 1800 cells): the H-solve
**never** fails; **every** failure is `Σ_oo`.

A singular Σ is *not* bad for the estimand — quadratic kappa is a **linear
functional of Σ** (`sum(Σ)`, `tr Σ`, mean spread; no `Σ⁻¹`), so it is well
defined on the boundary. The singularity only breaks the EM's matrix inversions.

## The rescue (experimental, gated OFF by default)

`em_options$soo_rcond > 0` routes the `Σ_oo` inversion through the rcond-truncated
PSD pseudo-inverse (minimum-norm conditional expectation), letting the EM
converge through the boundary. Estimand-preserving by construction. A degenerate
Σ then yields a singular H, which the companion sandwich pseudo-inverse keeps
finite, with `null_frac` reporting the truncated gradient mass.

Implementation parked on branch **`nt-fiml-soo-rescue`** (reverted from master).
`smoke.R` requires that branch to run.

## The smoke (smoke.R → smoke_raw.rds)

DGP: exchangeable latent-normal subject effect + per-rater bias (so Conger≠Fleiss),
5 raters × 5 categories, ρ=0.6. Truth from large-sample (μ,Σ). 1200 reps ×
n∈{10,20,40} × miss∈{0,0.2}. Stratify by `rcond(Σ) < 1e-8` = degenerate.

**Recovery:** rescue eliminates all failures — n=10/miss=0.2 went **49% → 0%**,
n=20 **7% → 0%**. Degeneracy is missingness-driven and **gone by n=40**.

**SE calibration (Conger; Fleiss identical in shape):**

```
 n  miss  stratum      nrep   bias  mc_sd  mean_se  se/sd  cov95  med_null_frac
10  0.2   clean         209  -0.056  0.176   0.136   0.78   0.78    0.000
10  0.2   degenerate    991  -0.102  0.178   0.066   0.37   0.37    0.539   ← rescued
20  0.2   clean        1087  -0.035  0.116   0.107   0.92   0.89    0.000
20  0.2   degenerate    113  -0.072  0.121   0.039   0.32   0.32    0.508   ← rescued
40  0.2   clean        1200  -0.013  0.080   0.077   0.96   0.92    0.000
                        (no degenerate fits remain at n=40)
```

Rescued degenerate fits are **more biased** (−0.10 vs −0.03) and their SEs are
**optimistic** — se/sd ≈ 0.35, so **coverage collapses to ~0.35** vs the clean
stratum's ~0.85–0.92. Shipping the rescue silently would report confident,
wrong intervals — worse than erroring. `null_frac` separates the strata cleanly
(0.000 clean vs ~0.5 degenerate), so the only defensible rescue shape is
"rescue the point, blank/flag the SE when null_frac is high" — but the point is
biased too, so erroring remains the more honest default.

## Decision

- **Package: no change.** Rescue + sandwich pinv + null_frac reverted to HEAD
  (parked on `nt-fiml-soo-rescue`). The H-pinv fires 0/1800 cells when the
  rescue is off — dead code by default — and only makes sense bundled with the
  rescue, which we are not shipping.
- **study-29 / manuscript:** keep the n≥40 guidance; correct the reason from
  "singular sandwich solve" to "singular fitted saturated covariance."
- Caveat: one DGP, one ρ. The mechanism (optimistic SE on a near-singular
  information) is robust enough that a broader DGP×ρ sweep was judged not worth
  it before concluding.

## Runner / local timing update (2026-06-12)

`run_experiment.R` now provides the maintained entry point:

```sh
Rscript experiments/studies/33-nt-fiml-degenerate-sigma/run_experiment.R --smoke --load-all
Rscript experiments/studies/33-nt-fiml-degenerate-sigma/run_experiment.R --screen --load-all --reps 5
Rscript experiments/studies/33-nt-fiml-degenerate-sigma/run_experiment.R --full --load-all
```

The compatibility `smoke.R` delegates to `run_experiment.R --smoke --load-all`.
Outputs are rectangular CSVs under `results/<mode>/`: `raw.csv`,
`recovery.csv`, `timing.csv`, `calibration.csv`, `truth.csv`, and
`metadata.csv`.

Timed on the parked `nt-fiml-soo-rescue` branch in a temporary worktree:

- cold `--smoke` (3 reps, n = 10/20, miss = 0.2, includes first `load_all`
  compile/load): 1:52.64 wall, 110.76 user, 8.62 sys.
- warm `--smoke`: 3.40 s wall, 3.20 user, 0.19 sys.
- warm `--screen --reps 5 --truth-n 50000` over the full cell grid
  n = 10/20/40 x miss = 0/0.2: 5.18 s wall, 5.00 user, 0.18 sys.

Median rescued fit times from the 5-rep screen were 0.010, 0.009, and 0.010 s
for complete-data n = 10, 20, and 40; and 0.146, 0.105, and 0.101 s for 20%
MCAR n = 10, 20, and 40. The original 1200-rep study is therefore a local run:
roughly 15--25 minutes sequential after the package is warm, not a week-scale
job. A 4-core desktop can shard this manually or via a small future runner
extension if needed; Modal is not necessary for Study 33.
