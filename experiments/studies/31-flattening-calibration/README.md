# 31 — Calibration of the cat_fiml flattening constant

**Question.** After demoting the gradient identification guard to a warning
and adding Dirichlet flattening to the raw cat_fiml EM
(`em_options = list(flatten = c)`, c = total pseudo-mass over the C^R
table), what value of c should the package recommend? The note
`papers/00_wip/kappa-missing/dev/notes/cat-fiml-mle-face-flattening.md`
predicts: any fixed c is asymptotically harmless (shift ≈ c·κ/n), the
window c ∈ [0.01, 1] should be statistically indistinguishable, and the
deliverable is a confirmation table, not a tuned constant.

**Design.** Latent-class DGP with exact truth (uniform latent class,
symmetric confusion with accuracy `acc`; all three coefficients share
kappa = (po − 1/C)/(1 − 1/C)). Grid:

- DGPs: lc3x3 (C=3, R=3, acc=.75), lc4x4 (C=4, R=4, acc=.70),
  lc5x4 (C=5, R=4, acc=.70) — the last two are the regimes where the old
  guard misfired.
- Mechanisms: `pairs` (planned missingness, one random rater pair per
  subject, MCAR) and `anchor` (MAR: rater 1 always observed, others
  observed w.p. 0.9/0.35 depending on rater-1 category).
- n ∈ {20, 40, 100}; flatten c ∈ {0, 0.01, 0.1, 0.5, 1}; 1000 replicates,
  common random numbers across c.

Metrics per cell: failure rate (now only design-guard failures — pair never
co-observed), bias, SD, RMSE, SE/SD, 95% Wald coverage, mean and p90 of the
`null_frac` diagnostic.

**Run:** `Rscript run_calibration.R [reps] [cores]` →
`results/calibration-raw.csv` + `results/calibration-summary.csv`.

**Findings (1000 reps, 2026-06-10).**

1. **The c ∈ [0.01, 0.1] window is statistically free.** RMSE ratio vs
   strict ML: median 1.000–1.001, max 1.005 across all 18 cells. Bias shift
   at n = 20 is ≤ 0.003 for c = 0.1. At c = 1 the attenuation becomes
   visible (bias shift up to −0.023 at n = 20, lc3x3; RMSE up to +2.9%) —
   matching the w·κ̂ = c/(n+c)·κ̂ prediction. **Recommend c = 0.1**; avoid
   c ≥ 0.5 as a default.
2. **Failures are now design-guard only.** With the face-tolerance stopping
   rule (see below), flatten > 0 has zero `not_converged` failures; the only
   errors are genuine "pair never co-observed" design failures (max 6% at
   n = 20, lc4x4 pairs). Strict ML (c = 0) retains a small residual of its
   own slow-EM non-convergences.
3. **Flattening inflates the SE (conservative).** se/sd jumps from ~0.96
   (strict) to ~1.47 (any c > 0); coverage rises from median 0.88 to 0.97.
   Mechanism (verified on the Louis spectrum): flattening retains the full
   data-supported table (206 vs 33 cells in a probe), tripling the count of
   eigenvalues just above the `info_rcond` cut; their reciprocals inflate
   the delta-method variance ~50%. Not an error — over-coverage — but a
   refinement target: prune harder (strict-ML-style support) for the
   variance step only, or retune `info_rcond` under flattening.
4. **Remaining undercoverage is FIML small-n point bias, not SEs.** The
   worst cells (cover ~0.70 at n = 20, C ≥ 4 pairs) are driven by the known
   negative finite-sample MLE bias (−0.1 to −0.19), unchanged across c.
   That is the separate "classify residual FIML bias" thread.

**Implementation note.** The first run exposed a convergence-cap bug: EM
drift along the likelihood-flat face is slower than the n/(n+c) linear
model in big sparse tables (up to 8% `not_converged` at C=5, c=0.1). Fixed
in `run_em_iterations`: iteration cap scaled by 300·(n+c)/c (bounded 1e6)
and the tolerance floored at `1e-4 · c/(n+c)` — the analytic-center
position only needs resolving to ~1e-4 because identified functionals are
flat along the face.
