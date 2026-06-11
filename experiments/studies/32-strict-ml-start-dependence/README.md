# Study 32: start-dependence of strict-ML cat_fiml

**Question.** With `flatten = 0` (strict ML), how much do the kappa point
estimate, its SE, and the `null_frac` diagnostic move per dataset when the
EM start changes? This is the empirical close-out of "should flattening be
the default": the design guard already pins kappa exactly on the
population MLE face (note §3 of
`papers/00_wip/kappa-missing/dev/notes/cat-fiml-mle-face-flattening.md`),
and study 31 showed face-point selection is free in aggregate (RMSE ratio
≤ 1.005 vs the flattened analytic center) while flattening breaks SE
calibration (se/sd 0.96 → ~1.5). What study 30/31 did not measure is the
per-dataset spread across starts, including the SE.

**Start perturbation.** `start_alpha ∈ {0.01, 0.1, 1, 10}` — the only
start perturbation the C++ supports: `initialise_theta()` is deterministic
(`theta_j ∝ start_alpha + Σ_g n_subs_g / n_comps_g`), no RNG anywhere in
the EM. Wider face exploration would need a jitter option (study 30 used
ad-hoc multi-start at one design point and found median kappa width
2.7e-3 at n = 40, 1-in-5 seeds ~9.8e-2).

**Design.** Study-31 DGPs (lc3x3 / lc4x4 / lc5x4) × {pairs-MCAR,
anchor-MAR} × n ∈ {20, 40} (face width is ~1e-12 by n = 160; study 30
script 04), 1000 reps, common random numbers, seed base `32000000 + k`.

**Metrics** (per cell and coefficient, over datasets where all four starts
converged): median / p90 / max of `spread_est / se(start = 0.1)`;
`frac_est_moves` = P(spread_est > 0.1·se); `frac_se_moves` =
P(max se / min se − 1 > 0.10); start-dependent-convergence rate; and
whether `null_frac > 0.01` at the reference start flags exactly the
datasets the start can move (`flag_rate_movers` vs `flag_rate_all`).

**Pre-registered reading.** If `frac_est_moves` and `frac_se_moves` are
negligible (< 1–2%) in every cell, start-dependence does not justify
changing the package default away from strict ML, and flattening stays a
documented opt-in uniqueness device (with its SE-conservatism caveat). If
`frac_se_moves` is material anywhere, reopen the parked SE work
(variance-support pruning / penalized-bread sandwich, note §7b) before
relying on strict-ML SEs in study 29.

Run: `Rscript run_experiment.R 1000 10` (smoke: `Rscript run_experiment.R 5 2`).

## Findings (1000 reps, 2026-06-10)

1. **The point estimate is start-independent in practice.**
   `frac_est_moves = 0` in every cell and coefficient: across 12 cells ×
   1000 datasets, no estimate ever moved by more than 0.1·SE across the
   four starts. The largest spread observed anywhere was 0.0098·SE
   (median spreads are 1e-8–1e-5·SE). The §4 face width exists, but the
   deterministic start lands at an equivalent point regardless of
   `start_alpha`.
2. **SE start-dependence is rare and flagged.** The SE moved by > 10%
   across starts in 0–3.3% of datasets (worst: lc5x4/lc4x4 anchor-MAR at
   n = 20: 3.3% / 2.7%; pairs-MCAR ≤ 0.25%; everything at n = 40 ≤ 1.1%).
   In the C ≥ 4 cells every such mover had `null_frac > 0.01` at the
   reference start (`flag_rate_movers = 1`), so the shipped diagnostic
   flags exactly the datasets where the rank-truncated SE is fragile.
3. **Convergence never depends on the start** (`start_dep_conv = 0`
   everywhere). All-start failures track the design guard (pair never
   co-observed), reaching ~19% at n = 20 under the pairs design with
   C ≥ 4 — genuine non-identification, not numerics.
4. **Bug found in study 31's summary:** its `fail_rate` averaged the
   failure flag over rows, but failed fits contribute 1 row vs 3 (one
   per coefficient) for successes, understating failure rates ~3×
   (row share r implies true rate p = 3r/(1+2r); reported 6.8% at
   lc4x4/pairs/20 is really ~18%, matching this study). Fixed in
   study 31's `summarize.R`, summary regenerated; bias/SD/RMSE/coverage
   columns were computed per coefficient on successful fits and are
   unaffected.

**Verdict (per the pre-registered reading).** Start-dependence does not
justify flattening as the default: the estimate never moves, convergence
never changes, and the only start effect is a > 10% SE wobble in ≤ 3.3%
of the sparsest MAR cells — datasets the `null_frac` warning already
flags. Flattening would trade that localized wobble for a uniform ~50%
SE inflation (study 31). Strict ML stays the default and the
recommendation; flattening remains a documented opt-in uniqueness
device. See note §7b.
