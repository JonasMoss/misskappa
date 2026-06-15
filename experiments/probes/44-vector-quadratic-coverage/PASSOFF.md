# Probe 44 — passoff note (next session: framing + sim mods)

**Status (2026-06-15).** All committed to `master`:
- `b3f8874` — `ranges=` fixed-range Brennan-Prediger for the vector quadratic estimators (C++→rcpp→R, tests, NEWS). Full testthat 536/0, C++ unit tests pass.
- `f18f9fa` — probe 44 itself: screen coverage (`parity_and_coverage.R`) + three-regime stress (`stress.R`), NOTES.md, summary RDS/CSV.
Memory: `project_probe44_vector_performance`. Nothing pushed to remote.

This note is the agenda for a **cleared session**: decide the **framing** and the **sim mods** before any of this becomes paper-grade. Read `NOTES.md` first (the full results); this is the "what's undecided" layer on top.

---

## What we now know (one-screen recap)

Estimators: `pairwise` (MCAR moment, Conger/Fleiss), `nt_fiml` (saturated), `nt_fiml_cs` (block-CS `Σ = J_R⊗Ψ + I_R⊗E`), on subjects×raters×features arrays, squared vector loss `(x−y)'W(x−y)`.

1. **MCAR:** all three coincide on Conger/Fleiss (constraint-invariant contractions `tB,tT,tG`); `nt_fiml_cs` is consistent **even when block-CS is false**, stress-tested to extreme misspecification. No efficiency gain over saturated for these coefficients — value is robustness, not variance.
2. **MAR:** `nt_fiml_cs` is **genuinely inconsistent**, severity graded in misspecification strength `s` (at `s=1`, Conger bias +0.022, coverage 0.74→0.34 as n 400→1500). Saturated `nt_fiml` is the consistent fix **but fails to converge** under strong covariance + missingness (singular per-pattern blocks). → a **consistency-vs-stability tradeoff**.
3. **Robustness (conv-rate):** in sparse/small-R designs `nt_fiml_cs` is the lone survivor (conv 0.95–1.0) while saturated collapses (0.00–0.40) — the connected-not-complete payoff.
4. **Brennan-Prediger:** broken with observed range (cov 0); **fixed** with `ranges=` (cov 0.93–0.97). Decide the default/UX (below).

---

## Open framing questions (the real work)

1. **Headline or appendix?** The MAR inconsistency argues for an honest "online-appendix comparator," but the MCAR consistency + robustness + connected-graph identifiability is a real selling point. Mirror the scalar decision in `papers/00_wip/quadratic/dev/notes/cs-reframing-memo.md` (and the structured-ICC framing, memory `project_ten_hove_icc`). **Decide: is vector block-CS a recommended estimator, or a studied-but-cautioned one?**
2. **How to present the tradeoff.** MCAR → recommend `nt_fiml_cs`. MAR → `nt_fiml_cs` only if block-CS≈right (graded bias), else saturated `nt_fiml` (consistent but may not converge). Is this a table, a decision flowchart, a proposition? What's the crisp rule a reader leaves with?
3. **Estimand / literature anchor.** Vector quadratic kappa == generalized concordance correlation (Lin / King–Chinchilli / Carrasco–Jover, already cited in `kappa_vector_quadratic.R`); the CS route is the repeated-measures-CCC variance-component model. Decide whether to frame the whole thing as "CCC under missing data" and connect to that world (relates to memory `project_concordance_study` / study 28).
4. **Where does it live?** Quadratic paper appendix vs a dedicated vector/CCC note vs alongside `alpha-missing`. Affects scope and which sims are "paper figures."

---

## Potential sim mods (before paper-grade)

- **Cleaner adversarial knob.** The two-camp `camps` is PSD-floored at `s≥2` (near-singular DGP), which confounds "CS bias" with "saturated fails because the truth is degenerate." Redesign a strength knob that stays **well-conditioned** at high `s` (e.g. tune off-diagonal camp coupling only, keep diagonal blocks PD) so the `nt_fiml_cs` bias curve is isolated from DGP degeneracy. Anchor headline claims on `s∈{0.5,1}` (already clean).
- **Is saturated's non-convergence fixable?** Probe whether a ridge / rank-aware fallback in the vector saturated EM (`kappa_vector_quadratic.R` docstring flags there's none) recovers the failed cells. If saturated could be made robust, the MAR tradeoff narrative weakens (saturated becomes the all-around answer). Relates to scalar `project_nt_fiml_degenerate_sigma` (that rescue went optimistic on SEs — watch for the same).
- **BP default/UX.** Now that `ranges=` exists: should `kappa()` **warn** when vector BP is read off the observed range (nudge users to set `ranges`)? Or document-only? Decide before exposing BP CIs in examples.
- **Breadth not yet run.** `parity_and_coverage.R --targeted` (t6 non-Gaussian, d∈{1,2,3}, coupling W, larger n, B=1000) is coded but unrun — needed for the paper simulation table and to confirm FIML-SE under-coverage under heavy tails. Stays local `mclapply`; graduate to `studies/34` if it grows.
- **BP truth.** Screen run used a large-n MC BP truth (range-dependent, hence the messy rows). With `ranges=` the BP truth is analytic (regime 3 does this for `blockcs`); redo any BP coverage cells with fixed range + analytic truth, drop the MC-truth crutch.
- **Intervals / n_eff.** Screen uses Wald + Fisher only. Study 29 also reports t-critical and arcsine intervals and an `n_eff` convention — decide whether to match for consistency across the paper family.
- **MCAR efficiency, once more.** We found `nt_fiml == nt_fiml_cs` (no variance gain) for Conger/Fleiss. Confirm there's no regime (very heavy MCAR, where the contraction estimates actually diverge) where CS wins — one targeted check would settle it.

---

## Re-run / artifacts

```
Rscript experiments/probes/44-vector-quadratic-coverage/parity_and_coverage.R           # screen (~45s)
Rscript experiments/probes/44-vector-quadratic-coverage/parity_and_coverage.R --targeted # breadth (heavy)
Rscript experiments/probes/44-vector-quadratic-coverage/stress.R                         # all 3 regimes (~6min)
Rscript experiments/probes/44-vector-quadratic-coverage/stress.R 1 3                      # pick regimes
```
Outputs: `coverage_summary.{rds,csv}`, `stress_regime{1,2,3}.rds` (committed); `replicates.csv` (gitignored, regenerable).

Key code: `r-package/R/kappa_vector_quadratic.R` (backend + `ranges`), `r-package/R/kappa.R` (`.kappa_vector_dispatch`), `src/estimate_normal_fiml.cpp` (`vector_bp_chance`, `estimate_vector_quadratic_normal_fiml`). Related memories: `project_ten_hove_icc`, `project_ten_hove_variance_component`, `project_concordance_study`, `project_quadratic_nt_fiml`.
