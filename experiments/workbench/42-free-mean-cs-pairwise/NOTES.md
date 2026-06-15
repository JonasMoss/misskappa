# Probe 42 — free-mean CS pairwise moment estimator ("pairwise + variance_components")

**Question.** Is the rater-centered free-mean compound-symmetry (CS) pairwise
moment estimator — the moment foil of the free-rater-mean CS FIML, estimand =
Conger's quadratic kappa — a *principled generalisation* of the (equal-mean)
pairwise, and does it work? Background:
`papers/00_wip/quadratic/dev/notes/cs-fixed-rater-pairwise.md`,
`dev/notes/counts-exchangeable-cs-direct.md`.

**Answer.** Yes on all counts. It nests the equal-mean pairwise by mean-model
restriction, shares the free-mean CS FIML's probability limit, the `t3`
mean-spread bias correction is necessary and correct, and it is MCAR-consistent
with calibrated (bootstrap) SE / MAR-biased — exactly the moment-estimator
profile, with FIML the robust route.

## Files
- `free_cs_pairwise.R` — `cs_free_pairwise()` (recipe: rater-center, pooled
  within-subject pairwise moments on residuals, `t3 = sum(mu_j-bar)^2 - tr(AV)`),
  `equal_mean=TRUE` collapses to the one-way pairwise. `..._boot()` = subject
  bootstrap SE (probe-grade).
- `verify.R`, `calibration.rds`.

## Recipe (the estimator)
    mu_j = available column mean ; r_ij = Y_ij - mu_j
    e    = sum_i sum_{j<k} (r_ij - r_ik)^2 / (2 sum_i C(m_i,2))
    psi  = sum_i sum_{j<k} r_ij r_ik / sum_i C(m_i,2)
    t3   = sum_j (mu_j - mubar)^2 - tr(A V),  A = I - 11'/R,
           V_jj=(psi+e)/n_j, V_jk = psi n_jk/(n_j n_k)
    Conger = psi/(psi + e + t3/(R-1)),  Fleiss = (R(R-1)psi - t3)/((R-1)(R(psi+e)+t3))

## Results (R=5, psi=0.8, e=1.2; equal-mean pop Conger=Fleiss=0.40, free-mean Conger=0.305 Fleiss=0.270)

**A. Nesting** (equal-mean DGP): free vs `equal_mean=TRUE` agree, gap shrinks with n
(|free−equal|_Conger = 3.0e-3 @ n=60 → 5.7e-5 @ n=2000); `t3_hat ≈ 0`. The
free-mean estimator reduces to the one-way pairwise — nests by mean-model
restriction, as claimed.

**B. vs free-mean CS FIML** (complete, 200 reps): mean |moment − FIML|_Conger =
**0.0009** (equal) / **0.0007** (free). Same probability limit (moment = MLE on
complete CS).

**C. t3 correction** (equal-mean DGP, true t3=0): naive (uncorrected) Conger is
spuriously below Fleiss at small n (0.387 vs 0.392 @ n=30), the corrected Conger
tracks Fleiss (0.393); the gap vanishes by n=500. Correction necessary + correct;
magnitude modest for typical R/n.

**D. Calibration** (free-mean DGP, n=200, B=300, truth Conger=0.305):

| mech | method | Conger bias | se_ratio | coverage |
|------|--------|-------------|----------|----------|
| complete | moment | −0.002 | 1.01 | 0.953 |
| complete | fiml   | −0.003 | 1.00 | 0.950 |
| mcar | moment | −0.008 | 1.00 | 0.940 |
| mcar | fiml   | −0.007 | 0.97 | 0.943 |
| mar  | **moment** | **−0.013** | 0.96 | **0.897** |
| mar  | fiml   | −0.001 | 0.99 | 0.940 |

Moment estimator: MCAR-consistent, bootstrap SE calibrated; MAR-biased (coverage
0.90 under benign anchor-MAR — magnitude is mechanism-dependent). FIML stays
consistent under MAR. The textbook moment-vs-MLE split.

## Verdict / fold-in
The estimator is **principled and sound** — it genuinely completes the
FIML-vs-moment grid (the `pairwise + variance_components` cell). It is *not* a
capability gain over FIML: MCAR-only, ~as efficient at best, less efficient at
equal means (it spends R means). Point estimate is trivial; the one non-trivial
piece for a package fit is an **analytic robust IF** to match the package `psi`
convention / `joint_vcov` (reuse the raw-pairwise IF machinery + mean-estimation
and `t3` terms) — the probe used a bootstrap SE instead. Recommend wiring it as
`kappa(estimator = "pairwise", covariance = "variance_components")` only if the
grid is wanted in the public API; otherwise this probe + the quadratic-paper note
stand as the validated record.

## Update (2026-06-15): analytic IF built, validated, and SHIPPED

`analytic_if.R` + `validate_if.R`. The estimator is a stacked M-estimator in
theta = (mu_1..mu_R, psi, e); the expected cross-Jacobian d g_{psi,e}/d mu = 0 at
the truth (E[r]=0), so A is block diagonal and the per-subject influence is closed
form (IF_mu_j = obs*(Y-mu_j)/n_j; IF_psi, IF_e = pooled-pair scores / Ptot; t3 via
d t3/d mu_j = 2(mu_j-mubar), bias correction O(1/n) negligible at first order).
bp_chance is a fixed function of the support (constant, as in the other quadratic
paths). Validation: analytic SE / bootstrap SE = 0.95-1.01 for Conger/Fleiss; for
BP it matches too once `values` is supplied so bp_chance is constant (the
values=NULL range-fallback inflates the bootstrap's BP SE only). MCAR coverage
0.93-0.94 at n=120-200. Shipped as `kappa(estimator = "pairwise_cs")`
(r-package/R/kappa_pairwise_cs.R); package output reproduces this probe to ~1e-18
and satisfies vcov == crossprod(psi)/n^2 exactly.
