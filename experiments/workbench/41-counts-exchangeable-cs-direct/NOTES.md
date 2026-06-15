# Probe 41 — direct-on-counts exchangeable CS quadratic kappa

**Question.** Can the exchangeable common-mean compound-symmetry (CS) normal
kappa be fit **directly on counts** — no realizing table — with a robust SE, and
in the (listwise, pairwise, fiml) trichotomy? Does it match the current
table-expansion `nt_fiml` and the existing count kappas?

**Answer.** Yes to all. The CS likelihood collapses onto the per-subject
sufficient statistics \((m_i,s_i,q_i)\); fiml/pairwise/listwise all estimate
\((\mu,\psi,e)\) from those and map to \(\kappa_F=\psi/(\psi+e)\),
\(\kappa_{BP}=1-2e/c_{BP}\). Robust (sandwich) SE for each. Derivation:
`dev/notes/counts-exchangeable-cs-direct.md`.

## Files
- `direct_cs.R` — `cs_fiml`, `cs_pairwise`, `cs_listwise`, `cs_counts()`.
- `verify.R` — comparison tables + calibration sim. `calibration.rds` — results.

## Results

**Q1 EXACT — direct fiml vs table-expansion `nt_fiml`.** Identical to 6 d.p. on
`dat.fleiss1971`; worst \(|{\rm direct}-{\rm table}|\) (point + SE) over 30
random partial-count datasets = **2.6e-7**. The table is provably unnecessary.

**Q2 COMPARE.**
- Complete (`dat.fleiss1971`): fiml = pairwise = listwise = `fleiss_cuzick` =
  `cat_fiml` = table-`nt_fiml`, all F=0.284072 (se 0.109311), BP=0.333889 (se
  0.101876).
- Partial (n=400, 30% MCAR): fiml=0.3675 (=table), pairwise=0.3655,
  `fleiss_cuzick`=0.3680, `cat_fiml`=0.3620, all se≈0.029; listwise=0.3526
  (se 0.061 — discards incomplete rows).

**Q3 CALIB (B=600, n=200, true-CS DGP; truth Fleiss=0.366, BP=0.497).**

| mech | estimator | Fleiss bias | se_ratio | coverage |
|------|-----------|-------------|----------|----------|
| complete | all four | -0.001 | 1.03 | 0.957 |
| mcar | fiml | -0.003 | 1.01 | 0.950 |
| mcar | pairwise | -0.004 | 1.00 | 0.943 |
| mcar | listwise | -0.015 | 0.99 | 0.925 |
| mcar | fleiss_cuzick | -0.003 | 1.01 | 0.952 |
| mar | fiml | -0.001 | 1.02 | 0.962 |
| mar | pairwise | -0.002 | 1.04 | 0.965 |
| mar | **listwise** | **-0.094** | 0.99 | **0.740** |
| mar | fleiss_cuzick | -0.002 | 1.03 | 0.963 |

Robust SE well-calibrated everywhere (se_ratio ≈ 1). Efficiency:
fiml ≈ pairwise ≈ fleiss_cuzick > listwise. Consistency: only **listwise**
breaks under MAR (it conditions on completeness, which correlates with values) —
the classic "listwise needs MCAR". Pairwise survives this *benign* anchor-MAR
but that is mechanism-dependent; fiml is the principled ignorable-missingness
MLE.

## Takeaways
1. Drop the realizing-table detour: fit `nt_fiml` counts directly.
2. Ship the (listwise, pairwise, fiml) trichotomy for counts; all robust-SE'd.
3. `fleiss_cuzick` already behaves like an efficient pairwise moment estimator
   here — the new pairwise (CS ANOVA covariance) matches it closely under MCAR.
4. Raw-data CS (variance_components) models lack a CS pairwise/listwise
   companion — follow-up todo.
