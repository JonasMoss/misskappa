# Probe 40 — variance-component ("ten Hove") alpha: oracle parity + CI coverage

**Question.** Commit `ff96f5c` ("Add variance-component normal FIML") added the
scalar compound-symmetry (CS) covariance option `Sigma = psi J + e I` to the
package's normal-FIML path, for both quadratic `kappa()` and continuous
`alpha()`. Two things were unverified after that commit:
  1. **Parity.** The package fits CS by a *saturated EM with a CS projection in
     the M-step* (`project_variance_components`). Probe 39 only ever verified a
     *different* implementation (Woodbury + R `optim`). Does the package's
     EM-projection estimate match the lme4 oracle on **missing** data?
  2. **Coverage.** No CI calibration of the package's analytic sandwich existed
     for alpha at all.

Runner: `parity_and_coverage.R` (+ a large-n MAR inconsistency check inline).

## Algebraic fact used throughout

Alpha depends on Sigma only through `tr(Sigma)` and `sum(Sigma)`, and the CS
projection **preserves both** (projected `total = T1`, projected `trace = T2`;
verified in the source comment). So on the *moments* the CS route and the
saturated route target the **same alpha** even when CS is false. The open
empirical question was whether the constrained **EM fixed point** keeps that
invariance once missingness couples the E-step to the wrong covariance.

## PART 1 — parity (DONE, PASS)

Package `covariance="variance_components"` alpha vs `lmer(value ~ 0 + item +
(1|subject), REML=FALSE)` -> alpha from `(psi, e)`. p=6, n=600, 30% missing.

| DGP | mech | worst \|pkg - lme4\| |
|---|---|---|
| cs / congeneric | mcar / mar | **5.4e-8** |

⇒ The package EM-projection estimate **is** the lme4 CS-constrained MLE on
missing data (MCAR and MAR, well- and mis-specified). Point estimate verified
against the established oracle — closes the gap probe 39 left.

## PART 2 — CI coverage (DONE)

B=800, n=250, p=6, 95% nominal. `truth` = population alpha from the true Sigma.
Methods: `pairwise` (available-case), `saturated` (nt_fiml), `varcomp`
(nt_fiml + CS). DGPs: `cs` (model true) and `congeneric` (CS false). Coverage
shown for the natural-scale Wald CI; Fisher is similar except where bias
dominates (then it is slightly worse, as expected).

| DGP | mech | method | bias | sd_emp | mean_se | cov_wald |
|---|---|---|---|---|---|---|
| cs | complete | all three | -.001 | .014 | .014 | .94 |
| congeneric | complete | all three | -.002 | .013 | .013 | .95 |
| cs | mcar | pairwise | -.003 | .019 | .018 | .94 |
| cs | mcar | saturated | -.001 | .016 | .016 | .94 |
| cs | mcar | **varcomp** | **-.001** | **.016** | **.016** | **.95** |
| congeneric | mcar | pairwise | -.002 | .018 | .017 | .94 |
| congeneric | mcar | saturated | .000 | .015 | .014 | .95 |
| congeneric | mcar | **varcomp** | **.000** | .015 | .015 | **.95** |
| cs | mar | pairwise | **-.028** | .023 | .022 | **.785** |
| cs | mar | saturated | -.001 | .018 | .017 | .95 |
| cs | mar | **varcomp** | **-.001** | .017 | .017 | **.95** |
| congeneric | mar | pairwise | **-.023** | .022 | .021 | **.826** |
| congeneric | mar | saturated | -.001 | .017 | .016 | .93 |
| congeneric | mar | **varcomp** | **-.018** | .017 | .017 | **.87** |

### Large-n MAR confirmation (n=60000, 3 seeds)

| DGP | varcomp bias | saturated bias |
|---|---|---|
| cs (correct) | +.0003 | +.0003 |
| congeneric (CS false) | **-.0177** | -.0003 |

The MAR+congeneric bias **persists at large n** ⇒ genuine asymptotic
inconsistency, not finite-sample.

## Verdict — same 2x2 (specification x missingness) as quadratic kappa (probe 39)

1. **MCAR ⇒ varcomp consistent AND best-calibrated, even under misspecification.**
   Bias ~0 and coverage ~.95 for both CS and congeneric. The CS-projection
   moment invariance survives MCAR. Efficiency: marginally tighter SD than
   saturated when CS is true, equal otherwise. The sandwich SE is well
   calibrated (se_ratio ~.99).
2. **MAR + correctly specified (CS true) ⇒ varcomp consistent, coverage ~.95.**
   MAR alone is not the problem.
3. **MAR + CS misspecified (congeneric) ⇒ varcomp INCONSISTENT** (bias -.018,
   coverage ~.87), while **saturated nt_fiml stays consistent** (~.93-.95). The
   saturated path is the fix here. Mirrors probe-39 quadratic kappa exactly.
4. **pairwise available-case is MAR-invalid** throughout (bias -.02 to -.03,
   coverage .73-.83) — the standard point. Fine under MCAR.

⇒ For the alpha-missing paper: the variance-component (ten Hove / parallel-items)
alpha CI is the efficient, well-calibrated choice under MCAR (consistent even
when items aren't parallel, because alpha = f(tr, sum) and CS preserves both),
and under MAR *when items are genuinely parallel*. Under MAR with non-parallel
(congeneric) items the saturated FIML alpha is the robust choice. Same moral as
the quadratic-kappa structured-vs-saturated story, transported to alpha.

**Open.** Non-Gaussian MAR (mixture/copula) should degrade saturated too (cf.
probe 39 point 4 and study 23); not run here. A boundary/small-R cell (psi near
0) is the place CS coverage would degrade most; not stressed here.
