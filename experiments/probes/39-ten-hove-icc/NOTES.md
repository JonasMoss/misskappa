# Probe 39 — ten Hove et al. (2025) MLE-RE ICC estimator in C++

**Question.** Can we implement ten Hove, Jorgensen & van der Ark (2025, *MBR*
60:5, "How to Estimate ICCs for Interrater Reliability from Planned Incomplete
Data") recommended estimator — **MLE-RE** — in C++, point estimate only, and
verify it? They are natural reviewer suggestions for the kappa-missing / quadratic
papers (ICC is the continuous analogue of quadratic kappa).

**What MLE-RE is.** `lmer(y ~ 1 + (1|subject) + (1|rater), REML = FALSE)`: FIML
for the two-way *crossed* random-effects model
`y_{sr} = mu + a_s + b_r + e_{sr}`, one obs per observed cell, missing cells
absent (ignorable planned-missing). ICCs are ratios of the three variance
components — `ICC(A,1)=s2s/(s2s+s2r+s2sr)`, `ICC(C,1)=s2s/(s2s+s2sr)`.
The traditional ANOVA mean-squares estimator (the "anova" cue) **cannot** produce
least-squares variance components once cells are missing — that gap is the whole
paper. MLE-RE was their winner on every criterion (bias, RMSE, coverage,
convergence ~100%, seconds to fit).

**Implementation** (`icc_fiml.cpp`, RcppEigen, ~70 lines). The marginal
covariance `V = s2sr·I_n + Z G Z'` is low-rank-structured, so Woodbury +
matrix-determinant-lemma collapse every likelihood eval to one Cholesky of the
`(S+R)×(S+R)` matrix `M = G^{-1} + (1/s2sr) Z'Z`. `mu` profiled by GLS; optim
over `log`-variances driven from R (the only non-C++ piece — trivially
internalizable later). Cheap even at S=200 (M is ~210×210).

**Verification** (`verify.R`, `sweep.R`) — oracle is `lme4`, exactly what the
paper used.
- Single design (S=60, R=12, 50% missing): variance components match lme4 to
  **~1e-6 relative**, logLik to 1e-13, both ICCs to ~1e-6.
- Sweep (S∈{30,200}, R∈{5,12}, 40–83% missing, 5 seeds each, 40 fits): worst
  relative disagreement **7.6e-4**, all from BFGS-vs-lme4 optimizer tolerance,
  not structure. 0 singular fits.

**Delta-method SEs** (`icc_se.R`) — DONE + verified. `V_theta` = inverse observed
information of the *profile* loglik (numeric Hessian via numDeriv of the C++
`icc_loglik`, mu profiled out so the mean's estimation is accounted for); analytic
ICC gradients; `SE = sqrt(grad' V_theta grad)`. Also a logit-transformed CI that
respects `[0,1]`. Validated against an **honest parametric bootstrap** (B=500,
refit each replicate) — the actual sampling SD the SE estimates (NB: not the
Krinsky-Robb "Monte-Carlo CI" the paper sells, which never refits):

| regime | design | ICC | SE_delta / SD_boot |
|---|---|---|---|
| interior | S=60, 50% miss, ICC≈.41/.56 | A,1 / C,1 | **1.03 / 0.99** |
| boundary | S=30, 83% miss, ICC≈.78/.80 | A,1 / C,1 | **0.89 / 0.95** |

i.e. delta-SE is ~exact in the interior and ~11% optimistic for ICC(A,1) near the
ceiling — exactly the boundary degradation predicted. There the Wald CI spills up
to 0.92–0.95 while the logit CI stays bounded and asymmetric. Confirms our read on
ten Hove: a model-based delta SE is enough; the MC-CI buys ~nothing it can't get
from a logit transform.

**Status.** Estimator core + delta SEs DONE + verified. Not yet folded into
package `src/` (would need the vendored-cpp path + an internal optimizer/Hessian
to replace R `optim`/`numDeriv`).

**Open / next.**
- Per Fig. 7, DM vs MC differ negligibly within MLE-RE (the coverage failures are
  MLE-CF/lavaan and ANOVA), so no need to add a Monte-Carlo CI. Optionally match
  `merDeriv`→`car::deltaMethod` numbers exactly (not installed; bootstrap is the
  stronger oracle anyway).
- A real coverage sim (vs nominal 95%) would quantify the boundary optimism;
  out of scope for "add SEs".

---

## Decision (2026-06-14): scope, packaging, SEs

**No C++, no core changes.** MLE-RE is a smooth 3-param profile likelihood; R
`optim` on the probe-39 Woodbury kernel is sufficient. It is sim-only and *not
recommended*, so it is NOT exported and NOT compiled into core — confirmed no
NLOPT in the package, and a 3-param profile needs none. **magmaan / misskappa
core stay untouched.** If it ever needs reuse across studies, expose as an
unexported R internal (cf. `joint_vcov`), never a `kappa()` estimator.

**Home.** Online appendix of the **quadratic/Conger paper**, next to the Gwet
comparison. Claim scoped: MLE-RE is the literature's estimator for
*continuous/quadratic* agreement from incomplete raw data *with inference*
(categorical lit = De Raadt / Warrens / Gwet, point-estimate only).

**Story = structured (MLE-RE) vs saturated (nt_fiml `fit_saturated_normal`),
shown two-sided so it's not a strawman:**
- *model true* (additive Gaussian, homoscedastic, parallel raters): both
  consistent; MLE-RE likely **more efficient** (correctly-specified nested
  3-param submodel of the saturated cov). Own it.
- *model false* (heteroscedastic / congeneric / genuinely categorical): MLE-RE's
  additive ICC ≠ true quadratic kappa → **inconsistent**; nt_fiml robust. Same
  moral as study 29's structured-cat_fiml bias, roles flipped.
- Recommendation rests on "parallel/homoscedastic/Gaussian raters implausible for
  real rating data," not on MLE-RE being bad.

**SEs — robust ARE needed.**
- Point estimate is normality-free (quasi-likelihood: consistent under correct
  2nd moments + MAR). The model-based delta SE built above (inverse profile
  Hessian) is correct ONLY under Gaussianity (validated == parametric bootstrap
  *on the Gaussian DGP*).
- ten Hove's *recommended* interval = merDeriv **robust (sandwich) acov** → delta.
  To quote their method honestly the appendix must use the robust SE, not ours.
- Correct robust SE for a CROSSED design is a **multiway (subject × rater)
  cluster-robust** sandwich (Cameron–Gelbach–Miller) / merDeriv casewise scores —
  NOT one-way (subjects share raters ⇒ not independent).
- Load-bearing in exactly one cell: **non-normal DGP, correct 2nd moments** →
  estimate good, naive CI under-covers, sandwich restores calibration (= the
  robustness lesson). Elsewhere (model true) model-based ≈ robust; (model false)
  bias dominates and SE choice is second-order.
- Hard limit to flag: σ²_r — hence **ICC(A,1)** — is governed by R (few rater
  clusters); no SE gives nominal ICC(A,1) coverage at small R (Fig. 7 + our
  boundary cell). ICC(C,1) is fine.
- Sim ground truth is Monte-Carlo empirical coverage regardless; robust SE only
  matters for MLE-RE's *reported* interval. At build time: call merDeriv (their
  tool) or hand-roll the two-way sandwich on the kernel.
- Decide estimand surface: expose ICC(A,1)/ICC(C,1) (single) and the `,k`
  averaged variants (Table 1) — all are functions of the same three components.
- Positioning (for the papers): MLE-RE is ICC-only and MCAR/ignorable-by-design;
  our contrast is kappa/alpha under genuine MAR. See reading-digest in
  `dev/refs/missing-data/`.
