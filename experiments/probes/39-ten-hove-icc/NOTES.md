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

**Status.** Estimator core DONE + verified. Not yet folded into package `src/`
(would need the vendored-cpp path + an internal optimizer to replace R `optim`).

**Open / next.**
- Variance/SE: paper uses `merDeriv` robust acov → delta-method *or* Monte-Carlo
  CI. NB their own Fig. 7 shows DM vs MC differ negligibly *within* MLE-RE; the
  coverage failures are MLE-CF (lavaan, fixed-rater) and ANOVA. So a delta-method
  SE on the variance components is enough; skip the MC-CI theater.
- Decide estimand surface: expose ICC(A,1)/ICC(C,1) (single) and the `,k`
  averaged variants (Table 1) — all are functions of the same three components.
- Positioning (for the papers): MLE-RE is ICC-only and MCAR/ignorable-by-design;
  our contrast is kappa/alpha under genuine MAR. See reading-digest in
  `dev/refs/missing-data/`.
