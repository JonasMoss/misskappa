# Alpha Cat-FIML triage

Date: 2026-06-09

This note records the current implementation audit for the categorical
maximum-likelihood coefficient-alpha estimator used in the alpha-missing paper.
The paper labels this route Cat-FIML / Cat-ML: a saturated multinomial likelihood
over discrete full item-response profiles, fitted by EM from incomplete rows,
with coefficient alpha read from the fitted scored distribution.

## Implemented algorithm

The R surface is `alpha(x, estimator = "cat_fiml", values = ..., em_options =
...)`. It routes through `alpha_cat_fiml()`, then `rcpp_alpha_raw(method =
"fiml")`, then `misskappa::estimate_alpha_fiml()`.

For a `n x p` item matrix with categories `0, ..., C - 1`, the full data model is
a saturated multinomial distribution
\[
  \pi_y,\qquad y \in \{0,\ldots,C-1\}^p,\qquad \sum_y \pi_y = 1.
\]
For subject `i`, let `A_i` be the set of full profiles compatible with the
observed entries in row `i`. The observed likelihood contribution is
\[
  L_i(\pi) = \sum_{y \in A_i} \pi_y.
\]
The EM update is therefore
\[
  q_i(y) = \frac{\pi_y}{\sum_{z \in A_i}\pi_z}\mathbf 1\{y\in A_i\},\qquad
  \pi_y^{new} = \frac{1}{n}\sum_i q_i(y).
\]
The C++ implementation groups identical observed patterns, enumerates
compatible completions for each group, applies this update to the grouped
counts, prunes only after convergence, and renormalizes the retained support.

Given score vector `v`, the fitted distribution is mapped to alpha through the
scored full-profile moments:
\[
  T_1 = \operatorname{Var}_\pi\left(\sum_{j=1}^p v_{Y_j}\right),\qquad
  T_2 = \sum_{j=1}^p \operatorname{Var}_\pi(v_{Y_j}),
\]
\[
  \alpha(\pi) = \frac{p}{p-1}\left(1 - \frac{T_2}{T_1}\right).
\]
This is exactly Cronbach's alpha for the covariance matrix implied by the
fitted multinomial distribution.

The standard error uses the observed information for the saturated multinomial
on the simplex. The implementation chooses the largest retained probability as
a reference cell, builds reduced observed scores for each observed-pattern
group,
\[
  s_g = \frac{\partial}{\partial\pi_*}\log\left(\sum_{y\in A_g}\pi_y\right),
\]
forms the reduced observed information as `sum_g n_g s_g s_g'`, pseudo-inverts
it with `info_rcond`, maps it back through the simplex Jacobian, and applies the
delta method with the alpha gradient. Per-subject influence rows are reconstructed
from the same grouped scores; the tests check that `crossprod(psi) / n^2`
reconstructs `vcov`.

## Audit results

No point-estimate bug was found in the Cat-FIML implementation.

Direct checks completed in this pass:

- Complete data: the EM degenerates to one compatible completion per row, and
  `estimate_alpha_fiml()` matches the pairwise/complete covariance alpha in the
  existing unit test.
- Independent EM: a standalone R implementation of the saturated observed-data
  EM, written without using the C++ code path, matched `alpha(estimator =
  "cat_fiml")` on three incomplete fixtures:
  - `tiny_missing`: package `0.896907216433`, independent `0.896907216440`,
    absolute difference `7.1e-12`.
  - `four_by_three`: package `-1.009293121056`, independent
    `-1.009293119408`, absolute difference `1.6e-9`.
  - `four_items`: package `-0.392970517177`, independent `-0.392970517146`,
    absolute difference `3.1e-11`.
- Regression test: `tests/unit/alpha_test.cpp` now pins the `tiny_missing`
  fixture to the independent EM value.
- Current API smoke runs completed into `/tmp` for the categorical smoke,
  categorical calibration, paper-facing alpha simulation, K&P categorical probe,
  and normal-FIML stress runner.

An implementation issue was found outside the estimator itself: several
alpha-missing simulation runners still called the old R API
`alpha(method = "available" | "fiml", type = ...)` and
`alpha_continuous(se_type = ...)`. The current package API is
`alpha(estimator = "pairwise" | "cat_fiml" | "nt_fiml")`. The paper-facing
and Cat-FIML-relevant runners were updated so the evidence can be regenerated
against the current package.

## Interpretation of poor Cat-FIML behavior

The poor Cat-FIML rows in the alpha-missing simulations are most consistent with
finite-sample behavior of the saturated multinomial MLE, not with a missed
implementation detail.

The paper preview uses six three-category items, so the saturated state space is
`3^6 = 729` cells at `n = 100`. Incomplete rows contribute only marginal
constraints over large compatible sets, and alpha is a nonlinear ratio of fitted
second moments. The saturated MLE is therefore high-dimensional relative to the
sample size even before moving to the larger `5^6 = 15625` diagnostic. Complete
data removes the missing-profile uncertainty and Cat-FIML collapses to the usual
empirical covariance alpha; the bad rows appear when sparse incomplete profiles
must be distributed over a large cell table.

Normal FIML can look much better in the same ordinal cells because it imposes a
low-dimensional Gaussian moment model on the scored data. That is useful as a
working model, but it is not evidence that the saturated discrete likelihood is
miscomputed. It is shrinkage/model structure buying finite-sample stability.

## Remaining guardrails

- Treat Cat-FIML as theoretically clean under discrete MAR but practically
  small-support. The manuscript should say this directly.
- Report Cat-FIML only for small `C^p` cells, and keep timing/support diagnostics
  in the reproducibility record.
- Do not present the `n = 100`, `3^6` preview rows as final estimator ranking.
  Add a larger `n` grid before using the magnitude of Cat-FIML bias as evidence.
- If Cat-FIML remains a paper focus, the natural follow-ups are explicit
  second-order bias correction, bootstrap/jackknife correction, mild smoothing or
  penalized multinomial likelihood, or replacing the saturated table with a
  parsimonious loglinear/latent-class/ordinal-copula model.
