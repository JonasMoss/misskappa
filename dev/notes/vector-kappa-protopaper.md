# Vector-valued agreement with separable feature weights — proto-paper

Status: **proto-paper, in preparation.** This is a short, paper-shaped overview of
*our* construction. It is not yet a manuscript and the API it documents is frontier.
The dense technical notes are the appendix: see
[`component-separable-vector-kappa.md`](component-separable-vector-kappa.md),
[`quadratic-vector-tensor-kappa.md`](quadratic-vector-tensor-kappa.md), and
[`crackles-vector-kappa.md`](crackles-vector-kappa.md).

## What this is, and what is new

A rating is a length-`p` vector `X_ij = (X_ij1, ..., X_ijp)` (subject `i`, rater `j`,
feature `l`). We want a chance-corrected multirater agreement coefficient for such
vector ratings, with **per-feature weights**, that also degrades gracefully when some
components are missing.

The relation to existing work:

- The general Fréchet-variance framing of multirater agreement (Moss 2024,
  *Psychometrika*) and the multirater / vector-valued-data treatment in the kappa-sd
  paper ("Chance-corrected measures of agreement with multiple raters and vector-valued
  data") give the *Fréchet* scaffolding and the complete-data vector idea. They do
  **not** cover the construction below.
- Janson & Olsson (2001) give an interval/nominal multivariate agreement coefficient;
  it is the closest precedent, but it is a single fixed multivariate loss, not a
  weighted component-separable family with a missing-data theory.

What is new here, and therefore what this proto-paper would claim:

1. The **component-separable weighted loss class** (below), with diagonal feature
   weights and a per-feature scalar loss, as the vector generalization of weighted
   kappa.
2. **Missing-component estimators** for it: pairwise-available (MCAR) and IPW (MCAR
   with feature- and rater-pair-varying observation rates), with influence-function
   standard errors.
3. The **full-feature-weight quadratic special case**, where a symmetric PSD weight
   matrix `W` couples features and the coefficient becomes a smooth function of the
   stacked mean and covariance, admitting a normal-theory FIML estimator.

## Loss class

The component-separable weighted vector loss is

```text
D(x, y) = h( sum_l v_l d_l(x_l, y_l) / sum_l v_l ),
```

with diagonal feature weights `v_l >= 0` (at least one positive), a scalar component
loss `d_l`, and a scalar transform `h`. The implemented factories:

- `hamming`:  `d_l(a,b) = 1[a != b]`,  `h(t) = t`.
- `absolute`: `d_l(a,b) = |a - b|`,    `h(t) = t`.
- `squared`:  `d_l(a,b) = (a - b)^2`,  `h(t) = t`.
- `rms`:      `d_l(a,b) = (a - b)^2`,  `h(t) = sqrt(t)`.

Diagonal weights are deliberate: they keep disagreement a sum of per-feature terms, so
a missing component only drops that component from the relevant pair, and the estimator
stays a ratio of feature-level moments. Non-diagonal weights break this reduction and
are handled only in the quadratic special case (next section).

The Conger- and Fleiss-type coefficients are the usual `kappa = 1 - d / d_chance`, with
observed disagreement `d`, distinct-rater-pair chance `d_C`, and pooled-margin chance
`d_F`. Estimands, pairwise/IPW moment forms, and the six-moment influence functions are
in [`component-separable-vector-kappa.md`](component-separable-vector-kappa.md).

## Quadratic special case (full feature weights)

When the loss is squared vector distance with a symmetric PSD weight matrix `W`,

```text
d_W(x, y) = (x - y)' W (x - y),
```

the coefficient depends on the data only through the stacked mean `mu` and covariance
`Sigma`. With `B_W = sum_{j,k} tr(W Sigma_jk)`, `T_W = sum_j tr(W Sigma_jj)`, and
`G_W = sum_j (mu_j - mu_bar)' W (mu_j - mu_bar)`,

```text
kappa_C^W = (B_W - T_W) / { (R-1) T_W + R G_W },
kappa_F^W = (B_W - T_W - G_W) / { (R-1)(T_W + G_W) }.
```

For `p = 1, W = 1` these reduce to the scalar quadratic kappas; for diagonal `W` they
reduce to the same-feature contraction of the separable squared loss. This is the only
place cross-feature couplings enter the estimand. Full derivation and the delta-method
gradient are in [`quadratic-vector-tensor-kappa.md`](quadratic-vector-tensor-kappa.md).

## Estimators and exposure

| Path | Weights | Missingness | Status |
| --- | --- | --- | --- |
| `kappa(array, estimator = "pairwise")` | diagonal `v_l` | complete / MCAR | public |
| `kappa(array, estimator = "ipw")` | diagonal `v_l` | MCAR, component-varying rates | public |
| `kappa_vector_quadratic(method = "pairwise")` | full `W` | MCAR | **internal** |
| `kappa_vector_quadratic(method = "nt_fiml")` | full `W` | ignorable, normal-theory | **internal** |

The diagonal-weight pairwise/IPW estimators are reachable through the public `kappa()`
3-D-array dispatch. The full-`W` quadratic estimators (`kappa_vector_quadratic()`) are
internal and not yet exposed.

## Open problems

- **Categorical full-profile FIML** for vector ratings is deferred: the natural
  saturated model has state space `prod_l C_l^R`, which needs support pruning and a
  documented full-profile estimand.
- **Nonregular RMS at perfect agreement:** `h'(t) = 1/(2 sqrt(t))` blows up at `t = 0`;
  the implementation sets the derivative to zero there, so RMS inference at perfect
  agreement is nonregular and should be flagged.
- **FIML on degenerate data:** the normal-theory FIML path inverts per-pattern
  covariance blocks and is singular when a feature has no within-pattern variance
  (e.g. a site all 28 observers score identically). The complete binary CRACKLES data
  triggers this; a regularized or rank-aware covariance step is needed before this path
  is robust enough to expose.
- **Estimand under MAR** for the diagonal IPW estimator: deletion should depend on
  observed labels/ratings, not latent truth; the CRACKLES missingness illustration is
  synthetic (the shipped data are complete).
