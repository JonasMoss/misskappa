# Quadratic vector kappa: tensor formulation

This note records the covariance-and-mean form of the vector-valued quadratic
kappas. It is the full-weight counterpart of the component-separable vector
note: once the loss is squared vector distance, the coefficient is a simple
smooth function of the stacked mean and covariance.

## Setup

For subject `i`, rater `j = 1, ..., R`, and feature vector dimension `p`, let

```text
X_ij in R^p,
X_i = (X_i1', ..., X_iR')' in R^(R p).
```

Let

```text
mu = E X_i,
Sigma = Cov(X_i),
W = W' >= 0,   W in R^(p x p).
```

The vector loss is

```text
d_W(x, y) = (x - y)' W (x - y).
```

The rater blocks of `mu` are `mu_j`; the block covariance between raters `j`
and `k` is `Sigma_jk`.

## Tensor contractions

Let `1_R` be the length-`R` all-ones vector and

```text
P_R = I_R - (1_R 1_R') / R.
```

The three quadratic summaries are

```text
B_W = tr{ (1_R 1_R' tensor W) Sigma },
T_W = tr{ (I_R tensor W) Sigma },
G_W = mu' (P_R tensor W) mu.
```

Equivalently, in block notation,

```text
B_W = sum_{j,k} tr(W Sigma_jk),
T_W = sum_j tr(W Sigma_jj),
G_W = sum_j (mu_j - mu_bar)' W (mu_j - mu_bar),
mu_bar = R^{-1} sum_j mu_j.
```

These are the only population objects used by the coefficients. Cross-feature
covariances enter through `W` and the block covariance matrix; no pairwise-loss
enumeration is needed once `(mu, Sigma)` is available.

## Coefficients

The vector quadratic Conger and Fleiss coefficients are

```text
kappa_C^W = (B_W - T_W) /
            { (R - 1) T_W + R G_W },

kappa_F^W = (B_W - T_W - G_W) /
            { (R - 1) (T_W + G_W) }.
```

For `p = 1` and `W = 1`, these reduce to the scalar quadratic formulas in the
quadratic-kappa paper. For diagonal `W`, they reduce to the weighted
same-feature contraction used by the first CRACKLES vector pilot.

## Estimators

The internal R backend `kappa_vector_quadratic()` implements two moment sources:

- `method = "pairwise"`: pairwise-available covariance moments. This is the
  distribution-free MCAR-oriented analogue of the quadratic paper's pairwise
  estimator.
- `method = "nt_fiml"`: saturated normal EM/FIML for `(mu, Sigma)`, with the
  same sandwich score/information machinery used by `kappa_quadratic_fiml()`.

Both estimators apply the same smooth map

```text
(mu, Sigma) -> (B_W, T_W, G_W) -> (kappa_C^W, kappa_F^W).
```

The delta method is correspondingly simple. With

```text
A_B = 1_R 1_R' tensor W,
A_T = I_R tensor W,
A_G = P_R tensor W,
```

the covariance-gradient components are the vech gradients of `A_B` and `A_T`,
and the mean-gradient component is

```text
dG_W / dmu = 2 A_G mu.
```

The coefficient Jacobian with respect to `(B_W, T_W, G_W)` is the same
three-summary Jacobian as the scalar quadratic estimator.

## Missingness boundary

This route is not the general component-separable IPW estimator. It is a
quadratic covariance plug-in. Full feature-weight matrices are natural here
because the full covariance matrix contains cross-feature terms.

For non-quadratic losses such as Hamming, absolute loss, max-component loss, or
whole-profile exact match, a full `W` does not reduce to `(mu, Sigma)` and needs
a separate missing-data estimand or a full categorical/profile model.
