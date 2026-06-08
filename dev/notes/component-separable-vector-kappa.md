# Component-separable vector kappa

Working note for the vector-valued pairwise agreement estimator. The class is
**component-separable vector losses**: a rater's rating is a length-`p` vector,
but disagreement is built from feature-level losses and diagonal feature
weights.

## Data and loss class

For subject `i`, rater `j`, and feature `l`, write the component rating as
`X_ijl` and its observation indicator as `M_ijl`. Rater vectors are
`X_ij = (X_ij1, ..., X_ijp)`.

The implemented loss class is

```text
D(x, y) = h( sum_l v_l d_l(x_l, y_l) / sum_l v_l ),
```

where `v_l >= 0` are diagonal feature weights, at least one `v_l` is positive,
`d_l` is a scalar component loss, and `h` is a scalar transform. The current
factories are:

- `hamming`: `d_l(a,b) = 1[a != b]`, `h(t) = t`.
- `absolute`: `d_l(a,b) = |a - b|`, `h(t) = t`.
- `squared`: `d_l(a,b) = (a - b)^2`, `h(t) = t`.
- `rms`: `d_l(a,b) = (a - b)^2`, `h(t) = sqrt(t)`.

The implementation uses numeric component values. Binary and finite
categorical components are represented numerically; `hamming` is exact
component mismatch.

Non-diagonal feature weights are intentionally excluded. They introduce
cross-feature interactions and break the component-wise missingness reduction
that makes the estimator simple.

## Pairwise component moments

For a within-subject rater pair `j < k`, define the available-component
observed-disagreement numerator and denominator

```text
A_d = n^{-1} sum_i sum_{j<k} sum_l
      M_ijl M_ikl v_l d_l(X_ijl, X_ikl),

B_d = n^{-1} sum_i sum_{j<k} sum_l
      M_ijl M_ikl v_l.
```

The observed disagreement is `d = h(A_d / B_d)`.

For the IPW estimator, replace each within-subject component-pair contribution
by

```text
M_ijl M_ikl / pi_jkl,
pi_jkl = P(M_ijl M_ikl = 1),
```

estimated empirically by the rater-pair-feature observation rate. The estimator
is a Hajek ratio, so the same inverse weights appear in numerator and
denominator.

## Chance moments

The Conger-type chance disagreement uses distinct rater pairs and independent
subjects:

```text
A_C = n^{-2} sum_{i,i'} sum_{j<k} sum_l
      M_ijl M_i'kl v_l d_l(X_ijl, X_i'kl),

B_C = n^{-2} sum_{i,i'} sum_{j<k} sum_l
      M_ijl M_i'kl v_l.
```

The Fleiss-type chance disagreement uses all ordered rater pairs:

```text
A_F = n^{-2} sum_{i,i'} sum_{j,k} sum_l
      M_ijl M_i'kl v_l d_l(X_ijl, X_i'kl),

B_F = n^{-2} sum_{i,i'} sum_{j,k} sum_l
      M_ijl M_i'kl v_l.
```

For IPW, chance contributions use marginal rater-feature observation rates:

```text
M_ijl M_i'kl / (pi_jl pi_kl),
pi_jl = P(M_ijl = 1).
```

Then

```text
d_C = h(A_C / B_C),   d_F = h(A_F / B_F),
kappa_C = 1 - d / d_C,
kappa_F = 1 - d / d_F.
```

The returned order is `(Conger, Fleiss)`.

## Influence functions

The estimator has the same six-moment shape as the scalar continuous
estimator:

```text
(A_d, B_d, A_C, B_C, A_F, B_F).
```

The within-subject observed moments use ordinary empirical IFs. The chance
moments are V-statistics and use the existing row/column kernel IF

```text
phi_i = (mean_i' K(i,i') - psi) + (mean_i' K(i',i) - psi).
```

The delta method first maps each numerator/denominator pair through
`h(A/B)`, then maps `(d, d_C, d_F)` to `(kappa_C, kappa_F)`. For `rms`, the
only extra derivative is `h'(t) = 1 / (2 sqrt(t))`; at `t = 0` the derivative
is treated as zero in the implementation, so RMS inference at perfect
agreement should be considered nonregular.

## Squared-loss connection

With complete data, diagonal weights, and squared component loss, the estimator
matches the mean/covariance contraction used in the CRACKLES pilot. The full
feature-weight version is recorded in
`dev/notes/quadratic-vector-tensor-kappa.md`. Stack
`X = (X_1', ..., X_R')'`, where `X_j` is length `p`, and write block
covariances as `Sigma^{jk}`. For diagonal `V = diag(v)`,

```text
T_V = sum_j tr(V Sigma^{jj}),
B_V = sum_{j,k} tr(V Sigma^{jk}),
G_V = sum_j (mu_j - mu_bar)' V (mu_j - mu_bar).
```

Then the complete-data squared-loss contractions are

```text
kappa_F = (B_V - T_V - G_V) / ((R - 1)(T_V + G_V)),
kappa_C = (B_V - T_V)     / ((R - 1)T_V + R G_V).
```

Cross-feature covariances do not enter the diagonal-weight estimand directly.
They matter only for model-based missing-data estimators that reconstruct a
full covariance matrix.

## Deferred categorical FIML

Categorical FIML for vector ratings is deferred. For finite components, the
natural model is a saturated multinomial distribution over full rater-feature
profiles. Missing components are latent entries of that full profile. This is
possible in principle, and is the direct categorical counterpart of the current
raw-rating FIML, but the state space is combinatorial:

```text
prod_l C_l^R
```

for `R` raters and feature category counts `C_l`. A practical cat-FIML pass
will need explicit support pruning, pattern indexing, and documentation of the
full-profile estimand. It is not part of the first component-separable IPW
implementation.
