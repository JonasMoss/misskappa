# Counts Weighting: F&C, GLS, Pooled Pairs, and irrCAC

Context: `kappa_counts(estimator = "pairwise")` currently uses a pooled-pairs
moment estimator. That is not the Fleiss--Cuzick unequal-judges convention and
should not be described as the principled count-data analogue.

## Estimators

Let `N_i` be an observed category-count row, `r_i = sum_c N_ic`, and `L` the
disagreement matrix. For rows with `r_i >= 2`, define

```text
d_i = (N_i' L N_i - diag(L)' N_i) / (r_i (r_i - 1)).
```

The competing count moment estimators mainly differ by row weight:

| Name | Row weight for `d_i` | Motivation |
|---|---:|---|
| Subject-wise / irrCAC | `1` | Average each subject's pair ratio equally. |
| Fleiss--Cuzick (F&C) | `r_i - 1` | Unequal interchangeable judges; binary kappa = one-way random-effects ICC. |
| GLS / IVW heuristic | approximately `r_i` | Row-level precision grows roughly linearly in the number of ratings. |
| Pooled pairs (current C++) | `r_i (r_i - 1)` | Treat every observed within-subject pair as a datum. |

When all row sums are constant these collapse to the same point estimate. When
row sums vary, pooled pairs over-weights long rows because pairs within a
subject share ratings and are correlated. Unit-weighted subject-wise averaging
goes the other way and ignores the extra precision in longer rows.

## F&C Convention

The local reference PDF is
`dev/refs/missing-data/Fleiss and Cuzick 1979 - The Reliability of Dichotomous Judgments - Unequal Numbers of Judges per Subject.pdf`.
Fleiss and Cuzick (1979, Eq. 3) define, for dichotomous counts with
`p_i = x_i / n_i`,

```text
kappa = 1 - sum_i n_i p_i (1 - p_i) /
            { N (nbar - 1) pbar (1 - pbar) }.
```

Since `N (nbar - 1) = sum_i (n_i - 1)`, this is the weighted average of row
pair-disagreement proportions with weights `n_i - 1`, divided by the usual
Fleiss chance disagreement from the pooled rating-token margin. In our
weighted multi-category setting, the direct extension is

```text
d_FC     = sum_i (r_i - 1) d_i / sum_i (r_i - 1),
p_FC     = sum_i N_i / sum_i r_i,
dF_FC    = p_FC' L p_FC,
kappa_F  = 1 - d_FC / dF_FC,
kappa_BP = 1 - d_FC / d_BP.
```

This should be the preferred public count moment estimator unless a stronger
reason emerges. It is citable, agrees with ordinary Fleiss for balanced counts,
and has the cleanest bridge to the unequal-judges literature.

## GLS Motivation

A GLS story starts from the row-level unbiased estimator `d_i` and asks how to
combine rows with different `r_i`. The exact inverse-variance weight is
model-dependent. A generic exchangeable-row variance decomposition has the form

```text
Var(d_i | r_i) ~= a + b / r_i + c / {r_i (r_i - 1)},
```

where `a` is between-subject heterogeneity in full-row disagreement and the
other terms are within-row sampling noise. If the first-order within-row term
dominates, inverse-variance weights are approximately proportional to `r_i`;
if between-subject heterogeneity dominates, the optimum moves toward equal row
weights. The key point is that information grows roughly linearly, not
quadratically, in `r_i`.

Thus `r_i` is a useful GLS/IVW heuristic, while F&C's `r_i - 1` is the
historical finite-sample convention. Quick simulations with `r_i in {2,20}`
showed `r_i` and `r_i - 1` essentially tied, pooled pairs slightly worse, and
unit-weighted subject-wise much worse.

## irrCAC Convention

The installed `irrCAC` 1.4 implementation uses the unit-weighted subject-wise
convention for distribution/count input. In `irrCAC::fleiss.kappa.dist()`:

```r
pa <- sum(sum.q[ri.vec >= 2] /
          ((ri.vec * (ri.vec - 1))[ri.vec >= 2])) / n2more
pi.vec <- t(t(rep(1/n, n)) %*%
            (agree.mat / (ri.vec %*% t(rep(1, q)))))
```

So the observed agreement is an equal average of row ratios, and the chance
margin is the equal average of row-normalized category proportions. This is
software parity to mention, not the estimator to adopt. It agrees with F&C and
pooled pairs only when all `r_i` are equal.

## Sufficient Conditions

Use these as the clean sufficient conditions for count-data moment estimators.
The first three conditions are the sufficiency claim: labeled rater data can be
collapsed to counts without losing information for the Fleiss/BP count target.
The remaining conditions are the missingness conditions needed for the observed
counts to estimate that same target.

1. Subjects are iid.
2. The target is a symmetric count target, not a rater-specific target. For a
   fixed full panel size `R`, the full-data estimand depends on subject `i`
   only through the full count vector `Z_i`, via

   ```text
   D(Z_i) = (Z_i' L Z_i - diag(L)' Z_i) / {R (R - 1)}
   P      = E[Z_i / R].
   ```

   Fleiss-type and Brennan--Prediger-type coefficients of the form
   `1 - E[D(Z_i)] / D_chance(P)` satisfy this. Rater-specific Cohen-type
   contrasts do not.
3. Rater labels carry no additional model information for this target. A
   sufficient concrete condition is exchangeability:
   `(X_i1*, ..., X_iR*)` is invariant under permutations of rater labels. In
   likelihood language, the conditional distribution of the labeled rating
   vector given `Z_i` is parameter-free. Then `Z_i` is sufficient for the
   full-data count target.
4. The observed row size satisfies `r_i >= 2` with positive probability. Given
   `Z_i` and `r_i`, the observed count vector `N_i` is generated by a
   simple random subset of `r_i` raters from the `R` full ratings:

   ```text
   P(N_i = n | Z_i = z, r_i = r)
     = prod_c choose(z_c, n_c) / choose(R, r).
   ```

   Equivalently, conditional on `r_i`, the missing subset is exchangeable over
   rater labels and does not depend on rating identities beyond `Z_i`.
5. For the simple moment estimators, `r_i` is independent of `Z_i` or fixed by
   design independently of ratings. More generally, for the chosen row weight
   `w(r_i)`, it is enough that the weighting moments factor:

   ```text
   E[w(r_i) d(Z_i)] / E[w(r_i)] = E[d(Z_i)]
   E[r_i Z_i / R] / E[r_i] = E[Z_i / R],
   ```

   but `r_i independent of Z_i` is the readable sufficient condition.

Under these conditions F&C, GLS/`r_i`, and pooled pairs are consistent for the
same exchangeable count target; they differ in efficiency. If `r_i` depends on
the unobserved full composition `Z_i`, the moment estimators generally target a
selection-weighted estimand. That is where the counts FIML model belongs.

## Implementation Decision

- Add C++ support for both F&C (`r_i - 1`) and GLS/`r_i` count weighting.
- Make the R-package public default F&C.
- Keep pooled pairs only as a legacy/internal comparison if needed; do not
  present it as the principled count estimator.
- Validate balanced-count parity against `irrCAC` / Fleiss 1971, and add an
  unequal-row test showing F&C differs from both current pooled pairs and
  `irrCAC` subject-wise.
