# Counts-FIML: model, assumptions, and what the estimator actually does

Working note for the appendix. We need to be crystal clear because the
counts-format estimator superficially "looks like" an MAR-handling FIML
but its assumptions are much narrower than the rater-identified FIML.

## Setup

Each of `n` subjects is rated by `R` raters. The full data for subject
`i` is `(X_{i,1}, ..., X_{i,R}) ∈ {1, ..., C}^R`. The counts
representation aggregates rater identity away:

```
Z_i(c) = |{j : X_{i,j} = c}|,    sum_c Z_i(c) = R.
```

Observed data is partial counts `N_i ≤ Z_i` (componentwise) with
`sum_c N_i(c) = r_i ≤ R`.

## Full-data model: exchangeable raters

**Assumption (E).** Raters are exchangeable: the joint distribution of
`(X_{i,1}, ..., X_{i,R})` is invariant under permutation of the rater
index. The composition `Z_i` is therefore a sufficient statistic for the
joint, and we model

```
Z_i ~ θ
```

where `θ` is a probability distribution over the simplex of compositions
of `R` into `C` non-negative integers (size `(R+C-1 choose C-1)`).

Note that `θ` need **not** be of the form `Multinomial(R, p)` for a
single `p` vector. Allowing general `θ` covers, e.g., mixtures of
multinomials (different subject populations with different rater
agreement levels). The parameter space is strictly larger than `{p ∈
Δ^C}`.

What **is** ruled out by (E): rater-specific category distributions.
Specifically, if rater 1 tends to give high categories and rater 2 tends
to give low ones, the full joint is not exchangeable and aggregating to
`Z_i` discards exactly the information needed to identify the per-rater
tendencies. For that case use rater-identified data and the FIML
estimator `estimate_fiml`.

## Observation model: uniformly-random subsampling

**Assumption (S).** Given full counts `Z_i` and observed total `r_i`,
the observed counts `N_i` are a uniformly-random sample without
replacement of `r_i` raters from the `R` available. The likelihood
factor is the multivariate hypergeometric

```
P(N_i | Z_i, r_i) = [prod_c (Z_i(c) choose N_i(c))] / (R choose r_i).
```

This is the natural sampling model when missingness is **MCAR** at the
rater level and raters are exchangeable: every subset of `r_i` raters is
equally likely to be the observed one. It is also the right likelihood
factor under the weaker MAR-on-counts condition

```
P(r_i, which-raters-missing | Z_i) = P(r_i, which-raters-missing | N_i)
```

because the second piece (which raters are missing) integrates out
uniformly under exchangeability, leaving only `P(r_i | Z_i)`. By Rubin's
ignorability theorem, the EM treats `P(r_i)` as a nuisance and we don't
need to model it.

What **is** ruled out by (S): non-uniform sampling. Specifically:

- Rater-specific dropout probabilities (e.g., rater 6 is always missing,
  raters 1-5 are always present).
- Value-dependent dropout (e.g., low-category raters drop out more
  often). This is MNAR.

If either holds, the multivariate hypergeometric likelihood is
misspecified, and the FIML estimate is generally biased.

## Likelihood and EM

Combining (E) and (S), the observed-data likelihood for one subject is

```
P(N_i | θ, r_i) = (1 / (R choose r_i)) · Σ_e θ[N_i + e] · prod_c (N_i(c) + e(c) choose N_i(c))
```

summed over compositions `e` of `R - r_i` into `C` parts. The
`(R choose r_i)` factor is constant in `θ` and dropped from the
optimisation.

**E-step.** Posterior over completions:

```
P(e | N_i, θ_t) = θ_t[N_i + e] · m(N_i, e) / Σ_{e'} θ_t[N_i + e'] · m(N_i, e')
```

where `m(N_i, e) = prod_c (N_i(c) + e(c) choose N_i(c))`.

**M-step.** Expected counts of each composition `z` accumulated over all
subjects, then renormalised:

```
θ_{t+1}[z] = (Σ_i E[1{Z_i = z} | N_i, θ_t]) / n.
```

## Identifiability

The parameter `θ` lives on a simplex of dimension `(R+C-1 choose C-1) - 1`.
For complete data (`r_i = R` for all `i`) this is fully identified — `θ`
is just the empirical distribution. For partial data the observed
marginal `P(N_i | θ, r_i)` maps `θ` into a smaller simplex of dimension
`(r_i + C - 1 choose C - 1) - 1`, so different `θ`'s can produce the
same observed marginal. `θ` is **not generally identified from partial
counts**.

However, the kappa functional is a function of *linear* moments of `θ`:

```
pa(θ) = Σ_z θ[z] · d_z,    d_z = (z^T L z - diag(L)^T z) / (R(R-1))
p_hat(θ)_c = Σ_z θ[z] · z(c) / R
pe_F(θ) = p_hat(θ)^T L p_hat(θ)
```

The expectations `pa(θ)` and `p_hat(θ)` are linear functionals that the
observed-data marginals **do** identify (they're functions of moments
that are estimable from observed counts under (S)). So **kappa is
identified even though `θ` is not.**

What this means in practice: different EM initialisations may give
different `θ_hat` for the same data, but they all give the same kappa.
The asymptotic variance via Louis observed information depends on `θ`
through the score functions; in the non-identified regime the
information matrix is rank-deficient and we use the pseudo-inverse on
the reference-removed parametrisation.

## MCAR vs MAR for counts: where the line actually sits

Under (E) and (S), the only remaining degree of freedom in the
missingness mechanism is `P(r_i | Z_i)`. Three regimes:

| Regime | Condition | Behaviour |
|---|---|---|
| **MCAR-counts** | `r_i ⊥ Z_i` | `estimate_available_counts` and `estimate_fiml_counts` consistent; asymptotically equivalent. |
| **MAR-counts** | `P(r_i \| Z_i) = P(r_i \| N_i)` | FIML still consistent (ignorability). Available-case also consistent: the marginal `p_c = Σ_i N_i(c) / Σ_i r_i` is unbiased because `E[N_i(c) \| r_i] = r_i · p_c` under exchangeable iid raters. |
| **MNAR-counts** | `P(r_i \| Z_i)` depends on unobserved `Z_i - N_i` | Both estimators biased; needs a joint model. |

The MAR / MCAR distinction is **technically real but practically
narrow** for counts: it requires that dropout count `r_i` depend on the
specific composition of observed counts (e.g., subjects whose first few
raters disagreed are more likely to lose later raters), which is unusual
in practice. The bigger driver of `r_i` variation is typically external
(time pressure, missing records), independent of values.

## What FIML buys over available-case

Under MCAR / MAR + exchangeable raters, both estimators are consistent
and have the same asymptotic limit. Differences in finite samples:

- When all `r_i = R` (complete counts), they coincide **exactly**.
- When `r_i` varies, FIML downweights short rows (whose `N_i` are
  consistent with many compositions) and upweights long rows (whose
  `N_i` pins down `Z_i` more tightly). This **can** reduce MSE.
- The variance estimators are different but asymptotically equivalent;
  FIML's Louis information uses the EM-implicit score, available-case
  uses influence-function decomposition of the moment-based statistic.

## Recommended language for the paper appendix

1. **State (E) explicitly.** "Counts data assumes exchangeable iid
   raters. Rater-specific tendencies cannot be identified from counts;
   for that setting use rater-identified data."
2. **State (S) explicitly.** "Missingness is modelled as uniformly-
   random rater dropout per subject, parameterised by the observed
   total `r_i`."
3. **Be clear about identifiability.** `θ` is not generally identified;
   kappa is. EM converges; the limit may depend on initialisation but
   the resulting kappa estimate does not.
4. **Be clear that this is not 'MAR for arbitrary missingness'.** For
   counts the MCAR-vs-MAR distinction is narrow. The MAR-handling
   power that justifies FIML over IPW in rater-identified data does
   not carry over.

## Open questions for the paper

- Should we report `θ_hat` to the user, or only the derived kappa?
  Reporting `θ_hat` is honest about what EM converged to, but it's
  not identified, which invites confusion. Recommendation: do not
  report `θ_hat` in the public API; document it as an internal
  intermediate.
- How aggressive should pruning be? `prune_tol = 1e-9` from the
  legacy drops vanishing-mass compositions; under-pruning gives a
  near-singular information matrix that the pseudo-inverse handles
  but may give noisy SEs.
- Does the simulation study need a counts-FIML cell? Probably not —
  the paper's MCAR / MAR story is in rater-identified data. Counts-
  FIML is a comparison estimator, mentioned in the appendix.
