# 36 — Is exchangeability-constrained cat-FIML the count-FIML?

**Question.** The saturated cat-FIML models the joint multinomial over
`{1..C}^R` for identified raters and integrates out missing raters. If we
impose rater exchangeability on that model, do we get the count-FIML of
`kappa_counts(estimator = "cat_fiml")`, or a different estimator that would
need its own implementation (with SEs)?

**Answer: exactly the count-FIML. Nothing new to implement.**

Argument (verified numerically in `check.R`):

1. Exchangeability means the joint pmf is constant on each permutation orbit
   of `{1..C}^R`; orbits are indexed by the count composition `z`
   (`sum z_c = R`), so the free parameter is `theta` on the composition
   simplex — the same parameter space the C++ backend
   (`estimate_fiml_counts.cpp`) uses.
2. For a subject with observed rater subset `S` (`|S| = m`, observed counts
   `c`), the rater-identified observed-data likelihood is
   `P(Y_S = y_S) = sum_z theta_z * multinom(R - m; z - c) / multinom(R; z)`,
   which is **proportional in theta** to the multivariate-hypergeometric
   counts likelihood `sum_z theta_z * prod_c choose(z_c, c_c) / choose(R, m)`
   that the backend maximises. MAR ignorability disposes of the missingness
   factor as usual (the mechanism may be rater-specific); the "uniformly
   random subsample" reading of the hypergeometric is a *consequence* of
   exchangeability, not an extra assumption.
3. Same likelihood up to a theta-free constant ⇒ same MLE, same observed
   information, same delta-method (Fleiss, BP) vcov.

`check.R` (R = 4, C = 3, n = 300, rater-specific missingness up to 55%):

- log-lik difference over 25 random theta: spread 4e-12 (constant).
- direct raw-data ML vs `kappa_counts(cat_fiml)`: Fleiss/BP agree to ~1e-7.
- observed-information + delta-method SEs vs the package's
  Louis-information SEs: agree to ~2e-6.

## R exposure (decided 2026-06-12)

The estimator is already shipped as `kappa_counts(estimator = "cat_fiml")`;
the only gap was the data-shape bridge from a rater-identified
subjects-by-raters matrix. Exposed as a **converter**,
`ratings_to_counts(x, categories = NULL)` (`R/kappa.R`, exported, man page +
tests), rather than a `kappa()` `estimator=` flag or a `kappa_exchangeable()`
wrapper. Rationale: minimal surface, composes with the existing counts entry
point, makes the *discard of rater identity* (= the exchangeability
assumption) visible in user code, also unlocks `fleiss_cuzick` from rater
data, and avoids conflating the fixed-rater (Conger/Fleiss) and exchangeable
(Fleiss/BP) estimands under one dispatcher.

**No `r_total` footgun (verified).** Fleiss and Brennan-Prediger are pairwise
functionals, so the count-FIML's estimate *and* SE are invariant to `r_total`
to 6-7 decimals wherever the EM converges; `r_total` above the max observed
row sum only makes the EM completion space larger and risks non-convergence
without changing the answer. So `kappa_counts`'s default
`r_total = max(rowSums)` is the correct, cheapest, most stable choice, and the
converter needs no `r_total` handling.

**Caveats.** The equivalence is a statement about the *exchangeable* model.
If raters are genuinely non-exchangeable (rater-specific marginals /
confusion), the count-FIML is misspecified and rater-specific MAR can bias
it — that is precisely the fixed-rater vs exchangeable regime split already
built into study 29. Also, `kappa_counts` consumes counts, so Conger (which
needs rater identity) is out of scope by construction; under exchangeability
the Conger and Fleiss estimands coincide anyway.
