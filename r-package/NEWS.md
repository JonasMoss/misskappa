# misskappa 3.0.1

* `agreement(estimator = "cat_fiml")` no longer warns that the saturated
  joint distribution is not uniquely identified. Simulations showed the
  point estimate does not move between likelihood maximizers by any
  practical amount, so the warning fired on most missing-at-random fits
  without meaning anything. The diagnostic stays in the fit as `null_frac`.
* New SE-stability check for `"cat_fiml"`: fits carry `se_cutoff_ratio`, the
  standard error recomputed with `info_rcond` lowered 100-fold divided by the
  reported one, and warn when it exceeds 1.25. This flags standard errors
  that lean on weakly identified directions of a sparse saturated table,
  which `null_frac` did not detect. The warning no longer recommends
  `flatten`.

# misskappa 3.0.0

First CRAN release.

* `agreement()` estimates Conger's (Cohen's for two raters), Fleiss', and
  Brennan-Prediger's agreement coefficients from a subjects-by-raters matrix
  with missing ratings. Estimators: inverse-probability weighting (`"ipw"`,
  the default, valid under MCAR), categorical FIML (`"cat_fiml"`, MAR), and
  for quadratic weights a pairwise moment estimator and robust normal-theory
  FIML, including a compound-symmetry variant for planned-incomplete designs.
* `agreement_counts()` handles counts-format data (subjects by categories),
  and `ratings_to_counts()` converts rater-identified data to that format.
* Ratings may be numeric, text, or factors; `categories` declares the full
  rating scale, including categories nobody used.
* Every fit carries influence-function standard errors. `confint()` gives
  Wald intervals (optionally Fisher-transformed) and exact jackknife
  intervals on complete data; `agreement_test()` tests equality of a
  coefficient across independent or paired fits.
* Multirater g-wise kernels (`g > 2`) and vector-valued ratings are
  supported.

Versions 1.0.0 to 2.0.0 were released on GitHub only; their history is in
the repository.
