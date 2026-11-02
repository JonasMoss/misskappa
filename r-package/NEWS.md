# misskappa 3.0.1

* `agreement(estimator = "cat_fiml")` no longer warns that the saturated
  joint distribution is not uniquely identified. Simulations showed the
  point estimate does not move between likelihood maximizers by any
  practical amount, so the warning fired on most missing-at-random fits
  without meaning anything. The diagnostic stays in the fit as `null_frac`.
* New silent diagnostic for `"cat_fiml"`: `se_cutoff_ratio`, the standard
  error recomputed with `info_rcond` lowered 100-fold divided by the reported
  one. Values well above 1 mean the standard error leans on weakly
  identified directions of a sparse table.
* The default `max_iter` for the EM estimators (`"cat_fiml"`, `"nt_fiml"`,
  count-format `"cat_fiml"`) is now 100000 (was 10000). Sparse fits often
  need more than 10000 iterations to meet the default `tol`, which the
  standard error needs; the estimate itself settles earlier.
* `"nt_fiml"` gains `fallback = "pairwise"` (the default): when the fit fails
  (singular or nearly singular fitted covariance, numerical failure, or no EM
  convergence) `agreement()` returns the `"pairwise"` estimate of the same
  coefficient with a warning and the reason in `fit$fallback`. Argument and
  identification errors are never caught. `fallback = "none"` restores the
  error. Other estimators reject `fallback`.
* `agreement_counts()`: the Fleiss-Cuzick moment estimator is now
  `estimator = "pairwise"` (was `"fleiss_cuzick"`), matching `agreement()`.
* Research entry points are renamed `research_agreement_*()` (were
  `research_kappa_*()`), matching the public names. They stay unexported
  and unsupported.
* `"nt_fiml"` on complete data no longer fails when the sample covariance is
  singular, for example when two raters agree perfectly. It now uses the
  sample-moment influence function directly, which equals the previous
  result whenever the covariance is invertible (and equals `"pairwise"`).
* When `"nt_fiml"` fails on incomplete data because the fitted covariance is
  singular, the error now says so and points to the moment estimators or
  `"cat_fiml"` (was "Numerical error."). Nearly singular fits (smallest over
  largest eigenvalue of the fitted covariance below 3e-5) are refused the
  same way: near that boundary the sandwich standard error explodes (up to
  tens of times the moment estimator's).

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
