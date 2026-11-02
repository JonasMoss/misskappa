# misskappa 3.0.0

First CRAN release.

* `agreement()` estimates Conger's (Cohen's for two raters), Fleiss', and
  Brennan-Prediger's agreement coefficients from a subjects-by-raters matrix
  with missing ratings. Estimators: inverse-probability weighting (`"ipw"`,
  the default, valid under MCAR), categorical FIML (`"cat_fiml"`, MAR), and
  for quadratic weights a pairwise moment estimator (`"pairwise"`) and robust
  normal-theory FIML (`"nt_fiml"`), including a compound-symmetry variant
  for planned-incomplete designs (`"nt_fiml_cs"`).
* `"nt_fiml"` falls back to `"pairwise"` by default (`fallback`) when its fit
  fails, for example when the fitted covariance is singular or nearly so
  because two raters agree perfectly; the fit records the reason. On
  complete data the two estimators coincide.
* `agreement_counts()` handles counts-format data (subjects by categories),
  with the Fleiss-Cuzick moment estimator (`"pairwise"`), categorical FIML,
  and compound-symmetry normal FIML; `ratings_to_counts()` converts
  rater-identified data to that format.
* Ratings may be numeric, text, or factors; `categories` declares the full
  rating scale, including categories nobody used.
* Every fit carries influence-function standard errors. `confint()` gives
  Wald intervals (optionally Fisher-transformed) and exact jackknife
  intervals on complete data; `agreement_test()` tests equality of a
  coefficient across independent or paired fits. Categorical FIML fits also
  carry two diagnostics, `null_frac` and `se_cutoff_ratio`.
* Multirater g-wise kernels (`g > 2`) and vector-valued ratings are
  supported.

Versions 1.0.0 to 2.0.0 were released on GitHub only; their history is in
the repository.
