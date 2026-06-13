# misskappa 1.1.0

## Breaking changes

* `kappa_counts()` now selects the count moment estimator with
  `estimator = "fleiss_cuzick"`. The old `"pairwise"` name is gone, though it
  stays valid for `kappa()` and `alpha()`. Count moments now use Fleiss-Cuzick
  unequal-judges weighting, so estimates and standard errors change for
  subjects with unequal rater counts. Balanced data is unchanged.

## New features

* New exported `ratings_to_counts()` collapses a rater-identified
  subjects-by-raters matrix to counts format. Feed it into
  `kappa_counts(estimator = "cat_fiml")` for the exchangeable count-FIML.

## Estimation and inference

* Identification guards across the missing-data estimators. The design guard
  (every rater pair co-observed by at least one subject) errors clearly.
* Saturated `cat_fiml` no longer hard-fails when the joint distribution is
  unidentified but the coefficient is. It reports a `null_frac` diagnostic and
  warns when that exceeds 0.01.
* New `em_options$flatten` for the raw categorical FIML gives a unique,
  start-independent fit. It leaves point estimates essentially unchanged but
  widens the standard errors, so strict ML (the default) stays recommended for
  inference.

## Performance

* The robust normal-theory FIML (`nt_fiml`) runs on a new C++ backend with a
  faster information pass. Results are unchanged.

## Bug fixes

* Declare `stats` in `Imports` and import the `influence` generic, so the
  `stats::influence()` method registers and the namespace loads cleanly.

# misskappa 1.0.0

Initial release: estimation and inference for weighted agreement coefficients
and coefficient alpha under missing data.

## Estimators

Each estimator returns a `misskappa_estimate` carrying named coefficients and
their asymptotic covariance matrix, selected with a single `estimator=`
argument.

* `kappa()` — weighted agreement coefficients (Conger, Fleiss,
  Brennan–Prediger) for raw ratings, any number of raters, nominal or quadratic
  loss. `estimator` is `"ipw"` or `"cat_fiml"` for categorical ratings;
  `"pairwise"` or `"nt_fiml"` for the quadratically weighted scored coefficient.
* `kappa_counts()` — counts-format input (subjects × categories); `estimator`
  is `"fleiss_cuzick"` or `"cat_fiml"`.
* `alpha()` — coefficient alpha for scored categorical and continuous item
  batteries; `estimator` is `"pairwise"`, `"cat_fiml"`, or `"nt_fiml"`.

The two FIMLs are distinct: `cat_fiml` (saturated multinomial) and `nt_fiml`
(robust normal-theory). `"pairwise"` is the MCAR moment estimator; `"ipw"` adds
inverse-probability weighting.

## Inference

* `coef()`, `vcov()`, `confint()`, `print()`, `as.data.frame()`, and
  `stats::influence()` methods. Per-subject influence functions are stored as
  `fit$psi`, an n-by-K matrix satisfying `vcov == crossprod(psi) / n^2`.
* `confint(transform = "fisher")` gives a delta-method Wald interval on the
  Fisher `atanh` scale, back-transformed with `tanh` (stays within (-1, 1),
  better small-sample coverage near the boundary).

## Hypothesis tests

* `kappa_test()` / `alpha_test()` — Wald tests that a coefficient is equal
  across fits, returning an `htest`: one-sample, independent two-sample, paired
  (same subjects), and G-way homogeneity.

## Documentation and data

* "Getting started with misskappa" vignette and a "Testing equality of
  agreement coefficients" how-to article, worked on real data.
* Bundled datasets, including the public-domain Holzinger & Swineford (1939)
  battery for always-runnable coefficient-alpha examples.
* Reference site at <https://jonasmoss.github.io/misskappa>.
