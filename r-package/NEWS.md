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
