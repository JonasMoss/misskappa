# misskappa (development version)

## Raw cat_fiml: identification guard relaxed, flattening added

* The saturated raw-data FIML estimators (`kappa(estimator = "cat_fiml")`,
  `alpha(estimator = "cat_fiml")`, and the g-wise variant) no longer hard-fail
  when the saturated joint distribution is not uniquely identified but the
  coefficient itself is. With every rater pair co-observed the coefficients
  are estimable functions of the identified pattern margins, so the fit now
  succeeds and reports a `null_frac` diagnostic (per-coefficient fraction of
  the delta-method gradient in the truncated null space of the Louis
  information), warning when it exceeds 0.01. The design-level guard — every
  rater pair must be co-observed by at least one subject — still errors.
* New `em_options$flatten` for the raw categorical FIML: total Dirichlet
  pseudo-mass spread over the complete pattern table. Any positive value
  makes the fitted table the unique interior posterior mode (the analytic
  center of the flat maximum-likelihood face), shrinking it toward uniform
  with weight `flatten / (n + flatten)`. Default `0` keeps strict ML and the
  legacy deterministic-start behaviour.

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
