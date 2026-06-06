# misskappa 2.0.0

First release of the rewritten package: a thin Rcpp interface over the
standalone C++23 `misskappa` library.

## Public API

Estimators return a `misskappa_estimate` object carrying named coefficients and
their asymptotic covariance matrix, keyed by input shape:

* `kappa()` — weighted agreement coefficients (Cohen/Conger, Fleiss,
  Brennan–Prediger) for raw categorical ratings; `method` is one of
  `"available"`, `"ipw"`, `"fiml"`, `"gwet"`.
* `kappa_counts()` — counts-format input (subjects × categories).
* `kappa_continuous()` — continuous ratings.
* `kappa_gwise()` — complete rectangular g-wise agreement.
* `alpha()` — coefficient alpha. `method = "available"` uses pairwise-available
  covariance; `method = "fiml"` selects a saturated EM fit via
  `type = "normal"` (Gaussian) or `type = "categorical"` (multinomial). This is
  the single alpha entry point.

## Inference

* `coef()`, `vcov()`, `confint()`, `print()`, `as.data.frame()`, and
  `influence()` methods on `misskappa_estimate`.
* `confint(transform = "fisher")` returns a delta-method Wald interval on the
  Fisher `atanh` scale, back-transformed with `tanh` (stays within
  (-1, 1); better small-sample coverage near the upper boundary). `print()`
  now reports confidence limits alongside the estimate and SE.

## Notes

* Hypothesis tests are deferred to a later release. The IF-based joint
  covariance / Wald-contrast engine (`joint_vcov()`, `wald_test()`) and the
  closed-form `kappa_quadratic()` estimator remain in the package as internal
  functions.
