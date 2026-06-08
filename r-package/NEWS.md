# misskappa 1.0.0

First release of the rewritten package: a thin Rcpp interface over the
standalone C++17 `misskappa` library.

## Public API

Estimators return a `misskappa_estimate` object carrying named coefficients and
their asymptotic covariance matrix. Each estimator is selected with a single,
unambiguous `estimator=` argument:

* `kappa()` — weighted agreement coefficients (Conger, Fleiss,
  Brennan–Prediger) for raw ratings. `estimator` is one of:
  `"ipw"` (inverse-probability-weighted) or `"cat_fiml"` (saturated-multinomial
  FIML) for categorical ratings; `"pairwise"` (pairwise-available moment) or
  `"nt_fiml"` (robust normal-theory FIML) for the quadratically weighted, scored
  coefficient (these two require `weight = "quadratic"`).
* `kappa_counts()` — counts-format input (subjects × categories);
  `estimator` is `"pairwise"` or `"cat_fiml"`.
* `alpha()` — coefficient alpha; `estimator` is `"pairwise"`, `"cat_fiml"`, or
  `"nt_fiml"`. (Replaces the former `method` + `type` arguments.)

The two FIMLs are named distinctly — `cat_fiml` (saturated multinomial) versus
`nt_fiml` (robust normal-theory) — so the same token never denotes two different
estimators. The unweighted loss is `weight = "nominal"`; the MCAR moment
estimator is `"pairwise"`.

Estimators dropped from the public surface (Gwet's reweighting, available-case
for categorical `kappa`, and g-wise agreement) remain available to simulations
through internal entry points but are no longer exported.

## Inference

* `coef()`, `vcov()`, `confint()`, `print()`, and `as.data.frame()` methods on
  `misskappa_estimate`. Per-subject influence functions are exposed as the
  documented `fit$psi` component (an n-by-K matrix satisfying
  `vcov == crossprod(psi) / n^2`) rather than via an `influence()` method.
* `confint(transform = "fisher")` returns a delta-method Wald interval on the
  Fisher `atanh` scale, back-transformed with `tanh` (stays within
  (-1, 1); better small-sample coverage near the upper boundary). `print()`
  now reports confidence limits alongside the estimate and SE.

## Hypothesis tests

* `kappa_test()` / `alpha_test()` — Wald tests that a coefficient is equal
  across fits, returning an `htest`. They cover one-sample, independent
  two-sample, paired (same-subject; the dependence is taken from the `fit$psi`
  influence functions), and G-way homogeneity, and are built only on the
  exported `coef`/`vcov`/`fit$psi`.
* The general linear-hypothesis engine they stand on — `joint_vcov()` and
  `wald_test()` (arbitrary contrasts and non-zero margins) — and the closed-form
  `kappa_quadratic()` remain internal, for the rare power-user case.

## Documentation and datasets

* New dataset `dat.holzinger1939`: the classic, public-domain Holzinger &
  Swineford (1939) nine-item mental-ability battery — a compact *continuous*
  item battery used for the always-runnable coefficient-alpha examples.
* Every exported function and dataset now carries a runnable example built on
  the bundled data, including paired and independent `kappa_test()` /
  `alpha_test()` calls. The legacy `dat.gwet2014` / `dat.klein2018` /
  `dat.zapf2016` help pages gained `format` and `source` details.
* The coefficient-alpha equality tests are additionally demonstrated on **real
  missing data with grouping structure** via `psych::bfi` (added to `Suggests`):
  a dependent test across two Big Five scales (same respondents), an independent
  test across gender, and G-way homogeneity across all five scales.
* New shipped vignette *"Getting started with misskappa"* (a load → fit →
  confint → test walkthrough), plus a *"Testing equality of agreement
  coefficients"* how-to article covering the one-sample, independent, paired,
  and G-way cases on real data.
