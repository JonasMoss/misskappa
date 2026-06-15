# misskappa 1.2.0

## New features

* New `estimator = "nt_fiml_cs"` for `kappa()` and `alpha()`: normal-theory
  FIML under the compound-symmetry (variance-component) covariance model
  `Sigma = psi J + e I`, the fixed-rater structure `Y_ij = mu_j + b_i + eps_ij`.
  It estimates free rater (item) means and is identified by a *connected* item
  co-observation graph, whereas the saturated `"nt_fiml"` needs every item pair
  jointly observed. The point estimate matches the fixed-rater maximum
  likelihood (`lmer(value ~ 0 + item + (1 | subject), REML = FALSE)`) and a
  one-factor SEM, with the usual sandwich standard error.
* New `estimator = "nt_fiml"` for `kappa_counts()` (with `weight = "quadratic"`):
  the exchangeable common-mean compound-symmetry quadratic kappa, fit directly on
  the per-subject count sufficient statistics without realizing the ratings
  table. It is `r_total`-invariant and equals the categorical-FIML exchangeable
  fit.
* New `estimator = "pairwise_cs"` for `kappa()`: the distribution-free MCAR
  moment companion to `"nt_fiml_cs"`. It rater-centers by the available column
  means and fits the compound-symmetry components from pooled within-subject
  pairwise moments, with an analytic influence-function sandwich standard error.
* `kappa()` exposes the full feature-metric quadratic vector estimators on the
  array (subjects-by-raters-by-features) surface: pass a full feature metric `W`
  to `estimator = "pairwise"`, `"nt_fiml"`, or `"nt_fiml_cs"` for the squared
  vector loss `(x - y)' W (x - y)`. `"nt_fiml_cs"` fits the block
  compound-symmetry covariance `Sigma = J_R %x% Psi + I_R %x% E`, the structured
  counterpart to the saturated `"nt_fiml"`. These report Conger, Fleiss, and
  Brennan-Prediger. Without `W`, `"pairwise"` / `"ipw"` stay the
  component-separable (diagonal-weight) estimators as before.

## API

* The covariance structure of the normal-theory FIML is now carried by the
  estimator name: `"nt_fiml"` is the saturated model and `"nt_fiml_cs"` the
  compound-symmetry one. The old `em_options$covariance` and `em_options$mean`
  are deprecated -- they warn for one release and are otherwise ignored on
  `kappa()` / `alpha()` (the coefficients are mean-free; use `kappa_counts()`
  for the exchangeable equal-mean model). The internal `"ten_hove"` covariance
  alias is removed.

## Documentation

* `kappa()`, `alpha()`, and `kappa_counts()` document the compound-symmetry
  estimators, and the getting-started vignette's estimator table gains a
  saturated-versus-compound-symmetry recommendation. Compound symmetry is
  consistent and best-calibrated under MCAR (even when the covariance is
  misspecified) and under correctly-specified MAR, but can be inconsistent
  under MAR when the true covariance departs from compound symmetry (for
  example rater subgroups with different bias); the saturated `"nt_fiml"` is
  the robust choice there.

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
