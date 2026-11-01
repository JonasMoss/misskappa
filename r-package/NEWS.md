# misskappa 3.0.0

misskappa is now an agreement-only package. Code written for 2.0.0 needs the
renames below.

## Breaking changes

* Coefficient alpha is removed: `alpha()`, `alpha_test()`, the C++ alpha
  estimators, and the `dat.holzinger1939` dataset. The last commit with alpha
  is tagged `alpha-final`.
* The public functions are renamed so the package no longer masks
  `base::kappa()`: `kappa()` is now `agreement()`, `kappa_counts()` is
  `agreement_counts()`, and `kappa_test()` is `agreement_test()`.
* Subjects with no observed rating are dropped by every estimator, so `n` and
  the influence rows count only rated subjects. Before, the normal-theory
  FIMLs dropped them and the other estimators kept them.
* Numeric ratings are matched as labels, so 0.5 and 1.2 are two categories.
  Before, the categorical estimators silently truncated them to integers.

## New features

* `agreement()` gains `categories`, which declares the full rating scale,
  including categories nobody used. The scale sets Brennan-Prediger's chance
  baseline and the range that scales the linear and quadratic weights. Before,
  both came from the observed categories only, so a 1-5 scale on which nobody
  chose 5 was treated as 1-4, and `values = 1:5` was an error.
* `agreement()` accepts character, logical, and factor ratings. A data frame
  whose columns are factors with shared levels takes its scale from the levels.
  Fits record the scale used as `fit$categories`. Text ratings have no order
  of their own, so linear and quadratic weights and the scored estimators
  need `categories` (or factors with shared levels) for them, and
  `ratings_to_counts()` orders its columns by shared factor levels.
* Influence rows are named by subject (row names, or row numbers), and
  `agreement_test(paired = TRUE)` matches fits by subject. A subject unrated
  in one fit contributes nothing to it, so fits with different missing rows
  can be compared. Duplicated row names are an error in paired tests.
* A coefficient whose chance disagreement is zero (every rating in one
  category) is reported as undefined with a warning, and its variance is `NA`
  instead of a misleading zero.

## Documentation

* The help pages are shorter and name the recommended default estimators
  (`"ipw"` for `agreement()`, `"fleiss_cuzick"` for `agreement_counts()`).
* The documentation is consolidated into three pages: the getting-started
  vignette (now including equality tests and custom contrasts), a
  coefficients, estimators, and validation article, and a beyond-pairwise
  article (g-wise kernels and vector-valued ratings).

## Minor

* Error messages no longer show internal function calls.
* `psych` and `coefficientalpha` leave Suggests.

## Internals

* Unexported code is split into plumbing (dot-prefixed helpers the public
  functions call) and research entry points (`research_*()` in
  `R/research.R`, and the C++ namespace `misskappa::research`), which no
  public function reaches and which carry no stability promise. Former
  internal names: `estimate_kappa_raw()` -> `research_kappa_raw()`,
  `estimate_kappa_fiml_multi()` -> `research_kappa_fiml_multi()`,
  `.kappa_pairwise_cs()` -> `research_kappa_pairwise_cs()`,
  `kappa_quadratic()` -> `research_kappa_quadratic()`, `kappa_continuous()`
  -> `research_kappa_continuous()`, `kappa_gwise()` ->
  `research_kappa_gwise()`, `joint_vcov()` -> `research_joint_vcov()`.
* The internal `wald_test()` and the C++-only SQUAREM option of the
  categorical EM are removed. Arbitrary contrasts are built from `fit$psi`
  (see the influence-functions article).
* The simulation experiments and research notes moved out of the public
  repository.

# misskappa 2.0.0

This release tightens the R interface. Several changes break code written for
1.x, so read the first section before upgrading.

## Breaking changes

* Brennan-Prediger is reported only when the rating scale is declared. The
  scored estimators `kappa(estimator = "nt_fiml")` and `"nt_fiml_cs"` report it
  when `values` is supplied (a uniform draw from the declared scores) and, for
  vector-valued ratings, when `ranges` is supplied. Before, they fell back to a
  continuous uniform on the observed score range, so the coefficient depended
  on whether `values` was passed (for codes 1 to 5 the two chance baselines
  differ by a factor of 1.5) and drifted with the sample range. The categorical
  estimators `"ipw"` and `"cat_fiml"` are unchanged.
* `kappa_counts(estimator = "nt_fiml")` is renamed `"nt_fiml_cs"`, so that the
  estimator name carries the covariance structure everywhere (it fits the
  exchangeable compound-symmetry model). It now defaults to
  `weight = "quadratic"`.
* `kappa()` errors on `...` arguments that do not apply to the input, such as
  typos, `W` on a matrix, `values` or `g` on an array, or `ranges` with
  `"pairwise"`. Every estimator also errors on `em_options` entries it does not
  read, including any `em_options` for the moment estimators.
* The `em_options$covariance` and `em_options$mean` shims deprecated in 1.2.0
  are removed (use `estimator = "nt_fiml_cs"`), and so is the no-op
  `em_options$fd_h`.
* `kappa(estimator = "pairwise")` no longer runs the complete-data g-wise
  estimator for non-quadratic `g > 2`. Use `"ipw"`, which is the same estimator
  when no ratings are missing.
* `kappa(estimator = "pairwise_cs")` is no longer offered. The estimator stays
  in the package as an internal function for simulations.
* The `sim` list of simulation helpers is removed.
* `kappa_test()` and `alpha_test()` default to `paired = FALSE` (independent
  samples), as `t.test()` does. `kappa_test()` compares the first coefficient
  of the first fit by default (`"Conger"` for `kappa()`, `"Fleiss"` for
  `kappa_counts()`), and `value` must be 0 when more than two fits are
  compared.
* Alpha fits no longer carry `weight = "score"`, and vector-valued fits record
  the `weight` the user asked for (the component loss stays in `fit$loss`).

## New features

* Complete-data small-sample intervals: supported `kappa()`, `kappa_counts()`,
  and `alpha()` fits carry exact leave-one-subject-out evaluations.
  `confint(method = "jackknife")` forms natural, Fisher, or log-complement
  transformed delete-one t intervals. Exact-boundary kappa uses a one-sided
  disagreement-event/chance fallback, and exact-boundary alpha is explicitly
  unsupported. Wald intervals remain the default, and missing-data inference
  is unchanged.
* Fits record the number of subjects `n`, which `print()` shows. `print()`
  also notes when a normal-theory EM did not converge.
* `as.data.frame()` gains `estimator`, `lower`, and `upper` columns, so fits
  can be stacked with `rbind()`. Extra arguments go to `confint()`.
* `kappa_test()` and `alpha_test()` results carry `null.value` and
  `alternative`, so the printed test states its hypothesis.
* New dataset `dat.tenhove2025`: a *synthetic* planned-incomplete inter-rater
  design (60 subjects, 6 raters, ordinal `1`--`5` scores) whose rater
  co-observation graph is connected but not complete, in the style of the
  planned-incomplete ICC designs of ten Hove et al. (2025). It shows
  `"nt_fiml_cs"` recovering the coefficient on a design where the saturated
  estimators are not identified.

## Documentation

* `?confint.misskappa_estimate` lists where jackknife intervals are available,
  and `as.data.frame()` has its own help page.
* Internal functions no longer have help pages.
* New "Working With Influence Functions" how-to article: the
  `vcov == crossprod(psi) / n^2` identity, and building custom contrasts and
  paired comparisons by hand from `fit$psi` / `influence()`.
* The "Missingness And Estimators" article gains a worked compound-symmetry
  example on `dat.tenhove2025` (saturated non-identified versus CS fit).
* `kappa()`, `alpha()`, and `kappa_counts()` reference the two methods papers
  (in preparation), and `inst/CITATION` lists them as placeholders to be
  updated on preprint and on publication.

## Internals

* The test suite no longer calls the sibling package magmaan. The normal-FIML
  parity tests compare against a pure-R reference EM and sandwich kept in a
  test helper. `utils` joins Imports and `coefficientalpha` joins Suggests.
* The C++ library drops `EmOptions::fd_h` and moves the research-only Louis
  diagnostic and the `build_psi_from_phi()` helper into private headers. Its
  headers now state that the C++ interface is the R package's engine and
  carries no separate stability guarantee.

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
* The full-`W` vector `"nt_fiml"` / `"nt_fiml_cs"` estimators accept a `ranges`
  argument (a features-by-2 matrix of `(lo, hi)` score bounds) that fixes the
  Brennan-Prediger chance baseline `sum_a W_aa (hi_a - lo_a)^2 / 6` instead of
  estimating each width from the observed data range, the vector analogue of the
  scalar `values` hook. Fixing the range makes the chance denominator
  deterministic, so the Brennan-Prediger standard error is exact and the estimand
  no longer drifts with the sample range.

## API

* The covariance structure of the normal-theory FIML is now carried by the
  estimator name: `"nt_fiml"` is the saturated model and `"nt_fiml_cs"` the
  compound-symmetry one. The old `em_options$covariance` and `em_options$mean`
  are deprecated. They warn for one release and are otherwise ignored on
  `kappa()` / `alpha()` (the coefficients are mean-free, so use `kappa_counts()`
  for the exchangeable equal-mean model). The internal `"ten_hove"` covariance
  alias is removed.

## Documentation

* `kappa()`, `alpha()`, and `kappa_counts()` document the compound-symmetry
  estimators, and the getting-started vignette's estimator table gains a
  saturated-versus-compound-symmetry recommendation. Compound symmetry is
  consistent and best-calibrated under MCAR (even when the covariance is
  misspecified) and under correctly-specified MAR, but can be inconsistent
  under MAR when the true covariance departs from compound symmetry (for
  example rater subgroups with different bias). The saturated `"nt_fiml"` is
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

* `kappa()`: weighted agreement coefficients (Conger, Fleiss,
  Brennan-Prediger) for raw ratings, any number of raters, nominal or quadratic
  loss. `estimator` is `"ipw"` or `"cat_fiml"` for categorical ratings, and
  `"pairwise"` or `"nt_fiml"` for the quadratically weighted scored coefficient.
* `kappa_counts()`: counts-format input (subjects × categories). `estimator`
  is `"fleiss_cuzick"` or `"cat_fiml"`.
* `alpha()`: coefficient alpha for scored categorical and continuous item
  batteries. `estimator` is `"pairwise"`, `"cat_fiml"`, or `"nt_fiml"`.

The two FIMLs are distinct: `cat_fiml` (saturated multinomial) and `nt_fiml`
(robust normal-theory). `"pairwise"` is the MCAR moment estimator. `"ipw"` adds
inverse-probability weighting.

## Inference

* `coef()`, `vcov()`, `confint()`, `print()`, `as.data.frame()`, and
  `stats::influence()` methods. Per-subject influence functions are stored as
  `fit$psi`, an n-by-K matrix satisfying `vcov == crossprod(psi) / n^2`.
* `confint(transform = "fisher")` gives a delta-method Wald interval on the
  Fisher `atanh` scale, back-transformed with `tanh` (stays within (-1, 1),
  better small-sample coverage near the boundary).

## Hypothesis tests

* `kappa_test()` / `alpha_test()`: Wald tests that a coefficient is equal
  across fits, returning an `htest`: one-sample, independent two-sample, paired
  (same subjects), and G-way homogeneity.

## Documentation and data

* "Getting started with misskappa" vignette and a "Testing equality of
  agreement coefficients" how-to article, worked on real data.
* Bundled datasets, including the public-domain Holzinger & Swineford (1939)
  battery for always-runnable coefficient-alpha examples.
* Reference site at <https://jonasmoss.github.io/misskappa>.
