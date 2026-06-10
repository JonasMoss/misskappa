# misskappa <img src="man/figures/logo.png" align="right" height="110" alt="" />

Estimation and inference for weighted agreement coefficients (Cohen, Fleiss,
Conger, Brennan-Prediger) with arbitrary numbers of raters and arbitrary
pairwise loss functions, plus coefficient alpha for scored categorical and
continuous item batteries. Supports incomplete ratings under MCAR and MAR.
Wraps the standalone C++17 `misskappa` library.

Each estimator is selected with a single `estimator=` argument:

* `kappa()` — raw rating matrices. `"ipw"` (inverse-probability-weighted) and
  `"cat_fiml"` (saturated-multinomial FIML) for categorical ratings; `"pairwise"`
  (pairwise-available moment) and `"nt_fiml"` (robust normal-theory FIML) for the
  quadratically weighted, scored coefficient.
* `alpha()` — item batteries. `"pairwise"`, `"cat_fiml"`, or `"nt_fiml"`.
* `kappa_counts()` — counts-format input (subjects × categories).
  `"fleiss_cuzick"` or `"cat_fiml"`.

The MCAR estimators (`"pairwise"`, `"ipw"`) are distribution-free; the FIML
estimators (`"cat_fiml"`, `"nt_fiml"`) are valid under ignorable missingness.
Coefficients carry a covariance matrix, so `coef()`, `vcov()`, and `confint()`
(with an optional Fisher `transform`) give estimates and Wald confidence
intervals. To test whether a coefficient is equal across fits — two groups, two
estimators, two timepoints, or several rater pairs — use `kappa_test()` /
`alpha_test()` (with `paired = TRUE` for same-subject fits); both return a
standard `htest`.

See the [function reference](reference/index.html) for estimators and inference
helpers, and the [mathematical guides](articles/index.html) for the loss-matrix
formulation, missing-data estimators, and validation strategy. The underlying
C++ library has its own [C++ API reference](cpp/).
