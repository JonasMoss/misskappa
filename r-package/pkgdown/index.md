# misskappa <img src="man/figures/logo.png" align="right" height="110" alt="" />

Estimation and inference for weighted agreement coefficients (Cohen, Fleiss,
Conger, Brennan-Prediger) with arbitrary numbers of raters and arbitrary
pairwise loss functions, plus coefficient alpha for scored categorical and
continuous item batteries. Supports incomplete ratings under MCAR
(available-case, IPW) and MAR (categorical or normal FIML / EM). Wraps the
standalone C++23 `misskappa` library.

The main entry points are `kappa()` for raw categorical rating matrices and
`alpha()` for numeric item batteries. Specialized helpers cover saturated
categorical alpha FIML (`alpha_cat_fiml()`), counts-format input, continuous
ratings, closed-form quadratic loss, and complete rectangular g-wise agreement.

See the [function reference](reference/index.html) for estimators and inference
helpers, and the [mathematical guides](articles/index.html) for the loss-matrix
formulation, missing-data estimators, and validation strategy. The underlying
C++ library has its own [C++ API reference](cpp/).
