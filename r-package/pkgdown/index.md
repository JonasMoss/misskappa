# misskappa <img src="man/figures/logo.png" align="right" height="110" alt="" />

Estimation and inference for weighted agreement coefficients (Cohen, Fleiss,
Conger, Brennan-Prediger) with arbitrary numbers of raters and arbitrary
pairwise loss functions. Supports incomplete categorical ratings under MCAR
(available-case, IPW) and MAR (FIML / EM). Wraps the standalone C++23
`misskappa` library.

See the [function reference](reference/index.html) for the estimators and
inference helpers, and the [mathematical guides](articles/index.html) for the
loss-matrix formulation, the missing-data estimators, and the validation
strategy. The underlying C++ library has its own [C++ API reference](cpp/).
