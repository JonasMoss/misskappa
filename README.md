

<!-- README.md is generated from README.qmd. Please edit that file and run `quarto render README.qmd`. -->

<img src="docs/logo.png" alt="misskappa" align="right" width="300">

# misskappa

[![CI](https://github.com/JonasMoss/misskappa/actions/workflows/ci.yml/badge.svg)](https://github.com/JonasMoss/misskappa/actions/workflows/ci.yml)
![Code: MIT](https://img.shields.io/badge/code-MIT-blue)
![Papers: CC BY 4.0](https://img.shields.io/badge/papers-CC_BY_4.0-lightgrey)
![missing data: handled](https://img.shields.io/badge/missing%20data-handled%20%F0%9F%98%8E-success)

A C++17 library and R package for weighted agreement coefficients with any number of raters,
arbitrary pairwise loss matrices, and support for missing ratings.

We support weighted variants of Cohen’s kappa ([1960](https://doi.org/10.1177/001316446002000104), [1968](https://doi.org/10.1037/h0026256)), its multirater variant Conger’s kappa ([1980](https://doi.org/10.1037/0033-2909.88.2.322)) Fleiss’ kappa ([1971](https://doi.org/10.1037/h0031619)) and the Brennan–Prediger coefficient ([1981](https://doi.org/10.1177/001316448104100307)).

Missing data is handled by several estimators. For categorical data, `Cat-FIML` is efficient under MAR and MCAR assumptions. For quadratic weights, pairwise is consistent under MCAR, while robust `NT-FIML` ([Yuan & Bentler, 2000](https://doi.org/10.1111/0081-1750.00078)) is consistent under MCAR and efficient under MCAR and MAR under normality. For general weights and arbitrary data, IPW is consistent under MCAR.

In addition to agreement coefficients, misskappa estimates Cronbach’s coefficient alpha ([1951](https://doi.org/10.1007/BF02310555)) under missing data, with pairwise-available, robust `NT-FIML` ([Zhang & Yuan, 2016](https://doi.org/10.1177/0013164415594658)), and categorical `Cat-FIML` estimators.

## Installation

Install the R package from GitHub with [remotes](https://remotes.r-lib.org/).
It builds from source, so a C++17 toolchain is required (Rtools on Windows, the
Command Line Tools on macOS).

``` r
# install.packages("remotes")
remotes::install_github("JonasMoss/misskappa", subdir = "r-package")
```

## Usage (R)

``` r
library(misskappa)
```

Weighted kappa for five raters and three categories, with some ratings
missing ([Klein, 2018](https://doi.org/10.1177/1536867X1801800408)). The IPW
estimator is consistent under MCAR:

``` r
kappa(dat.klein2018, estimator = "ipw")
#> misskappa: estimator=ipw, weight=nominal
#>                  estimate     se  lower  upper
#> Conger             0.4301 0.1050 0.2244 0.6358
#> Fleiss             0.4054 0.1196 0.1710 0.6399
#> Brennan-Prediger   0.4204 0.1136 0.1978 0.6430
```

Coefficient alpha under missing data, on the Neuroticism items of
[`psych::bfi`](https://CRAN.R-project.org/package=psych) — 2800 respondents
from the [SAPA project](https://www.sapa-project.org/) ([Revelle, Wilt &
Rosenthal, 2010](https://doi.org/10.1007/978-1-4419-1210-7_2)), 106 with at
least one missing answer. Robust normal-theory FIML uses every respondent:

``` r
data(bfi, package = "psych")
N <- paste0("N", 1:5)
alpha(as.matrix(bfi[, N]), estimator = "nt_fiml")
#> misskappa: estimator=nt_fiml, weight=score
#>       estimate    se lower  upper
#> alpha   0.8138 0.006 0.802 0.8256
```

`alpha_test()` compares reliabilities. Is the Neuroticism scale equally
reliable for men and women? The two groups are different respondents, so the
samples are independent (`paired = FALSE`):

``` r
g <- split(seq_len(nrow(bfi)), bfi$gender)
alpha_test(
  men   = alpha(as.matrix(bfi[g[["1"]], N]), estimator = "nt_fiml"),
  women = alpha(as.matrix(bfi[g[["2"]], N]), estimator = "nt_fiml"),
  paired = FALSE)
#> 
#>  Independent-sample test of equal alpha across 2 fits
#> 
#> data:  men, women
#> X-squared = 3.3766, df = 1, p-value = 0.06613
#> sample estimates:
#>       men     women 
#> 0.7959317 0.8209147
```

Each result is a `misskappa_estimate` object with the estimates and their
covariance matrix. Extract them with `coef()`, `vcov()`, and `confint()`.

## Documentation

The [package website](https://jonasmoss.github.io/misskappa/) has the full
reference and the worked-example articles:

- [Getting started](https://jonasmoss.github.io/misskappa/articles/misskappa.html)
- [Agreement coefficients and loss matrices](https://jonasmoss.github.io/misskappa/articles/agreement-coefficients.html)
- [Missingness and estimators](https://jonasmoss.github.io/misskappa/articles/missingness-estimators.html)
- [Testing equality of agreement coefficients](https://jonasmoss.github.io/misskappa/articles/equality-tests.html)
- [Validation strategy](https://jonasmoss.github.io/misskappa/articles/validation.html)
- [C++ API reference](https://jonasmoss.github.io/misskappa/cpp/)

## Contents

- **C++17 library** in `include/misskappa/` and `src/`
- **R package** in `r-package/`, wrapping the library via Rcpp
- **Manuscripts** in `papers/`, WIP papers building on this package

## Quick build

``` sh
cmake --preset dev
cmake --build --preset dev
ctest --preset dev
```

Or via `just`:

``` sh
just test            # dev build + ctest
just r-install       # build opt + install R package
just r-check         # R CMD check + testthat
just paper <slug>    # build the manuscript PDF for papers/<slug>/
```

See `AGENTS.md` for the project contract and each paper’s `AGENTS.md` for
manuscript-specific direction.

## License

**Code** under the MIT License, **Papers** under CC BY 4.
