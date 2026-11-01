

<!-- README.md is generated from README.qmd. Please edit that file and run `quarto render README.qmd`. -->

<img src="r-package/man/figures/logo.png" alt="misskappa" align="right" width="192">

# misskappa

[![CI](https://github.com/JonasMoss/misskappa/actions/workflows/ci.yml/badge.svg)](https://github.com/JonasMoss/misskappa/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/JonasMoss/misskappa/graph/badge.svg)](https://codecov.io/gh/JonasMoss/misskappa)
![Code: MIT](https://img.shields.io/badge/code-MIT-blue)

An R package, backed by a C++17 library, for weighted agreement coefficients
with any number of raters and missing ratings: Cohen’s kappa
([1960](https://doi.org/10.1177/001316446002000104),
[1968](https://doi.org/10.1037/h0026256)), Conger’s kappa
([1980](https://doi.org/10.1037/0033-2909.88.2.322)), Fleiss’ kappa
([1971](https://doi.org/10.1037/h0031619)), and the Brennan–Prediger
coefficient ([1981](https://doi.org/10.1177/001316448104100307)).

Missing ratings are handled by inverse-probability weighting (`"ipw"`,
consistent under MCAR), categorical full-information maximum likelihood
(`"cat_fiml"`, efficient under MAR), and, for quadratic weights, a pairwise
moment estimator (`"pairwise"`, MCAR) and robust normal-theory FIML
(`"nt_fiml"`, [Yuan & Bentler, 2000](https://doi.org/10.1111/0081-1750.00078)).
Every estimate comes with a standard error.

## Installation

The package builds from source, so it needs a C++17 toolchain (Rtools on
Windows, the Command Line Tools on macOS).

``` r
# install.packages("remotes")
remotes::install_github("JonasMoss/misskappa", subdir = "r-package")
```

## Usage

``` r
library(misskappa)
```

Five raters sort subjects into three categories, and some ratings are missing
([Klein, 2018](https://doi.org/10.1177/1536867X1801800408)):

``` r
agreement(dat.klein2018, estimator = "ipw")
#> misskappa: estimator = ipw, weight = nominal, n = 10
#>                  estimate     se  lower  upper
#> Conger             0.4301 0.1050 0.2244 0.6358
#> Fleiss             0.4054 0.1196 0.1710 0.6399
#> Brennan-Prediger   0.4204 0.1136 0.1978 0.6430
```

Ratings may be numbers, strings, or factors. Pass `weight = "linear"` or
`"quadratic"` for ordered categories, and `categories = 1:5` to declare the
full rating scale when some categories went unused.

`agreement_test()` compares coefficients across fits. These two studies have
different subjects, so the samples are independent (use `paired = TRUE` for
fits on the same subjects):

``` r
agreement_test(
  gwet  = agreement(dat.gwet2014,  estimator = "ipw"),
  klein = agreement(dat.klein2018, estimator = "ipw")
)
#> 
#>  Independent-sample Wald test of equal Conger across 2 fits
#> 
#> data:  gwet, klein
#> X-squared = 0.0019646, df = 1, p-value = 0.9646
#> alternative hypothesis: true difference in Conger is not equal to 0
#> sample estimates:
#>      gwet     klein 
#> 0.4233051 0.4301049
```

The package also handles counts-format data (`agreement_counts()`), *g*-wise
kernels that compare `g > 2` raters at once, and vector-valued ratings. Fits
work with `coef()`, `vcov()`, `confint()`, and `as.data.frame()`.

## Documentation

The [package website](https://jonasmoss.github.io/misskappa/) has the
reference and articles:

- [Getting started](https://jonasmoss.github.io/misskappa/articles/misskappa.html)
- [Missingness and estimators](https://jonasmoss.github.io/misskappa/articles/missingness-estimators.html)
- [Agreement coefficients and loss matrices](https://jonasmoss.github.io/misskappa/articles/agreement-coefficients.html)
- [Testing equality of agreement coefficients](https://jonasmoss.github.io/misskappa/articles/equality-tests.html)
- [Working with influence functions](https://jonasmoss.github.io/misskappa/articles/influence-functions.html)
- [g-wise agreement](https://jonasmoss.github.io/misskappa/articles/gwise-agreement.html)
- [Vector-valued agreement](https://jonasmoss.github.io/misskappa/articles/vector-valued-agreement.html) *(in preparation)*
- [How misskappa is validated](https://jonasmoss.github.io/misskappa/articles/validation.html)
- [C++ API reference](https://jonasmoss.github.io/misskappa/cpp/)

## Development

The C++ library (`include/misskappa/`, `src/`) is the engine of the R package
(`r-package/`) and has no separate stability promise. `papers/` holds the
completed papers that accompany the package. `AGENTS.md` has the project
conventions.

``` sh
just test        # C++ dev build + tests
just r-install   # vendor the C++ into the R package and install it
just r-check     # fast R CMD check (no vignettes)
```

## License

Code under the MIT License, papers under CC BY 4.0.
