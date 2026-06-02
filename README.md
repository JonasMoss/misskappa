<img src="docs/logo.png" alt="misskappa" align="right" width="300">

# misskappa

[![CI](https://github.com/JonasMoss/misskappa/actions/workflows/ci.yml/badge.svg)](https://github.com/JonasMoss/misskappa/actions/workflows/ci.yml)
![C++23](https://img.shields.io/badge/C%2B%2B-23-00599C?logo=cplusplus&logoColor=white)
![R package](https://img.shields.io/badge/R-package-276DC3?logo=r&logoColor=white)
![Code: MIT](https://img.shields.io/badge/code-MIT-blue)
![Papers: CC BY 4.0](https://img.shields.io/badge/papers-CC_BY_4.0-lightgrey)
<br>
![κ: substantial agreement](https://img.shields.io/badge/%CE%BA-substantial%20agreement-brightgreen)
![raters: ∞](https://img.shields.io/badge/raters-%E2%88%9E-blueviolet)
![missing data: handled](https://img.shields.io/badge/missing%20data-handled%20%F0%9F%98%8E-success)
![Kappa](https://img.shields.io/badge/Kappa-%F0%9F%98%8F-9146FF?logo=twitch&logoColor=white)

A C++23 library for weighted agreement coefficients with any number of raters,
arbitrary pairwise loss matrices, and support for missing ratings.

## Coefficients

- **Cohen's kappa** — [Cohen 1960](https://doi.org/10.1177/001316446002000104), [1968](https://doi.org/10.1037/h0026256)
- **Fleiss' kappa** — [Fleiss 1971](https://doi.org/10.1037/h0031619)
- **Conger's kappa** — [Conger 1980](https://doi.org/10.1037/0033-2909.88.2.322)
- **Brennan–Prediger** — [Brennan & Prediger 1981](https://doi.org/10.1177/001316448104100307)

## Missing data

Incomplete categorical ratings are handled under two mechanisms:

- **MCAR** — available-case analysis and inverse-probability weighting
- **MAR** — full-information maximum likelihood and EM

## Contents

- **C++23 library** — `include/misskappa/` and `src/`
- **R package** — `r-package/`, wrapping the library via Rcpp
- **Manuscripts** — `papers/`, with the lead manuscript targeting Psychometrika

## Quick build

```sh
cmake --preset dev
cmake --build --preset dev
ctest --preset dev
```

Or via `just`:

```sh
just test            # dev build + ctest
just r-install       # build opt + install R package
just r-check         # R CMD check + testthat
just paper <slug>    # build the manuscript PDF for papers/<slug>/
```

See `AGENTS.md` for the project contract and each paper's `AGENTS.md` for
manuscript-specific direction.

## License

This repository is dual-licensed:

- **Code** under the MIT License — see `LICENSE`. Covers the C++ library, the
  R package, and the build tooling.
- **Papers** under CC BY 4.0 — see `papers/LICENSE`. Covers the manuscripts,
  figures, and text in `papers/`, in line with Plan S and the Norwegian
  rights-retention strategy.
