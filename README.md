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

A C++23 library for agreement coefficients (Cohen / Fleiss / Brennan-Prediger)
with arbitrary numbers of raters, arbitrary pairwise loss matrices, and
missing-data support under MCAR (available-case, IPW) and MAR (FIML / EM).
Gwet's estimator is retained for comparison.

The repository ships three things:

- A standalone C++23 static library (`include/misskappa/`, `src/`).
- An R package wrapping it via Rcpp (`r-package/`).
- Research manuscripts using both (under `papers/`, with the lead manuscript targeting Psychometrika).

The original C++17 + Armadillo implementation that grew alongside the
manuscript is preserved under `dev/legacy/` for reference.

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

See `AGENTS.md` for the project contract and each paper's own
`papers/<slug>/AGENTS.md` for manuscript-specific direction (the three
spinoffs are under `papers/{ipw,fiml,quadratic}/`).

## License

This repository is dual-licensed:

- **Code** — the C++ library (`include/`, `src/`), the R package (`r-package/`),
  and the build tooling are released under the **MIT License** (see `LICENSE`).
- **Papers** — the manuscripts, figures, and text under `papers/` are released
  under **CC BY 4.0** (see `papers/LICENSE`), in line with Plan S / the
  Norwegian rights-retention strategy.
