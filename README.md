<img src="docs/logo.png" alt="misskappa" align="right" width="300">

# misskappa

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
