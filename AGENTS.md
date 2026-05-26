# AGENTS.md

Working rules for coding agents in the misskappa repo.

## What misskappa is

A C++23 library for **agreement coefficients** — Cohen / Fleiss / Brennan-Prediger
and friends — with arbitrary numbers of raters, arbitrary pairwise loss matrices,
and support for missing categorical ratings under MCAR and MAR. The primary
audience is methods developers. The R package under `r-package/` is a thin
Rcpp wrapper over the prebuilt static library, intended for empirical examples
and the accompanying manuscript under `paper/`.

The first paper (`paper/kappa-missing.tex`, targeting Psychometrika) studies
the MCAR estimators (available-case and IPW) and shows that Gwet's earlier
inferential method is inconsistent.

## Non-negotiables

- **C++23**, built with `-fno-exceptions -fno-rtti`. Eigen is included under
  `EIGEN_NO_EXCEPTIONS`. Failures are values: `std::expected<T, Error>` aliased
  as `Result<T>` in `include/misskappa/result.hpp`. The C++23 floor is carried
  by exactly that one feature — `std::expected` is the error model — a
  deliberate single-feature dependency (GCC 13 / Clang 17).
- **No virtual functions, no `concept`/`requires`, no inheritance hierarchies
  on the hot path.** Extension is via free function templates over duck-typed
  callables (e.g., the loss function). Plain aggregate structs for data
  (`Estimation`, `EmResult`, option bags like `EmOptions`).
- **No `std::function` on the hot path.** Continuous loss functions pass as
  small POD wrappers around function pointers.
- **`irrCAC` is the oracle for closed-data agreement.** C++ golden tests read
  JSON fixtures written by `tests/tools/regen_oracle.R` (which calls installed
  irrCAC) and assert agreement to a documented tolerance. CI does not invoke R.
- **Legacy is reference, not foundation.** The original C++17 + Armadillo
  implementation lives under `dev/legacy/misskappa/` for math reference. Do not
  build it, do not edit it. The new tree is rewritten on Eigen + `Result<T>`.

## Where things live

- `include/misskappa/` — public headers (stable surface).
- `src/` — implementations plus private `detail_*.hpp`.
- `tests/unit/` — focused unit tests (doctest, header-only).
- `tests/golden/` — fixture-based parity checks against irrCAC.
- `tests/fixtures/` — checked-in JSON. Regenerate via `tests/tools/regen_oracle.R`.
- `tests/tools/` — maintainer-only fixture-generation scripts (R).
- `r-package/` — Rcpp bindings; consumes the prebuilt `libmisskappa.a`,
  separate from and not part of the C++ build.
- `paper/` — manuscript, bibliography, figures, tables, simulation scripts,
  curated results. Each subdir documented in `paper/AGENTS.md`; the prose /
  table / figure style guide is `paper/STYLE.md`.
- `dev/legacy/` — frozen reference: original R package, C++ implementation,
  analysis scripts, LyX manuscript, supporting notes. Unbuilt, do not edit.
- `dev/notes/` — repo-level development notes (port plan, validation plan, todo).
- `external/` (ignored) — optional upstream source mirrors for reading.
- `resources/` (ignored) — local data and scratch.

## Build

```sh
cmake --preset dev
cmake --build --preset dev
ctest --preset dev

cmake --preset opt
cmake --build --preset opt
ctest --preset opt
```

`dev` is the local debug build (AddressSanitizer + UBSan). `opt` is the release
build (`-O3 -DNDEBUG -march=native`) and is the artifact the R package links.

Repo-root `justfile` wraps the common loops: `just build`, `just test`
(`dev` build + ctest), `just opt`, `just test-opt`, `just r-install`,
`just r-check` (reinstall + R-level tests), `just paper` (delegates to
`paper/justfile`), and `just regen-oracle` (regenerates `tests/fixtures/`).
`just` with no recipe lists them.

The R bindings link the prebuilt non-sanitized `opt` `libmisskappa.a`;
`r-package/src/Makevars` makes the package objects depend on it, so a C++
header change correctly forces the R glue to recompile against the new ABI.

## R package direction

The R package is a methods-developer interface over the C++ library, not a
second implementation. Prefer thin exported R wrappers around a single C++
entry point per estimator, with C++ argument structure kept visible. Small R
helpers are fine when they compose existing wrappers, validate R-shaped inputs,
or preserve names for inspection; they should not contain parallel kappa logic.

Public R surface:

- `kappa(x, method, weight, ...)` — single user entry point. `method` is one of
  `"available"`, `"ipw"`, `"fiml"`, `"gwet"`. `weight` is a string naming the
  loss or a user-supplied matrix.
- S3 generics on `misskappa_estimate`: `print`, `coef`, `vcov`, `as.data.frame`,
  `confint` (Wald CIs from `vcov`).
- `sim` — a list of simulation closures (`sim$mcar`, `sim$mar`, `sim$jsm`).
  Exposed as a single object so simulation helpers occupy one manual entry
  rather than many.
- Datasets carried over from the legacy package: `dat.fleiss1971`,
  `dat.gwet2014`, `dat.klein2018`, `dat.zapf2016`.

## Namespace layout

The C++ surface lives in one namespace, `misskappa`. Sub-namespaces are
introduced only when there is a clear boundary:

- `misskappa::loss` — weight matrix and loss-function factories.
- `misskappa::em` — EM machinery for FIML.
- `misskappa` (top-level) — the four estimator entry points and shared types.

No deep nesting. No per-estimator namespaces.

## Conventions

- Lowercase `snake_case` for filenames and free functions; `CamelCase` for
  types; constants are `snake_case` (`version_major`, `na_code`, etc.).
- Public headers include with `#include "misskappa/foo.hpp"`. Private detail
  headers under `src/detail_*.hpp` include with relative paths.
- Comments only when the why is non-obvious. Headers carry the contract.
- Keep `dev/notes/todo.md` as the single active backlog. Fold finished planning
  docs into it or remove them; do not create parallel roadmaps.
- Keep `paper/AGENTS.md` and `paper/dev/todo.md` current when the manuscript or
  simulation scope changes.
- Commit every finished user request as a coherent completed change before
  handing back, unless the user explicitly asks not to commit. Keep work on the
  current branch unless explicitly asked to branch.

## Working with the legacy reference

`dev/legacy/misskappa/` contains the original Armadillo + Rcpp implementation.
When porting an estimator:

1. Read the math out of the legacy source (`kappanp.cpp`, `emdiscrete.{h,cpp}`,
   `kappaml.cpp`, `common.cpp`).
2. Rewrite on Eigen + `Result<T>` in the new tree.
3. Add a parity test under `r-package/tests/testthat/test-parity-legacy.R` that
   loads both packages and asserts numerical agreement.
4. Add the irrCAC fixture comparison under `tests/golden/`.

Do not copy code mechanically — the legacy code uses exceptions, `Rcpp::stop`,
and `arma::` types throughout. The port is a rewrite that preserves the
algorithm.
