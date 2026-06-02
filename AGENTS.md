# AGENTS.md

Working rules for coding agents in the misskappa repo.

## What misskappa is

A C++23 library for **agreement coefficients** — Cohen / Fleiss / Brennan-Prediger
and friends — with arbitrary numbers of raters, arbitrary pairwise loss matrices,
and support for missing categorical ratings under MCAR and MAR. The primary
audience is methods developers. The R package under `r-package/` is a thin
Rcpp wrapper over the prebuilt static library, intended for empirical examples
and the accompanying manuscripts under `papers/`.

The combined draft was split into three papers (see
`papers/split-plan.md`) and then deleted on 2026-06-02 once fully
migrated (recoverable from git history):
`papers/ipw/` (Psychometrika; IPW + AC/Gwet comparison),
`papers/fiml/` (Psychometrika or Biometrika; FIML / semiparametric efficiency
under MAR), and `papers/quadratic/` (Biometrics; closed-form quadratic =
Lin's CCC on pairwise-available data).

The three spinoffs are intended as a simultaneous submission family, not
strictly sequential papers. Cross-references among Papers A/B/C are allowed and
expected; use stable provisional labels ("Paper A", "Paper B", "Paper C" or
"companion paper") until final titles, names, and bibliography keys are chosen.

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
  implementation was deleted on 2026-06-02; it lives on in git history (last
  at `6b30f98`, retrievable with `git checkout 6b30f98 -- dev/legacy/misskappa`)
  and is the source of the oracle values the unit tests are frozen against. The
  new tree is rewritten on Eigen + `Result<T>`.

## Where things live

- `include/misskappa/` — public headers (stable surface).
- `src/` — implementations plus private `detail_*.hpp`.
- `tests/unit/` — focused unit tests (doctest, header-only).
- `tests/golden/` — fixture-based parity checks against irrCAC.
- `tests/fixtures/` — checked-in JSON. Regenerate via `tests/tools/regen_oracle.R`.
- `tests/tools/` — maintainer-only fixture-generation scripts (R).
- `r-package/` — Rcpp bindings; consumes the prebuilt `libmisskappa.a`,
  separate from and not part of the C++ build.
- `papers/` — research manuscripts. `papers/ipw/`, `papers/fiml/`,
  `papers/quadratic/` are the three spinoffs (plus the independent
  `papers/ac1-paper/`). Each carries its own `AGENTS.md`, `STYLE.md`,
  `justfile`, bibliography, figures, tables, scripts, and results. The
  cross-paper plan and section mapping live in `papers/split-plan.md`;
  cross-paper todos in `papers/todo.md`. The combined draft they were
  split from was deleted on 2026-06-02 (recoverable from git history).
- `docs/` — repo-level documentation assets, including shared artwork such as
  the project logo. Generated docs output is under `docs/site/` and is ignored;
  do not commit generated HTML.
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
`just r-check` (reinstall + R-level tests), `just paper <slug>` (delegates
to `papers/<slug>/justfile`, e.g. `just paper ipw pdf`), and
`just regen-oracle` (regenerates `tests/fixtures/`). Documentation recipes are
`just docs-r` (pkgdown), `just docs-cpp` (Doxygen), `just docs` (combined
local site), and `just docs-clean`. `just` with no recipe lists them.

The R bindings link the prebuilt non-sanitized `opt` `libmisskappa.a`;
`r-package/src/Makevars` makes the package objects depend on it, so a C++
header change correctly forces the R glue to recompile against the new ABI.

## Documentation

The documentation site is generated, not hand-edited. `pkgdown` is the front
door for the R package under `r-package/`, Quarto articles under
`r-package/vignettes/articles/` carry the mathematical exposition, and Doxygen
generates the C++ API reference from `include/misskappa/` into `/cpp/` inside
the built site.

- Put long-form math, estimator definitions, missingness assumptions, and
  validation narratives in Quarto/pkgdown articles.
- Keep Doxygen/header comments focused on C++ API contracts: inputs, output
  order, dimensions, error conditions, and short formulas only.
- Do not commit generated documentation output under `docs/site/` or
  `docs/doxygen/`. R `.Rd` files under `r-package/man/` are also generated
  during docs builds and ignored for now. Build locally with `just docs`;
  GitHub Actions publishes generated artifacts.
- GitHub Pages deployment is gated by the repository variable
  `ENABLE_PAGES_DEPLOY=true`. Private-repository Pages requires a GitHub plan
  that supports it; otherwise the workflow still builds docs without deploying.
- When docs tooling changes, update `justfile`, `.github/workflows/docs.yml`,
  `r-package/_pkgdown.yml`, and this section together as needed.

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
- `misskappa` (top-level) — all estimator entry points and shared types.
  Input shape is encoded in the function name (`_counts`, `_continuous`)
  rather than via per-shape sub-namespaces; this keeps callers from
  having to switch namespaces when comparing estimators across shapes.

No deep nesting. No per-estimator or per-shape namespaces. Private
implementation helpers live in `misskappa::detail` (in `src/detail_*.hpp`)
or in anonymous namespaces inside the `.cpp` that owns them.

## Conventions

- Lowercase `snake_case` for filenames and free functions; `CamelCase` for
  types; constants are `snake_case` (`version_major`, `na_code`, etc.).
- Public headers include with `#include "misskappa/foo.hpp"`. Private detail
  headers under `src/detail_*.hpp` include with relative paths.
- Comments only when the why is non-obvious. Headers carry the contract.
- Agent instructions live in `AGENTS.md` files and installed Codex skills; do
  not keep copied agent-skill snapshots in the tree.
- Keep `dev/notes/todo.md` as the single active backlog. Fold finished planning
  docs into it or remove them; do not create parallel roadmaps.
- Keep each `papers/<slug>/AGENTS.md` and `papers/<slug>/dev/todo.md` current
  when the corresponding manuscript or simulation scope changes.
- Commit every finished user request as a coherent completed change before
  handing back, unless the user explicitly asks not to commit. Keep work on the
  current branch unless explicitly asked to branch.

## Git hygiene

Multiple agents may be working in this repo at the same time. To avoid
clobbering each other's pending work:

- **Stage explicitly.** Run `git add <specific paths>` for the files your task
  touched. Do not use `git add -A`, `git add .`, or `git commit -a` — those
  bulk forms silently pull another agent's in-flight changes into your commit,
  producing a commit whose message no longer matches its contents.
- **Sanity-check before committing.** Run `git diff --cached --name-only` and
  confirm every file listed belongs to the task you just finished. If you see
  files you didn't touch, `git restore --staged <path>` them out.
- **Do not `git reset HEAD~N` to undo a commit you regret.** Mixed reset (the
  default) leaves the working tree alone, so HEAD and the working tree diverge
  — every change in the discarded commit reappears as unstaged, often after
  the next agent has already moved on from that state. Prefer:
    - `git revert <sha>` to back out a committed change cleanly, or
    - `git reset --soft HEAD~1` followed by `git restore --staged <paths>` for
      the files that shouldn't have been in the commit, then re-commit. Only
      do this if the commit is the tip of `HEAD` and is yours alone.
- Avoid amending or rebasing a commit unless you authored it and it has not
  been built on by another agent's work.

## Working with the legacy reference

The original Armadillo + Rcpp `misskappa` (the source of the unit-test oracle
values) was deleted from the tree on 2026-06-02 but remains in git history. To
consult it when porting or debugging an estimator, retrieve it with
`git checkout 6b30f98 -- dev/legacy/misskappa` (then `git rm --cached` /
`rm -rf` when done) and read the math out of `kappanp.cpp`, `emdiscrete.{h,cpp}`,
`kappaml.cpp`, `common.cpp`. Rewrite on Eigen + `Result<T>` in the new tree and
back the result with an irrCAC fixture under `tests/golden/`.

Do not copy code mechanically — the legacy code uses exceptions, `Rcpp::stop`,
and `arma::` types throughout. The port is a rewrite that preserves the
algorithm.
