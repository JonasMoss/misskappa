# AGENTS.md

Working rules for coding agents in the misskappa repo.

## What misskappa is

A C++17 library for **agreement coefficients** — Cohen / Fleiss / Brennan-Prediger
and friends — with arbitrary numbers of raters, arbitrary pairwise loss matrices,
and support for missing categorical ratings under MCAR and MAR. The primary
audience is methods developers. The R package under `r-package/` is a thin
Rcpp wrapper over the library (vendored into the package for installs),
intended for empirical examples and the accompanying manuscripts. The C++
library is the R package's engine: its interface follows the package's needs
and carries no separate stability guarantee (stated in `estimate.hpp`).

The combined draft was split into a family of companion papers — IPW + AC/Gwet
comparison, FIML / semiparametric efficiency under MAR, and the closed-form
quadratic = Lin's CCC route — plus standalone papers on the AC1 guessing model
and on coefficient alpha under missing data. These manuscripts live in a
separate private repository and are out of scope for this public tree; use
stable provisional labels ("Paper A", "Paper B", "Paper C", or "companion
paper") in notes that reference them.

## Non-negotiables

- **C++17.** The library target builds with `-fno-exceptions -fno-rtti`
  (applied to the library only, so the R package's Rcpp glue still compiles with
  the exceptions R expects); Eigen is included under `EIGEN_NO_EXCEPTIONS`.
  Failures are values: a hand-rolled `Result<T>` in
  `include/misskappa/result.hpp` stands in for `std::expected<T, Error>` (which
  is C++23), keeping the floor at C++17 so the package builds on the widest
  toolchain range.
- **No virtual functions, no `concept`/`requires`, no inheritance hierarchies
  on the hot path.** Extension is via free function templates over duck-typed
  callables (e.g., the loss function). Plain aggregate structs for data
  (`Estimation`, `EmResult`, option bags like `EmOptions`).
- **No `std::function` on the hot path.** Continuous loss functions pass as
  small POD wrappers around function pointers.
- **`irrCAC` is the oracle for closed-data agreement.** Its complete-data
  estimates and variances are frozen as literals in
  `r-package/tests/testthat/test-parity-oracle.R`, so CI never needs irrCAC.
  Refresh them from a local irrCAC install when a convention changes.
- **Legacy is reference, not foundation.** The original C++17 + Armadillo
  implementation was deleted on 2026-06-02; it lives on in git history (last
  at `6b30f98`, retrievable with `git checkout 6b30f98 -- dev/legacy/misskappa`)
  and is the source of the oracle values the unit tests are frozen against. The
  new tree is rewritten on Eigen + `Result<T>`.

## Where things live

- `include/misskappa/` — public headers (stable surface).
- `src/` — implementations plus private `detail_*.hpp`.
- `tests/unit/` — focused unit tests (doctest, header-only).
- `r-package/tests/testthat/` — R tests, including the frozen irrCAC parity
  oracle and a pure-R normal-FIML reference (`helper-normal-fiml-oracle.R`).
- `r-package/` — Rcpp bindings; consumes the prebuilt `libmisskappa.a`,
  separate from and not part of the C++ build.
- `papers/` — the completed, published manuscripts that ship with this repo
  (`papers/guessing/` and `papers/kappa-sd/`), each with its own sources,
  bibliography, figures, tables, and reproduction scripts. Work-in-progress
  manuscripts (the kappa companion papers and the standalone coefficient-alpha
  paper) live in a separate private repository, not here.
- `docs/` — repo-level documentation assets, including shared artwork such as
  the project logo. Generated docs output is under `docs/site/` and is ignored;
  do not commit generated HTML.
- `dev/notes/` — repo-level development notes; `todo.md` is the single backlog.
- `external/` (ignored) — optional upstream source mirrors for reading.
- `resources/` (ignored) — local data and scratch.

## Build

```sh
cmake --preset dev
cmake --build --preset dev --parallel "$(nproc)"
ctest --preset dev

cmake --preset opt
cmake --build --preset opt --parallel "$(nproc)"
ctest --preset opt
```

`dev` is the local debug build (`-O0 -g`, no sanitizers — the separate `asan`
preset adds AddressSanitizer + UBSan). `opt` is the release build
(`-O3 -DNDEBUG`, no `-march`, matching R's own flags); it produces
`libmisskappa.a`, which the C++ unit tests link and which `just r-dev` links the
R glue against (see "R package build" below).

**Always build in parallel.** Pass `--parallel "$(nproc)"` (the `just` recipes
do this for you); without it the Unix Makefiles generator compiles the ~26
`-O3` Eigen-heavy translation units *serially*, which on a multi-core box is the
difference between ~3 min and ~20 min for a clean `opt` build. Compiles also
route through **ccache** when it is installed: `CMakeLists.txt` sets
`CMAKE_CXX_COMPILER_LAUNCHER` via a guarded `find_program` (a no-op if ccache is
absent), so a cache-warm rebuild after `just clean` or a branch switch is
near-instant (only edited files miss). This is the CMake build's own ccache
wiring, independent of the R-side `~/.R/Makevars` setup described below.

Repo-root `justfile` wraps the common loops: `just build`, `just test`
(`dev` build + ctest), `just opt`, `just test-opt`, `just r-install`,
`just r-dev` (fast dev install — see "R package build" below),
`just r-check` (fast `R CMD check` — no vignettes; `just r-check-release` for
the full vignette-building release check), `just paper <slug>` (delegates
to `papers/<slug>/justfile`). Coverage recipes are
`just cpp-cov` (clang source-based instrumentation + `llvm-cov` over the
doctest suite; per-file summary plus HTML under `build-cov/coverage-html/`),
`just r-cov` (`covr` over the testthat suite, which also reports the vendored
C++ as hit through Rcpp), and `just cov` for both. Documentation recipes are
`just docs-r` (pkgdown), `just docs-cpp` (Doxygen), `just docs` (combined
local site), and `just docs-clean`. `just` with no recipe lists them.

### R package build: vendored (release) vs linked (dev)

The R package builds in two modes; pick by what you're doing.

- **`just r-install` / `just r-check` / `just r-check-release` — the
  release/self-contained form.**
  `just vendor` copies the canonical C++ (`src/` + `include/`) into
  `r-package/src/` so the package compiles every source itself, with no reach-up
  to `../..` and no prebuilt library. This is the form that ships and that
  `remotes::install_github` / CRAN build, so `r-package/src/Makevars` must stay
  self-contained — never make it reference `../../build-opt/libmisskappa.a` or
  any out-of-tree path (it would fail on CRAN's machines, which only unpack
  `r-package/`). `just vendor-check` guards against the `@generated` vendored
  copies drifting; never hand-edit them — edit the canonical sources and
  re-`just vendor`.
- **`just r-dev` — the fast dev loop.** Compiles only the Rcpp glue
  (`RcppExports.cpp` + `rcpp_glue.cpp`) and links it against the prebuilt `opt`
  `libmisskappa.a`, so an install does a lib-only `cmake` build plus a two-file
  compile instead of recompiling the whole library. It works from a throwaway
  `build-rdev/` mirror with the dev-only `dev/r-makevars-dev` swapped in for
  `src/Makevars` (R only honours `OBJECTS=` set in the *package* Makevars, so the
  override has to be a package-level file in a copy), which keeps the committed
  self-contained Makevars untouched and keeps anything dev-only out of the
  tarball. It deliberately skips `r-docs`: roxygenise → `pkgbuild::compile_dll`
  would recompile every source, the exact cost we're avoiding. Run
  `just r-install` after changing the R-facing interface (new `[[Rcpp::export]]`,
  `@export`, NAMESPACE, man pages) — `r-dev` only relinks existing exports.

**Fast vs. full `R CMD check`.** Use `just r-check` in the dev loop: it skips the
Quarto vignettes and downgrades absent Suggests to a NOTE (mirroring CI), so it is
much cheaper than the old behaviour. Run `just r-check-release` **only before a
release/tag** — that is the one that renders the vignettes and checks the full
tarball as `remotes::install_github` / CRAN and the Docs workflow do. Vignette
breakage is caught there and in the Docs workflow, not in the inner loop. Skipping
vignettes leaves two expected WARNINGs (`checking files in 'vignettes'` and
`checking package vignettes`, both because `inst/doc` is absent) — benign and
identical to CI; treat anything beyond those two as a real finding.

For faster iteration, route R's compiler through ccache in `~/.R/Makevars`
(`CXX17 = ccache g++`, and the other `CXX*`/`CC`): repeated R compiles — the
vendored `r-install` path included — then return cached objects instead of
rebuilding the ~40M-per-file Eigen objects. This is a personal/machine setup,
not a repo file, and is invisible to CRAN. The fast `just r-check` additionally
sets `CCACHE_NOHASHDIR`/`CCACHE_BASEDIR` so even `R CMD check`'s throwaway temp
build dir hits the cache (in-place `r-install` always did; the temp dir otherwise
busts ccache's path-sensitive hash).

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
  `docs/doxygen/`. R `.Rd` files under `r-package/man/` are generated by
  roxygen but tracked as package documentation; regenerate them with
  `just r-docs` or `just docs-r` and commit intentional changes. Track pkgdown
  source assets under `r-package/pkgdown/` (for example favicons), but not the
  built site. Build locally with `just docs`; GitHub Actions publishes
  generated artifacts.
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

Public R surface (3.0.0, GitHub-first). The package is agreement-only:
coefficient alpha was removed in 3.0.0 (last present at tag `alpha-final`,
which the alpha-missing paper pins).

- Estimators select the method with a single `estimator=` argument:
  - `agreement(x, estimator, weight, g, values, em_options, ...)` — raw ratings.
    `"ipw"` / `"cat_fiml"` for categorical ratings (integer codes, any
    `weight`); `"pairwise"` / `"nt_fiml"` / `"nt_fiml_cs"` for the
    quadratically weighted scored coefficient (numeric scores, weight defaults
    to `"quadratic"`). `g > 2` selects the g-wise kernels (`"ipw"`,
    `"cat_fiml"`); a 3-D array selects vector-valued ratings, with `W`,
    `ranges`, `feature_weights`, `loss` passed through `...`.
  - `agreement_counts(x, estimator, weight, values, r_total, em_options)` —
    counts-format input; `"fleiss_cuzick"`, `"cat_fiml"`, or `"nt_fiml_cs"`
    (exchangeable compound symmetry, quadratic only).
  The two FIMLs are named distinctly (`cat_fiml` = saturated multinomial,
  `nt_fiml` = robust normal-theory) so one token never denotes two estimators;
  a `_cs` suffix always means the compound-symmetry covariance; the unweighted
  loss is `"nominal"`; the MCAR moment estimator is `"pairwise"`.
- `ratings_to_counts()` converts rater-identified ratings to counts.
- `agreement_test()` — `t.test`-style Wald equality tests
  returning `htest` (independent samples by default, `paired = TRUE` uses the
  stacked influence rows).
- S3 methods on `misskappa_estimate`: `print`, `coef`, `vcov`,
  `as.data.frame` (with CI columns), `confint` (Wald by default,
  `method = "jackknife"` on supported complete-data fits, `transform =
  "fisher"` / `"log_complement"`), and `stats::influence` (the per-subject
  influence rows stored as `fit$psi`). Fits also carry `n`.
- Datasets: `dat.fleiss1971`, `dat.gwet2014`, `dat.klein2018`,
  `dat.zapf2016`, `dat.mcduff2019`, `dat.vanbelle2019`,
  and the synthetic `dat.tenhove2025`.

Interface rules:

- Every argument either applies to the call or errors. Inapplicable `...`
  arguments and unknown `em_options` entries are errors, never silently
  ignored.
- Brennan-Prediger needs a declared rating scale. Categorical estimators take
  it from the categories; the scored estimators report it only when `values`
  (or `ranges` for vector ratings) declares it.
- One constructor, `.new_estimate()` in `R/estimate_object.R`, builds every
  fit object.

Kept in the tree but **internal** (`@noRd`, reachable with `:::` for
simulations): Gwet's reweighting and available-case for categorical `kappa`
(`estimate_kappa_raw()`), the full continuous dispatch (`kappa_continuous()`),
g-wise agreement (`kappa_gwise()`), the compound-symmetry moment estimator
(`.kappa_pairwise_cs()`), the multi-weight categorical FIML
(`estimate_kappa_fiml_multi()`), the closed-form `kappa_quadratic()`, and the
IF-based joint covariance / Wald contrast engine (`joint_vcov()`,
`wald_test()`).

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
  types; constants are `snake_case` (e.g. `na_code`).
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
back the result with a frozen irrCAC parity value in
`r-package/tests/testthat/test-parity-oracle.R`.

Do not copy code mechanically — the legacy code uses exceptions, `Rcpp::stop`,
and `arma::` types throughout. The port is a rewrite that preserves the
algorithm.
