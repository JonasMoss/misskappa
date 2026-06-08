# misskappa public-API redesign + equality tests — handoff (2026-06-08)

A continuation brief for a fresh session. The redesign and the equality-test
verbs are **done and green**; this note records the final surface, the decisions
behind it (so they are not relitigated), and the open items.

## Done (this work)

- **Unified estimator selector.** All estimators take one `estimator=` argument:
  - `kappa(x, estimator, weight, values, ...)` — `"ipw"`/`"cat_fiml"` for
    categorical ratings (any `weight`); `"pairwise"`/`"nt_fiml"` for the
    quadratically weighted scored coefficient (treat `x` as numeric scores,
    require `weight = "quadratic"`). Absorbs the former `kappa_continuous`.
  - `kappa_counts(x, estimator, weight, ...)` — `"pairwise"`/`"cat_fiml"`.
  - `alpha(x, estimator, values, ...)` — `"pairwise"`/`"cat_fiml"`/`"nt_fiml"`
    (replaces `method` + `type`).
- **Equality-test verbs.** `kappa_test()` / `alpha_test()` return an `htest`;
  cover one-sample, independent two-sample, paired (same-subject), and G-way
  homogeneity. Built only on `coef`/`vcov`/`fit$psi`.
- **Influence functions exposed as `fit$psi`** (documented component;
  `vcov == crossprod(psi)/n^2`). The `influence()` S3 method was **removed**
  (clashes with `stats::influence` diagnostics; clutters the reference index).
- **Internal-only (kept for sims, reachable via `misskappa:::`):** `joint_vcov()`
  + `wald_test()` (the arbitrary-contrast engine), `kappa_quadratic()`,
  `kappa_continuous()`, `kappa_gwise()`, `estimate_kappa_raw()` (the latter
  exposes the dropped Gwet / available-case-categorical methods).
- **Vocabulary:** `available`→`pairwise`, `identity`→`nominal`; `cat_fiml`
  (saturated multinomial) vs `nt_fiml` (robust normal-theory) named distinctly.
- **`se_type="normal"` removed** package-wide; the sandwich/IF covariance is the
  only SE path, so every estimator exposes `fit$psi`.
- **Worked examples (public-API demos for the papers):**
  `experiments/studies/24-alpha-equal-cocron` (paired alpha on `cocron::knowledge`) and
  `experiments/studies/25-kappa-equal-examples` (anxiety homogeneity, Westlund–Kurland MS
  two-sample, FEES `multiagree` benchmark, mcduff dependent-on-missing-data).
  Both use the public verbs only. Headline numbers: alpha 0.764/0.804 χ²≈4.36;
  MS quadratic κ 0.525/0.626 z≈−1.02; FEES exact match to `multiagree`;
  mcduff 0.643/0.646 χ²≈0.016.
- **Docs:** `NEWS.md`, `_pkgdown.yml`, `pkgdown/index.md`, `AGENTS.md` updated;
  roxygen/NAMESPACE regenerated; full `testthat` suite passes; in-place
  `R CMD INSTALL` works.

Public surface: `alpha`, `kappa`, `kappa_counts`, `kappa_test`, `alpha_test`,
`sim` + datasets + S3 (`print`/`coef`/`vcov`/`confint`/`as.data.frame`).

## Decisions (settled — do not relitigate)

- `estimator=` over `method`+`type`; `cat_fiml`/`nt_fiml`; `nominal`/`pairwise`.
- `kappa_quadratic` merged into `kappa` (no separate quadratic front door).
- Dropped from the public surface (kept internal for sims): Gwet, available-case
  for *categorical* `kappa`, `kappa_gwise`, the full `kappa_continuous` dispatch.
- `se_type="normal"` not exposed (sandwich strictly dominates).
- IFs are the `fit$psi` component, **not** an `influence()` function
  ("exposed, not named"). The `rownames`-on-`$psi` alignment guard was
  **declined** (subject ids on rating matrices are rare; equal-`n` + documented
  same-row-order is enough).
- Two thin verbs (`kappa_test`/`alpha_test`) over one internal worker — kept for
  tab-completion discoverability and the per-family `coef` default.
- The Tier-2 contrast engine (`joint_vcov`/`wald_test`) stays internal.

## Open items (next)

1. **C++ vendoring (release blocker).** `src/Makevars` uses
   `MISSKAPPA_INCLUDE ?= $(CURDIR)/../../include`; the headers live at
   `repo/include/misskappa/*.hpp`, *outside* the package, so an in-place
   `R CMD INSTALL` works but `R CMD build`/`R CMD check` (a tarball) cannot
   compile. Fix: vendor `include/misskappa/` into the package (e.g.
   `r-package/inst/include/` or `src/`) so headers travel with the tarball.
   Tracked on the release checklist in `dev/notes/todo.md`.
2. **Experiments framing/cleanup.** Done: `experiments/` is now organized by
   lifecycle, with canonical paths in `experiments/INDEX.md` and root-level
   compatibility README stubs. Every experiment **except 24/25** still uses the
   pre-redesign API, so treat frozen runners as records rather than migration
   targets. Findings to carry in:
   - duplicate `12-` prefix is resolved at the top level by the lifecycle
     buckets (`studies/12-clean-mar-dgp` versus
     `archive/pre-redesign/12-quadratic-rare-disagreement`);
   - `18`/`23` (`*-paper-simulation`) are *fuller* than the papers' own
     `scripts/*_simulation.R` — distill, don't blind-delete;
   - `22` (`quadratic-nt-fiml-validation`): its `se_type` sandwich-vs-normal
     framing is dead post-redesign, but it is the cited NT-FIML-vs-`magmaan`
     validation record;
   - `26-crackles-vector-kappa` is the user's WIP (leave it);
   - several experiments are cited in memories (`project_experiment_findings` →
     01–04; `project_equality_test_examples` → 09/14/22) — don't dangle those.
3. **Manuscript fold-in.** Add short "Illustration" subsections to the
   alpha-missing and quadratic papers using experiments 24/25, plus bib entries
   (Feldt 1969, Hakstian–Whalen 1976, Donner et al. 1996, Vanbelle 2017,
   Diedenhofen & Musch `cocron`).

## Pointers

- Plan file: `~/.claude/plans/i-m-eager-to-do-merry-codd.md`.
- Auto-memories: `project_equality_test_examples`, `project_r_package_release`,
  `feedback_api_minimal_footprint`.
- Key files: `r-package/R/{kappa.R, agreement_test.R, joint_vcov.R}`,
  `r-package/{NEWS.md,_pkgdown.yml,pkgdown/index.md}`, `AGENTS.md`,
  `experiments/studies/{24-alpha-equal-cocron,25-kappa-equal-examples}`.
