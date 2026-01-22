---
name: misskappa-maintainer
description: Maintain the misskappa paper + code repo. Use when working in this repository to (1) edit/build the manuscript (kappa-missing.lyx/bib/pdf), (2) develop/check the R package in code/misskappa (Rcpp, tests, R CMD check via Justfile), (3) run paper-facing analysis scripts in code/analysis that emit LaTeX-ready outputs into notes/, and (4) keep llm/ and todo.md in sync with project scope and next actions.
---

# Misskappa Maintainer

Keep work aligned across the manuscript, package implementation, and paper-facing analysis outputs. Prefer the repository’s existing workflows and avoid churn in generated artifacts.

## Quick orientation

- Next actions: `todo.md`
- Scope and goals: `llm/project.md` and `llm/paper-scope.md`
- Repo workflows: `Justfile`
- Package: `code/misskappa/` (see `code/misskappa/README.md`)
- Analysis scripts: `code/analysis/` (typically write LaTeX-ready outputs into `notes/`)
- Notes index: `notes/README.md`

## Workflow decision tree

1. If the task is “package correctness” (API, estimators, inference code, tests), work under `code/misskappa/` and run `just test` and (when needed) `just check`/`just check-cran`.
2. If the task is “paper output” (tables, figures, LaTeX snippets), work under `code/analysis/`, run the relevant `Rscript`, and verify/update the corresponding `notes/*.tex`.
3. If the task is “manuscript text”, edit `kappa-missing.lyx` and keep references in `kappa-missing.bib` consistent. Avoid touching `kappa-missing.pdf` unless explicitly required.
4. If the task is “project memory / planning”, update `todo.md` and `llm/*` (keep them short and high-signal).

## R package workflow (`code/misskappa/`)

Prefer `just` recipes from the repo root:

- Run tests: `just test`
- Full check (tarball + `R CMD check`): `just check`
- CRAN-ish check: `just check-cran`
- Clean build artifacts: `just clean`

Avoid editing generated files unless intentionally regenerating them:

- `code/misskappa/src/RcppExports.cpp`
- `code/misskappa/R/RcppExports.R`
- `code/misskappa/NAMESPACE`

If you change roxygen/Rcpp interfaces, regenerate using the project’s established tooling (don’t hand-edit export stubs).

## Analysis workflow (`code/analysis/`)

Run analysis scripts from the repo root, e.g.:

- `Rscript code/analysis/fleiss1971_countform_irrCAC_report.R`

Expected behavior:

- Scripts use `code/misskappa` and write LaTeX-ready outputs into `notes/` (see `notes/README.md` for what files mean).
- Prefer deterministic outputs: set seeds inside scripts when simulations are involved.

## Manuscript workflow (`kappa-missing.lyx`)

- Prefer small, targeted edits; LyX diffs are noisy.
- Don’t regenerate/overwrite large binaries (especially `kappa-missing.pdf`) unless the task explicitly requires it.

## Bundled references

- Repo map and “do-not-edit” list: `references/repo-map.md`
