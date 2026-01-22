# misskappa (paper + R package)

This repo contains:

- A manuscript source: `kappa-missing.lyx` (compiled PDF: `kappa-missing.pdf`) and bibliography `kappa-missing.bib`.
- An R package with Rcpp: `code/misskappa/`.
- Paper-facing analysis scripts: `code/analysis/` (typically write LaTeX-ready outputs into `notes/`).
- Working notes/derivations: `notes/`.
- LLM-facing “project memory”: `llm/` (high-level scope + goals).

## Quick orientation

- Project scope/goals: `llm/project.md`, `todo.md`.
- Notes index: `notes/README.md`.
- Package docs: `code/misskappa/README.md`.

## Common workflows (preferred)

Use `just` from the repo root (see `Justfile`):

- `just test` — run package tests via `testthat::test_local()`.
- `just check` / `just check-cran` — run `R CMD check` on the built tarball.
- `just clean` — remove `*.Rcheck`, tarballs, and compiled artifacts under `code/misskappa/src/`.

Run analysis scripts from the repo root, e.g. `Rscript code/analysis/<script>.R`.

## Conventions / hygiene

- Avoid editing generated files unless you are intentionally regenerating them:
  - `code/misskappa/src/RcppExports.cpp`
  - `code/misskappa/R/RcppExports.R`
  - `code/misskappa/NAMESPACE` (roxygen-generated)
- Prefer small, surgical diffs in `kappa-missing.lyx` (LyX diffs are noisy).
- Don’t modify large binary artifacts (e.g. `kappa-missing.pdf`) unless the task explicitly requires it.
- `notes/_autosave/` is not meant for edits/commits; keep generated note outputs in `notes/` where scripts expect them.

