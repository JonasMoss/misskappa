# dev/legacy

**Frozen reference. Do not edit. Not built.**

This directory holds the original `misskappa` implementation as it stood
when we started the C++23 restart. Files are preserved verbatim so the
new tree can read the math out of them while porting:

- `misskappa/` — original R package (C++17 + Armadillo + Rcpp + exceptions).
- `code-analysis/` — original analysis scripts (e.g. simulation runner).
- `notes/` — derivation and design notes that fed into the manuscript.
- `notes/kappa-missing.lyx` — original LyX manuscript (the source of
  truth at the time of the restart). Superseded first by a combined
  `kappa-missing.tex` rewrite (2026-05-26) and then by the three active
  papers under `papers/{ipw,fiml,quadratic}/`. The intermediate combined
  draft was deleted on 2026-06-02 once fully migrated (recoverable from
  git history); the cross-paper plan it produced lives on at
  `papers/split-plan.md`.
- `Justfile` — the old repo-root Justfile.
- `chunks/`, `pdf/`, `llm/`, `skills/`, `todo.md` — supporting material.

If you need to update an old behaviour, port the math into the new tree
under `include/misskappa/` + `src/` rather than editing here.
