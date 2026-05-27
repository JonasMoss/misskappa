# dev/legacy

**Frozen reference. Do not edit. Not built.**

This directory holds the original `misskappa` implementation as it stood
when we started the C++23 restart. Files are preserved verbatim so the
new tree can read the math out of them while porting:

- `misskappa/` — original R package (C++17 + Armadillo + Rcpp + exceptions).
- `code-analysis/` — original analysis scripts (e.g. simulation runner).
- `notes/` — derivation and design notes that fed into the manuscript.
- `kappa-missing.lyx` — original LyX manuscript (the source of truth at
  the time of the restart). Replaced by `papers/combined/kappa-missing.tex`
  (which is itself being split into `papers/{ipw,fiml,quadratic}/`).
- `kappa-missing.pdf` — last compiled PDF snapshot.
- `Justfile` — the old repo-root Justfile.
- `chunks/`, `pdf/`, `llm/`, `skills/`, `todo.md` — supporting material.

If you need to update an old behaviour, port the math into the new tree
under `include/misskappa/` + `src/` rather than editing here.
