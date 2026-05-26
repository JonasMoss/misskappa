# misskappa repo TODO

Active repo-level backlog. Paper-specific TODO lives at `paper/dev/todo.md`.
The port plan with the eight-step roadmap is at `dev/notes/port-plan.md`.

## Currently in flight

- [x] Step 1: skeleton + legacy move into `dev/legacy/`.
- [ ] Step 2: port loss layer (`common.cpp` -> `src/loss.cpp` on Eigen + Result).
- [ ] Step 3: available-case estimator + variance.
- [ ] Step 4: IPW + Gwet estimators.
- [ ] Step 5: FIML / EM.
- [ ] Step 6: inference consolidation.
- [ ] Step 7: R package wiring (Rcpp glue + S3 + sim list).
- [ ] Step 8: paper conversion (LyX -> .tex) + scripts wiring.

## Backlog

- [ ] Wire CI (compile-only job under multiple compilers; R CMD check job).
- [ ] Decide whether to enable the `asan` preset by default once Fedora's
      libasan / libubsan runtime packages are installed locally.
- [ ] Add `.clang-format` and a `just fmt-cpp` recipe.
- [ ] Re-evaluate whether counts-format input and continuous ratings should
      come back into the library after the paper is submitted.
