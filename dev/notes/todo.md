# misskappa repo TODO

Active repo-level backlog for the public library + R package. Paper-specific
TODOs live in the separate private papers repo, not here. Finished work is
recorded in `r-package/NEWS.md` and git history; keep this file to what is
still open.

## 3.0.1 release

- [x] **cat_fiml warning removed.** The `null_frac` warning (fired on 70-93%
      of fits at package defaults) is gone; `se_cutoff_ratio` is a new silent
      diagnostic. A warning on it was tried and dropped: in the companion
      paper's main simulation it flagged hard designs, not bad fits within a
      design.
- [x] **Default `max_iter` 100000** (was 10000) for the EM estimators. The
      estimate settles early but the SE needs the default `tol = 1e-8`: on
      24 slow fits the SE was up to 2.5x too large at 1e-6 and 2x at 1e-7,
      and within 4% at 1e-8 (estimates within 3e-4 throughout).
- [x] **`nt_fiml` singular covariance.** Complete data use the sample-moment
      influence rows (equal to the old sandwich to 5e-15; works with
      perfectly agreeing raters). Incomplete data with a singular fitted
      covariance give a `singular_covariance` error that says so.
- [ ] **Release.** CRAN rejected 3.0.0 (HTML issues; being fixed
      separately), so 3.0.0 never shipped. For the resubmission, decide the
      version number and fold the 3.0.1 NEWS entries into the first-release
      notes.
- [ ] **`nt_fiml` on incomplete data with a singular fitted covariance** still
      fails (now with a clear error). In the companion paper's simulation this
      is common on balanced incomplete block designs with high agreement
      (up to 90% of fits at n = 100). The pseudo-inverse rescue (tag
      `archive/nt-fiml-soo-rescue`) was rejected in June: rescued fits were
      biased with SE/SD about 0.35.
- [ ] **`nt_fiml` near-singular edge cases (from the companion paper's
      rerun).** (a) "Numerical error." remains on 0-11% of small-n fits
      (n = 40, rater-specific MCAR worst): the observed-information solve in
      the sandwich fails, same near-singular geometry as the singular-
      covariance error; give it an informative message too. (b) About 2% of
      the fits that do return on a balanced incomplete block design at n = 40
      have absurd SEs (> 0.5, up to ~100) from a nearly singular fitted
      covariance; consider flagging or refusing when rcond(Sigma) is tiny.
- [ ] **Optional estimate-side check.** If a selection warning is ever wanted
      back, bound the move directly: step from theta-hat along the null-space
      component of the gradient until a support cell hits zero (ratio test)
      and warn only if the coefficient changes by more than 0.1 SE.

## 3.0.0 release

- [x] **Agreement only.** Coefficient alpha removed (C++ and R); last present
      at tag `alpha-final`, which the alpha-missing paper pins. A future alpha
      or reliability package starts from that tag.
- [x] **Rename.** `kappa()` / `kappa_counts()` / `kappa_test()` ->
      `agreement()` / `agreement_counts()` / `agreement_test()` (no more
      masking `base::kappa`).
- [x] **Rating scale.** `categories =` declares the full scale; character,
      logical, and factor labels; numeric labels never truncated.
- [x] **Subjects.** Empty subjects dropped everywhere; `fit$psi` rows named;
      paired tests align by subject.
- [x] **Repo pruning.** `experiments/` and research notes moved to the private
      papers repo; `docs/logo.png`, C++ SQUAREM, and `wald_test()` dropped.
- [ ] **Private papers on the new names.** The WIP paper scripts call
      `kappa()` / `kappa_counts()` / `kappa_test()` and internals such as
      `estimate_kappa_fiml_multi()` (now `research_kappa_fiml_multi()`; the
      full map is in NEWS.md); switch them to the new names (or pin them to
      `alpha-final`).
- [ ] **Completed papers out of the public tree.** Drop `papers/guessing` and
      `papers/kappa-sd` after making their own repos public or archiving them
      on Zenodo (needs the author), then fix the README Contents and badge,
      `just paper`, and the AGENTS.md papers bullets.
- [x] **Release check.** Full `R CMD check` with the Quarto vignettes:
      Status OK (2026-11-01).
- [ ] **Tag `v3.0.0`** and push it with `alpha-final` (the author pushes).
- [ ] **CRAN prep (if submitting).** Windows/macOS CI (r-lib
      check-r-package), win-builder/rhub, and the arXiv ids for the two
      in-preparation citations.

## 2.0.0 release (never tagged; superseded by 3.0.0)

- [x] **Fix R CMD check.** The magmaan-oracle tests now use a pure-R reference
      EM and sandwich (test helper); `utils` in Imports, `coefficientalpha` in
      Suggests. (ab8f000)
- [x] **Breaking R API cleanup.** Brennan-Prediger only for a declared scale;
      counts `nt_fiml` -> `nt_fiml_cs`; strict `...` and `em_options`; 1.2.0
      shims and `fd_h` removed; no g > 2 `"pairwise"`; `pairwise_cs`
      internal; `sim` removed; `kappa_test()` defaults (`paired = FALSE`,
      first coefficient); fits carry `n`; `as.data.frame()` with CIs; one
      `.new_estimate()` constructor; internal help pages dropped. (8ad0187)
- [x] **C++ engine pruning.** `EmOptions::fd_h` removed; the Louis diagnostic
      and `build_psi_from_phi()` moved to private detail headers; the
      estimate header states the C++ API is the R package's engine.
- [x] **Docs pass.** Vignette, articles (validation rewritten), README,
      AGENTS.md.
- [x] **Release check.** Full `R CMD check` with the Quarto vignette:
      Status OK (0 errors, 0 warnings, 0 notes).
## Open (library)

- [ ] **MASS-RENAME the two methods-paper citations on preprint + publication.**
      WIP scaffold LANDED 2026-06-17: the two papers are now forward-referenced
      under working titles in `inst/CITATION` (keys `moss-kappa-missing`,
      `moss-quadratic`), in the DESCRIPTION `Description` field, and in the
      roxygen `@references` of `agreement()` / `agreement_counts()` — all
      marked "(in preparation)". What remains is to **update those same spots
      TWICE** as the papers progress:
        * **On arXiv:** add the arXiv id + URL and switch the wording from "in
          preparation" to "preprint".
        * **On publication:** swap in the journal name / volume / DOI and drop the
          preprint wording.
      Touch points to grep when renaming: `inst/CITATION`, `DESCRIPTION`, and
      `@references` blocks in `R/agreement.R` (search `in preparation`). Lock the
      final titles first; the quadratic paper is still untitled (coefficient-led,
      Biometrics-target), so its working title "Quadratically weighted agreement
      coefficients" is provisional.
- [ ] **POD-pointer overloads — only when a non-Eigen consumer shows up.**
      The Eigen API is already thin for the two consumers we have (in-tree
      C++ tests and R via Rcpp+RcppEigen). If we ever add a Python / Julia /
      standalone-CLI binding, add overloads in
      `include/misskappa/estimate.hpp` that take
      `(const int* ratings, int n, int R, const double* W, int C)` and write
      results into caller buffers, with the Eigen versions becoming
      one-liners over `Eigen::Map`. Until then, not worth the surface
      doubling.
- [ ] **Re-evaluate counts-format and continuous-rating scope after submission.**
      Counts-format and continuous-rating estimators are currently present in
      the library and R package. Revisit their public-surface priority after
      the paper submission path is settled.

## Post-2.0 / research (not release blockers)

- [ ] **Complete the small-sample interval coverage matrix.**
      Long-run target: every meaningful estimator x input shape x weighting
      combination exposed by the public R API (plus the retained simulation
      helpers) should either carry exact delete-one summaries or document why
      that interval is mathematically inapplicable. Keep the C++ layer
      coefficient-agnostic and composable; the R layer should only select and
      label the appropriate summaries.

      - Add complete-data refits for the structured scalar routes:
        `pairwise_cs`, `nt_fiml`, and `nt_fiml_cs` kappa;
        and count-format quadratic `nt_fiml`. Verify each delete-one point
        functional against its fitted estimator instead of silently attaching
        an empirical substitute.
      - Cover the research surfaces: `research_kappa_quadratic()`,
        `research_kappa_fiml_multi()`, available-case/Gwet continuous kappa,
        and all supported `g > 2` nominal/Hubert/linear functionals.
      - Generalize delete-one summaries to vector-valued ratings: the
        component-separable pairwise/IPW routes and the full-`W` quadratic
        `pairwise`/`nt_fiml`/`nt_fiml_cs` routes, preserving feature metrics,
        fixed score ranges, and coefficient ordering.
      - Finish boundary handling: define the weighted/discrete boundary
        construction for count matrices with unequal rater totals.
      - Add a table-driven R test matrix over estimator, input shape, weight,
        coefficient family, transform, and boundary status so future API
        additions cannot create silent holes.
      - In the later missing-data pass, implement genuine delete-one refits for
        IPW, categorical FIML, normal FIML, and count FIML while preserving the
        missingness pattern; specify how all-missing rows and the effective
        subject count map to degrees of freedom and influence rows. Keep this
        phase separate from the current complete-data paper.
- [ ] **Ten Hove / variance-component SEM route for quadratic agreement.**
      Decision 2026-06-14: start with the variance-component SEM target, not
      the looser unconstrained exchangeable covariance projection. The
      variance-component model is the stronger assumption: it imposes the
      latent random-effect interpretation and nonnegative/PSD component
      constraints, while exchangeability only constrains the marginal
      covariance pattern. Treat `magmaan` as an external/off-limits validation
      oracle for this repo: the misskappa work should not edit or depend on
      that sibling package. Use `experiments/probes/39-ten-hove-icc/` (now in the private papers repo) as the
      validation seed: lmer/SEM point-estimate parity, robust-vs-naive SE
      checks, estimand-matched comparison to saturated `nt_fiml`, and the MAR
      rater-camps failure mode. If a misskappa shortcut is added, name and
      document whether it is the true variance-component estimator or only the
      weaker compound-symmetry covariance estimator. Repeated-measure/vector
      profiles are the main follow-up requirement: the natural target is a
      fixed-rater variance-component SEM with PSD subject-profile and
      residual-profile components, not an unconstrained block average. Alpha
      can reuse the same covariance fit after the quadratic-kappa path is
      validated. Status 2026-06-14: branch `ten-hove` has the scalar
      variance-component covariance option for quadratic `kappa()` and
      continuous `alpha()`; vector-valued ratings explicitly reject that option
      until the block-structured repeated-measure model is implemented and
      validated.
- [ ] **Categorical FIML for vector profiles.**
      Deferred. Component-wise missingness makes full-profile categorical
      equality and general finite-vector losses latent. A saturated
      multinomial FIML over rater-feature profiles is the natural extension,
      but the support is combinatorial and needs a separate design pass.

## Done (details in NEWS.md and git history)

- Step 1: skeleton + legacy move into `dev/legacy/`.
- Step 2: port loss layer (`common.cpp` -> `src/loss.cpp` on Eigen + Result).
- Step 3: available-case estimator + variance.
- Step 4: IPW + Gwet estimators.
- Step 5: FIML / EM.
- Step 6: inference consolidation.
- Step 7: R package wiring (Rcpp glue + S3 + sim list).
- **Complete-data small-sample interval layer.**
- **R package release (GitHub-first): unified-API refactor shipped at v1.0.0.**
- **GitHub-installability: vendor the C++ library into the R package.**
- **Hypothesis-testing front door (shipped).**
- **Fold direct-on-counts CS kappa in + add listwise/pairwise companions.**
- **Normal-FIML covariance structure moved into the estimator name; `mean` dropped.**
- **Example the CS estimator where it earns its keep (incomplete design).**
- **Free-mean CS moment foil — SHIPPED as `pairwise_cs` in v1.2.0.**
- **Add component-separable vector-valued pairwise agreement.**
- **Add full-weight quadratic vector covariance route.**
- **Prototype coefficient alpha estimators.**
- **Settle the inference-surface policy.**
- **Expose influence functions for FIML estimators.**
- **Expose influence functions for counts and continuous estimators.**
- **Add Wald tests for misskappa estimates.**
- **Keep ASan opt-in.**
- **Add `.clang-format` and a `just fmt-cpp` recipe.**
- **Wire CI.**
- **Add closed rectangular g-wise / Frechet agreement.**
- **Remove the O(n^2) raw-estimator chance kernels.**
- **Fleiss-Cuzick count weighting as the counts moment default** (shipped in 1.1.0).
