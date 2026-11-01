# misskappa repo TODO

Active repo-level backlog for the public library + R package. Paper-specific
TODOs live in the separate private papers repo, not here. Finished work is
recorded in `r-package/NEWS.md` and git history; keep this file to what is
still open.

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

- [ ] **Rework the cat_fiml non-identification warning (next release).**
      `.warn_null_frac()` warns when more than 1% of the delta-method gradient
      lies in the dropped null space of the Louis information. Evidence from
      the companion paper's simulation probe (private papers repo, C = 5,
      R = 5-6, n = 1000) and study 32 (June):
        * It fires too often to mean anything: about half of all fits, 97-100%
          under MAR, 67-93% under MCAR with rater-specific rates.
        * The estimate it warns about does not move: refits with
          `start_alpha` 0.01 / 0.1 / 1 differ by <= 3e-5 and `flatten = 0.1`
          by <= 2e-4 (sampling SD 0.012-0.025); study 32 never saw > 0.01 SE.
        * The real fragility is the SE, and `null_frac` does not detect it.
          Between `info_rcond` 1e-4 and 1e-6 the SE changes by a factor of
          1.0-3.0 (weakly identified directions just above the cutoff). The
          cell with the largest `null_frac` (0.065) moved at most 1.33x; the
          one with the smallest (0.012) moved up to 3x.
        * The message recommends `flatten`, which inflates SEs 1.6-3.4x.
      Proposal:
        1. Stop warning on `null_frac` alone; keep `fit$null_frac` in the
           object and document it. If an estimate-side check is still wanted,
           bound the movement directly: step from theta-hat along the
           null-space component of the gradient until a support cell hits
           zero (ratio test) and warn only if the resulting change in the
           coefficient exceeds 0.1 SE. (Sketch; check the reduced-coordinate
           bookkeeping.)
        2. Add an SE-stability diagnostic from the eigendecomposition the
           variance step already does (no refit): recompute the SE with the
           cutoff at `info_rcond / 100`, store the ratio as
           `fit$se_cutoff_ratio`, and warn when it exceeds a threshold (1.25
           to start; calibrate on the companion paper's simulation). Message:
           the SE depends on the information cutoff (x vs y), typical of a
           sparse saturated table, so treat Wald intervals with caution.
        3. Drop the `flatten` recommendation from the message (flatten is a
           uniqueness device; its SEs are conservative).
      Touch points: `R/agreement.R` (`.warn_null_frac`), the variance step in
      `src/estimate_fiml.cpp` (`apply_fiml_variance_workspace`), the shared
      warning call in `R/research.R`, docs for `em_options`, NEWS, tests.

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
