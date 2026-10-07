# misskappa repo TODO

Active repo-level backlog for the public library + R package. Paper-specific
TODOs live in the separate private papers repo, not here. Finished work is
recorded in `r-package/NEWS.md` and git history; keep this file to what is
still open.

## 2.0.0 release

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
- [ ] **Repo pruning (needs the author's go-ahead for cross-repo moves).**
      Move `experiments/` and the paper-specific `dev/notes/*` into the
      private papers repo (fix the scripts' package paths and the private
      papers' references); drop `papers/guessing` and `papers/kappa-sd` from
      the public tree after making their own repos public or archiving them
      on Zenodo, then fix the README badge/Contents, `just paper`, and the
      AGENTS.md papers bullets; drop `docs/logo.png` (README now uses
      `r-package/man/figures/logo.png`).
- [x] **Release check.** Full `R CMD check` with the Quarto vignette:
      Status OK (0 errors, 0 warnings, 0 notes).
- [ ] **Tag `v2.0.0`** once the repo pruning above is settled (the author
      pushes).

## Open (library)

- [ ] **MASS-RENAME the two methods-paper citations on preprint + publication.**
      WIP scaffold LANDED 2026-06-17: the two papers are now forward-referenced
      under working titles in `inst/CITATION` (keys `moss-kappa-missing`,
      `moss-quadratic`), in the DESCRIPTION `Description` field, and in the
      roxygen `@references` of `kappa()` / `alpha()` / `kappa_counts()` — all
      marked "(in preparation)". What remains is to **update those same spots
      TWICE** as the papers progress:
        * **On arXiv:** add the arXiv id + URL and switch the wording from "in
          preparation" to "preprint".
        * **On publication:** swap in the journal name / volume / DOI and drop the
          preprint wording.
      Touch points to grep when renaming: `inst/CITATION`, `DESCRIPTION`, and
      `@references` blocks in `R/kappa.R` (search `in preparation`). Lock the
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

- [ ] **Test projected Satterthwaite calibration for pivotal robust scores.**
      Decision note: `dev/notes/projected-satterthwaite-robust-scores.md`.
      Generalize the Pearson numerator--variance projection to the full
      sandwich-studentized rank-\(q\) score and match its projected meat
      fluctuation to an approximate Hotelling reference. First audit the
      multivariate formula, raw-versus-centered score meat, and fitted-nuisance
      contribution. Then create one focused sibling-magmaan probe: a correctly
      specified five-indicator
      one-factor model, normal control plus a Vale--Maurelli/Fleishman
      skewness-2/excess-kurtosis-7 condition, \(n=30,50,100,200\), and true
      scalar plus global score nulls. Compare ordinary robust
      \(\chi_q^2\), projected working-normal Hotelling, unprojected diagnostic,
      and an empirical/delete-one diagnostic; keep MV/FMG procedures secondary
      because they correct a different nonpivotal weighted-\(\chi^2\) limit.
      Report type-I error, effective df, Wishart-shape diagnostics, meat
      conditioning/singularity, failures, and runtime before adding power or
      broader nonnormal families. Separately verify the conjectured
      \((n-k)/2\) partial-correlation rule and the preliminary
      \(n/[3+2/\{P(P-1)\}]\) tau-equivalent-alpha rule only after the SEM null
      probe fixes the general calculation.
- [ ] **Complete the small-sample interval coverage matrix.**
      Long-run target: every meaningful estimator x input shape x weighting
      combination exposed by the public R API (plus the retained simulation
      helpers) should either carry exact delete-one summaries or document why
      that interval is mathematically inapplicable. Keep the C++ layer
      coefficient-agnostic and composable; the R layer should only select and
      label the appropriate summaries.

      - Add complete-data refits for the structured scalar routes:
        `pairwise_cs`, `nt_fiml`, and `nt_fiml_cs` kappa; `nt_fiml_cs` alpha;
        and count-format quadratic `nt_fiml`. Verify each delete-one point
        functional against its fitted estimator instead of silently attaching
        an empirical substitute.
      - Cover the standalone/internal surfaces: `kappa_quadratic()`,
        `estimate_kappa_fiml_multi()`, available-case/Gwet continuous kappa,
        and all supported `g > 2` nominal/Hubert/linear functionals.
      - Generalize delete-one summaries to vector-valued ratings: the
        component-separable pairwise/IPW routes and the full-`W` quadratic
        `pairwise`/`nt_fiml`/`nt_fiml_cs` routes, preserving feature metrics,
        fixed score ranges, and coefficient ordering.
      - Finish boundary handling: derive and validate an alpha-specific
        perfect-agreement interval, and define the weighted/discrete boundary
        construction for count matrices with unequal rater totals.
      - Add a table-driven R test matrix over estimator, input shape, weight,
        coefficient family, transform, and boundary status so future API
        additions cannot create silent holes.
      - In the later missing-data pass, implement genuine delete-one refits for
        IPW, categorical FIML, normal FIML, and count FIML while preserving the
        missingness pattern; specify how all-missing rows and the effective
        subject count map to degrees of freedom and influence rows. Keep this
        phase separate from the current complete-data paper.
- [ ] **Reliability family (omega): closed-form core + fit-based tier; rename to `missrel`.**
      Add McDonald's omega and the multidimensional family as misskappa coefficients
      (omega is alpha's congeneric sibling; keep all internal-consistency + inter-rater
      reliability here, SEM proper stays in magmaan). Method + paper live in
      `magmaan/papers/closed-form-omega`: the closed forms are already implemented in
      pure R there (`r-package/R/harness-h.R`: `closed_form_general_loadings`,
      `closed_form_h_general`, ...) and in the magmaan C++ core
      (`measures::frontier::reliability`); the fit-free path is ~100 lines of
      covariance algebra and needs NO SEM machinery. **RENAME to `missrel` while
      pre-CRAN** (continues the `miss*` brand; "reliability" covers both inter-rater
      [kappa] and internal-consistency [alpha/omega]).

      Three tiers by dependency:
      - **Core (no SEM dep):** closed-form `omega_total` / `omega_h` / `H` /
        `H_general` as covariance functionals from the (missing-data) covariance, with
        ADF / missing-data analytic SE. REUSES the existing FIML/EM covariance +
        sandwich/influence-function scaffold (same machinery as `alpha`; omega is just
        a different `g(Sigma)`). Headline novelty = the omega family **under missing
        data**, fit-free, with an analytic interval — genuinely unserved (gap below).
      - **Fit-based (`Suggests: lavaan`):** categorical/polychoric omega, free
        (non-proportional) bifactor `omega_h`, correlated-error omega, general-structure
        (Raykov-Shrout) reliability. Fit in lavaan, delta SE over `vcov()` (robust
        under MLR/WLSMV). Value-add over `semTools::compRelSEM` = the SE/CI
        (compRelSEM is POINT-ONLY for every coefficient) + missing data. Point
        estimates are not novel (RS2002 / Green-Yang / semTools); the interval +
        missing-data handling are.
      - **Speed backend (magmaan, optional, NOT declared):** magmaan is not an
        installable R package yet, so it CANNOT go in `Suggests` today. Wire via a
        runtime flag `backend = c("lavaan","magmaan")`, default `"lavaan"`; reach
        magmaan only on explicit opt-in through a SOFT reference
        (`asNamespace("magmaan")`, guarded by `requireNamespace`), never `magmaan::`
        — keeps `R CMD check` clean while undeclared and is forward-compatible:
        promote to `Suggests` once the Saga packaging makes magmaan
        `install_github`-able. DEFAULT lavaan for reproducibility (do NOT
        auto-pick-magmaan-if-present, or published numbers become backend-dependent).
        magmaan is purely a sim/throughput accelerator (~50x).

      Dividing line (which tier a coefficient lands in): simple structure + diagonal
      residuals -> closed form, no fit (core); cross-loadings / correlated errors /
      polychoric / free bifactor -> fit (lavaan tier).

      SOFTWARE GAP this fills (verified 2026-06-29): `semTools::compRelSEM` (the
      confirmatory reliability tool) returns POINT ESTIMATES ONLY for every coefficient
      (no SE/CI for anything); analytic intervals are packaged only for single-factor
      omega (`MBESS::ci.reliability`); multidimensional `omega_total` / bifactor
      `omega_h` have only bootstrap (psych) or Bayesian (Pfadt 2022 / `Bayesrel`). No
      frequentist analytic distribution-free CI for multidim/bifactor omega exists in
      software, and none under missing data. NOTE package scope ⊃ paper scope: the
      paper is the fit-free closed form; the lavaan tier is package completeness, not a
      paper claim.
- [ ] **Measure alpha FIML feasibility before expanding scope.**
      Exact EM over the full fixed-category pattern space is plausible for
      small Likert batteries but scales as `C^R`. Add guardrails/documentation
      around feasible `C, R`, benchmark memory/time on representative `C x R`
      grids, and only then consider approximate or sparse EM variants
      (active-support EM, low-rank/composite-likelihood starts, or other
      approximations) as a separate follow-up.
- [ ] **Ten Hove / variance-component SEM route for quadratic agreement.**
      Decision 2026-06-14: start with the variance-component SEM target, not
      the looser unconstrained exchangeable covariance projection. The
      variance-component model is the stronger assumption: it imposes the
      latent random-effect interpretation and nonnegative/PSD component
      constraints, while exchangeability only constrains the marginal
      covariance pattern. Treat `magmaan` as an external/off-limits validation
      oracle for this repo: the misskappa work should not edit or depend on
      that sibling package. Use `experiments/probes/39-ten-hove-icc/` as the
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
- [ ] **Add an alpha experiment suite.**
      Simulate simple complete-data, MCAR, and MAR categorical item batteries;
      compare complete-data alpha, pairwise/available covariance alpha, and
      categorical FIML alpha for bias, SE calibration, and coverage. Add an
      applied example on `psych::bfi`, reporting one alpha per Big Five
      personality scale under the available-case and FIML paths, with clear
      notes about item scoring/reversal and the fixed-category support used.
      Starter categorical smoke/calibration runners now live in
      `experiments/probes/15-alpha-categorical-smoke/` and
      `experiments/probes/16-alpha-calibration-sweep/`. Experiment 16 has a
      capped `B = 200` calibration pass: `5^6, n = 4000` was too slow
      for an ordinary run, and `5^6` FIML retains downward bias through
      `n = 1000`. Independent saturated-EM checks agree with the package
      estimate, and the quadrature alpha truth agrees with large complete-data
      Monte Carlo, so treat the residual high-dimensional categorical-FIML
      bias as a finite-sample fact rather than an implementation blocker.
      Manuscript notes can state this matter-of-factly and leave debiasing,
      penalized/smoothed EM, bootstrap/jackknife corrections, or approximate
      sparse variants as optional follow-up. The normal-FIML comparison,
      applied `psych::bfi` example, alpha-specific Louis spectrum diagnostic,
      and explicit high-dimensional FIML guardrails remain open. The proposed
      paper-facing simulation grid is in
      `dev/notes/alpha-missing-simulation-study.md`: use an essential-tau
      setup derived from Zhang-Yuan's six-item calibration, include
      Zhang-Yuan's congeneric cell, add stronger congeneric and two-factor
      non-congeneric cells, borrow Enders/Savalei missingness rates, and use
      clean observed-anchor MAR rather than truth-dependent missingness.
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
