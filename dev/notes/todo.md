# misskappa repo TODO

Active repo-level backlog for the public library + R package. Paper-specific
TODOs and the split/section mapping live in the separate private papers repo,
not here. The port plan with the eight-step roadmap is at
`dev/notes/port-plan.md`.

## Roadmap status

- [x] Step 1: skeleton + legacy move into `dev/legacy/`.
- [x] Step 2: port loss layer (`common.cpp` -> `src/loss.cpp` on Eigen + Result).
- [x] Step 3: available-case estimator + variance.
- [x] Step 4: IPW + Gwet estimators.
- [x] Step 5: FIML / EM.
- [x] Step 6: inference consolidation.
- [x] Step 7: R package wiring (Rcpp glue + S3 + sim list).
- [x] Step 8: paper conversion (LyX -> .tex) + scripts wiring — moved to the
  separate private papers repo; out of scope for this public tree.

## Active backlog

- [x] **Complete-data small-sample interval layer.**
      Added exact delete-one summaries as a separate composable C++ inference
      surface for raw categorical kappa with arbitrary agreement matrices,
      count-format Fleiss/Brennan--Prediger kappa, continuous-loss
      Conger/Fleiss kappa, and score-based coefficient alpha. The generic
      `jackknife_covariance()` is coefficient-independent; kappa additionally
      exposes item and chance disagreement for boundary methods. Supported
      complete R fits carry the summaries automatically, and
      `confint(method = "jackknife")` provides natural, Fisher, and
      log-complement transformed t intervals while retaining Wald as the
      default. Exact-boundary kappa uses the current binomial/chance hybrid;
      exact-boundary alpha deliberately returns `NA` pending an alpha-specific
      construction. Missing-data estimators remain unchanged and unsupported
      by this layer.

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

- [x] **R package release (GitHub-first): unified-API refactor shipped at v1.0.0.**
      Decided 2026-06-06, landed and committed check-clean as `v1.0.0`. The
      public surface is one entry point per input shape driven by a shared
      `estimator =` selector `pairwise | ipw | cat_fiml | nt_fiml` on `alpha()`,
      `kappa()`, `kappa_counts()`; `se_type` dropped (sandwich SE only);
      available-case + Gwet pulled from public `kappa()` but reachable for sims
      via internal `estimate_kappa_raw()`; default weight is `nominal`;
      `kappa_continuous()` + `kappa_gwise()` demoted to internal (public
      estimators are just `alpha`/`kappa`/`kappa_counts`); `kappa_quadratic()`
      simplified to `(x, values)` empirical-IF and now returns IFs;
      `kappa_quadratic_fiml()` (robust NT-FIML) backs `estimator = "nt_fiml"`;
      `confint(transform = "fisher")` delta-method CI + `print()` CI columns.
      The test suite, roxygen cross-refs, and `_pkgdown.yml` were resynced to
      the new API and `R CMD check` is clean. Hypothesis tests are no longer
      deferred — see the shipped front-door item below.

- [x] **GitHub-installability: vendor the C++ library into the R package.**
      `dev/vendor-cpp.sh` (= `just vendor`) copies the canonical `src/*.cpp`,
      `src/*.hpp`, and `include/misskappa/*.hpp` into `r-package/src/`, each
      stamped with a `@generated` banner (and cleaned up by that banner on the
      next run, so renames/deletes don't leave stale copies). `Makevars` is now
      self-contained — `PKG_CPPFLAGS = -I. -Imisskappa -DEIGEN_NO_EXCEPTIONS`, C++17, no
      `PKG_LIBS`, no `../../include` — and R compiles every vendored `.cpp`.
      Eigen comes from `LinkingTo: RcppEigen`. Dropped `-fno-exceptions/-rtti`
      (Rcpp needs exceptions; the library never throws, so it's a no-op). The R
      justfile recipes depend on `vendor`, so editing only the canonical sources
      is enough; `just vendor-check` is the CI/pre-commit drift guard. Verified
      by `R CMD build` + `R CMD INSTALL` of the tarball into a throwaway lib with
      `MISSKAPPA_*` env scrubbed: `alpha()`/`kappa()` run, full testthat suite
      201/201. `opt`/`test-opt` still build the static lib for the C++ ctest.

- [x] **Hypothesis-testing front door (shipped).**
      `t.test`-style `kappa_test()` / `alpha_test()` return `htest` over the
      already-built `joint_vcov()` / `wald_test()` engine, exported in v1.0.0.
      Use snake_case, not
      dotted names: `base::kappa` is a generic, so `kappa.test` would
      mis-dispatch as an S3 method. One-sample (`theta = theta0`) and
      independent two-sample (variances add) work for any estimator;
      paired/dependent needs stacked IFs via `joint_vcov()` and works for any
      fit whose `fit$psi` is non-NULL. Same-data method/weight/coefficient
      contrasts stay `wald_test()` territory. Exploratory runners:
      `experiments/studies/24-alpha-equal-cocron/`, `experiments/studies/25-kappa-equal-examples/`.
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
- [ ] **Replace counts-format pairwise weighting with an F&C default.**
      Decision note: `dev/notes/counts-weighting-fc-gls.md`. Add C++ support
      for Fleiss--Cuzick row weights (`r_i - 1`) and a GLS/IVW option
      (`r_i`), use F&C as the public R default, and keep pooled-pairs only as
      a legacy/internal comparison if it remains useful. Validate balanced
      parity against Fleiss/irrCAC and add an unequal-row test showing the
      difference from both current pooled pairs and irrCAC's unit-weighted
      count convention.
- [x] **Fold direct-on-counts CS kappa in + add listwise/pairwise companions.**
      Done 2026-06-15. `r-package/R/kappa_counts_cs.R` fits the exchangeable
      common-mean compound-symmetry (CS) normal kappa **directly on the counts**
      from the per-subject sufficient statistics (m, s, q) — no realizing table,
      no `r_total`, no EM. `kappa_counts(estimator = "nt_fiml")` now uses it
      (reproduces the old table expansion to ≤2.6e-7, point + robust SE; the
      `.counts_exchangeable_scores()` helper is removed). The moment/available-
      case slot stays `fleiss_cuzick`. Prototyped `"pairwise"` and `"listwise"`
      were both **dropped**: pairwise redundant with `fleiss_cuzick` (same moment
      slot, different weighting); listwise is caller-side `na.omit()` row
      dropping, not a baked-in estimator. Shipped counts surface:
      `{fleiss_cuzick, cat_fiml, nt_fiml}`. Tests + regenerated
      `kappa_counts.Rd`; full testthat suite green. Probe 41 /
      `dev/notes/counts-exchangeable-cs-direct.md` are the record.

- [x] **Normal-FIML covariance structure moved into the estimator name; `mean` dropped.**
      Done 2026-06-15. `kappa()`/`alpha()` gained `estimator = "nt_fiml_cs"` (the
      compound-symmetry covariance variant); `"nt_fiml"` stays saturated. The
      `mean` axis is removed from the public surface entirely — alpha is a
      mean-free covariance functional, and for kappa "equal means" *is* the
      counts/exchangeable case (`kappa_counts`), so it was redundant. `covariance`
      and `mean` in `em_options` are deprecated (honoured one release with a
      warning via `.resolve_nt_fiml()`); the internal backends
      (`kappa_quadratic_fiml`, `alpha_continuous`, `kappa_continuous`) keep their
      `covariance`/`mean` params. Rationale: with `mean` gone, `covariance` is a
      lone binary knob on a single estimator, so an orthogonal argument
      over-advertised impossible/redundant cells — the name encodes it instead.
      Tests + regenerated man; suite green.

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

- [x] **Example the CS estimator where it earns its keep (incomplete design).**
      Done 2026-06-17. Added the *synthetic* `dat.tenhove2025` (60×6, connected-
      but-incomplete planned-missing design; `data_raw/tenhove2025.R` regenerates
      it from a fixed seed) and used it to show "saturated non-identified → CS
      fits": worked section in `articles/missingness-estimators.qmd`, an
      `@examples` entry in `kappa.Rd`, and `test-cs-estimators.R` (CS fits +
      IF-identity + saturated errors). Also closed the smaller gaps —
      `kappa_counts.Rd` now examples `estimator = "nt_fiml"`, the getting-started
      vignette names `pairwise_cs` and links the new worked example. Still open
      (minor, parked): the counts `nt_fiml` print header surfaces the deprecated
      `covariance=…, mean=…` labels.

- [x] **Free-mean CS moment foil — SHIPPED as `pairwise_cs` in v1.2.0.**
      Decision flipped from "FIML only": the analytic robust IF (the one
      non-trivial bit flagged below) was built, so `kappa(estimator =
      "pairwise_cs")` ships as the distribution-free MCAR moment companion to
      `nt_fiml_cs` (NEWS 1.2.0, exampled in `kappa.Rd`). History below kept as the
      probe record.
      The estimator story pairs each FIML with a moment foil: saturated/free-mean
      = raw `nt_fiml` vs raw `pairwise` (available-case cov); CS/equal-mean =
      counts `nt_fiml` vs `fleiss_cuzick`; CS/**free-mean** = raw
      `kappa(estimator="nt_fiml_cs")` (Conger) vs the open cell. **Probe 42**
      (`experiments/workbench/42-free-mean-cs-pairwise/`)
      built the foil (rater-center → pooled within-subject pairwise moments on
      residuals → `t3 = sum(mu_j-bar)^2 - tr(AV)`) and confirmed it: nests the
      equal-mean pairwise by mean-model restriction (|free−equal|→5.7e-5),
      matches the free-mean CS FIML to <0.001 (shared limit), the `t3` correction
      removes a spurious small-n Conger bias, MCAR-consistent / MAR-biased (FIML
      robust). Sound but no capability gain over FIML (MCAR-only, less efficient
      at equal means). **Decision still: ship FIML only**, unless the grid is
      wanted in the public API — then wire as `kappa(estimator = "pairwise_cs")`
      (parallel to the `nt_fiml`/`nt_fiml_cs` naming); the one non-trivial bit is
      an analytic robust IF (probe used a bootstrap SE) reusing
      the raw-pairwise IF + mean/`t3` terms. Notes:
      `papers/00_wip/quadratic/dev/notes/cs-fixed-rater-pairwise.md`,
      `dev/refs/icc-stuff/refs.md`.

- [x] **Step 8: paper conversion (LyX -> .tex) + scripts wiring.**
      Moved out of this repo: the manuscripts now live in the separate private
      papers repo, each with its own paper-local todo. Nothing paper-side
      remains in the public tree.

## Deferred

- [x] **Add component-separable vector-valued pairwise agreement.**
      C++ now has internal pairwise/IPW estimators for rectangular
      `subjects x raters x features` data, passed as rater-major flattened
      matrices. The implemented loss class is component-separable with
      diagonal feature weights: Hamming, absolute, squared, and RMS. The R
      wrapper `kappa_vector()` is internal/unexported for experiments. The
      math and FIML boundary are in
      `dev/notes/component-separable-vector-kappa.md`; the CRACKLES pilot uses
      `experiments/workbench/26-crackles-vector-kappa/`.
- [x] **Add full-weight quadratic vector covariance route.**
      The internal R backend `kappa_vector_quadratic()` estimates Conger and
      Fleiss from `(mu, Sigma)` with a full symmetric PSD feature metric `W`,
      using either pairwise-available covariance moments or saturated
      normal-theory FIML. Tensor formulas live in
      `dev/notes/quadratic-vector-tensor-kappa.md`; the CRACKLES experiment
      compares pairwise covariance and NT-FIML under richer component
      missingness.
- [ ] **Categorical FIML for vector profiles.**
      Deferred. Component-wise missingness makes full-profile categorical
      equality and general finite-vector losses latent. A saturated
      multinomial FIML over rater-feature profiles is the natural extension,
      but the support is combinatorial and needs a separate design pass.
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

## Done

- [x] **Prototype coefficient alpha estimators.**
      `estimate_alpha_available()` computes the scored-categorical
      pairwise-covariance alpha baseline with per-subject influence rows, and
      `estimate_alpha_available_continuous()` exposes the same MCAR sandwich
      path for numeric item scores. `estimate_alpha_fiml()` reuses the raw
      categorical EM / Louis machinery to fit the saturated full-response
      distribution and map theta to alpha. The R package exposes numeric alpha
      through `alpha(method = c("available", "fiml"))` and the finite-category
      multinomial EM path through `alpha_cat_fiml()`. Counts input,
      feasibility benchmarking, approximate EM, and experiments remain separate
      follow-ups.
- [x] **Settle the inference-surface policy.**
      `misskappa_estimate` now has a consistent R-facing inference contract:
      all estimators with subject-level influence rows carry them in `fit$psi`
      for `joint_vcov()` and multi-fit `wald_test()`, and
      `stats::influence(fit)` is registered as the public S3 accessor for the
      same matrix.
- [x] **Expose influence functions for FIML estimators.**
      Raw FIML and counts-format FIML now return per-subject influence
      functions from the same reduced observed-score / Louis-information
      path used for their vcov. C++ and R tests assert that
      `crossprod(psi) / n^2` reconstructs the reported covariance, and
      `joint_vcov()` / multi-fit `wald_test()` can include FIML fits.
- [x] **Expose influence functions for counts and continuous estimators.**
      Counts-format available-case and continuous MCAR estimators now return
      per-subject IF matrices through the shared `misskappa_estimate` surface,
      so `joint_vcov()` and `wald_test()` can combine them with other
      same-subject fits.
- [x] **Add Wald tests for misskappa estimates.**
      `wald_test()` now tests single-fit and cross-fit linear hypotheses
      against `misskappa_estimate` objects. Single-fit tests use `vcov()`;
      multi-fit tests use `joint_vcov()` and therefore require aligned
      per-subject influence functions.
- [x] **Keep ASan opt-in.**
      The default `dev` preset stays a portable debug build without
      sanitizers. The `asan` preset remains available for targeted local runs
      once the required libasan / libubsan runtimes are installed; making it
      the default would make ordinary build/test loops too host-dependent.
- [x] **Add `.clang-format` and a `just fmt-cpp` recipe.**
      The formatter config follows the active C++ style with two-space
      indentation and a 100-column limit. The recipe formats active C++
      surfaces only (`include/`, `src/`, `tests/unit/`, `r-package/src/`),
      leaving legacy and vendored third-party code untouched.
- [x] **Wire CI.**
      `.github/workflows/ci.yml` now runs a C++ dev build/test job and an
      R package check job on Ubuntu 24.04. The R job builds the opt static
      library and runs `R CMD check --no-manual` against `r-package/`. The
      irrCAC parity oracle is frozen into the testthat suite, so no external
      oracle package is installed.
- [x] **Add closed rectangular g-wise / Frechet agreement.**
      `estimate_gwise` and `kappa_gwise()` now cover complete
      subjects-by-raters designs with Frechet nominal, Frechet absolute,
      Frechet quadratic, and Hubert disagreement kernels. They return
      Cohen-type and Fleiss-type chance-corrected coefficients with
      influence-function vcov. Categorical chance terms use exact finite
      support enumeration; continuous distances remain direct over item
      tuples. Missing-rater support remains intentionally out of scope.
- [x] **Remove the O(n^2) raw-estimator chance kernels.**
      `estimate_available`, `estimate_ipw`, and `estimate_gwet` now compute
      chance-disagreement V-statistic row sums, column sums, and totals from
      weighted category/rater totals in `src/estimate_raw.cpp`, preserving the
      existing influence-function covariance without looping over all subject
      pairs.
