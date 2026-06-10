# 30 — Geometry of the raw cat_fiml MLE face and the identification guard

**Question.** The reduced-gradient identification guard
(`assert_reduced_gradient_identified`, src/estimate_fiml.cpp:448) hard-fails
raw cat_fiml in study-29 planned-missingness designs. Is the non-identified
nuisance benign for kappa? Is there a canonical point estimate when the MLE
is a face, or do we need centroid-style machinery?

**Findings (2026-06-10).**

1. `01-margins-pin-kappa.R` — With every rater pair co-observed, all pairwise
   joint margins of theta are linear functionals lying in the row space of the
   pattern-margin map A, so Conger/Fleiss/BP kappa are *exactly constant* on
   any set {theta : A theta = const}. Verified on the classic 2×2×2
   three-way-interaction direction: theta moves, all three kappas frozen to
   1e-10. Population identification of kappa under the design guard is exact,
   not asymptotic.

2. `02-reproduce-guard-overfiring.R` — The guard fires erratically in the
   sparse regime (n=40–100, C=4–5, R=4–5, pair-rotation designs) on data that
   pass the design guard, i.e. where kappa is population-identified. This is
   the study-29 blocker.

3. `03-face-geometry-and-smoothing.R` — For a firing case (n=40, C=4, R=4):
   - The MLE set is a convex polytope (observed-data log-lik is concave in
     theta since pattern-cell probabilities are linear). Its dimension is
     ≥ 189 of 255 — enormous; explicit centroid computation is hopeless and
     unnecessary.
   - Hand-rolled saturated EM from 6 random starts: identical log-likelihood,
     but kappa varies by 4.7e-3 across the face. So finite-sample width is
     nonzero (empirical pattern margins are mutually incompatible, the MLE
     compromises, zero-count cells are not all pinned) — but it is ~20×
     smaller than the SE (0.11) and vanishes asymptotically.
   - Dirichlet(1+delta) MAP smoothing (add delta pseudo-count per cell in the
     M-step) selects the approximate analytic center of the face: at
     delta=1e-4 the start-to-start kappa spread collapses 4.7e-3 → 1e-7 while
     moving the estimate by only 2e-4. delta=1e-3 over-smooths (shift 4e-3,
     comparable to face width).

4. `04-face-width-vs-n.R` — Multi-start kappa spread (8 starts, 5 seeds) at
   n = 40, 160, 640 for the C=4, R=4 pair design. Median width 2.7e-3 at
   n=40, ~1e-14 at n=160, ~3e-12 at n=640: the width is a sparse-table
   small-n phenomenon driven by zero-count pattern cells, and it vanishes
   abruptly (much faster than n^(-1/2)) once those cells fill in. Caveat: one
   n=40 seed had width 9.8e-2 — the same order as the SE — so occasionally
   the face is genuinely wide and a point estimate is not honestly unique.
   That is the case a *diagnostic* (not a hard error) should flag.

**Conclusion.** The design guard (complete pair co-observation) is the correct
and complete hard identification gate for kappa. The gradient-projection guard
conflates benign finite-sample face width (and numerical noise at an
approximate MLE) with fatal non-identification and should be demoted from a
hard error to a diagnostic. The canonical point selection, if wanted, is the
analytic center via delta-smoothed EM — one line in the M-step — not a
polytope centroid (#P-hard, parameterization-dependent).
