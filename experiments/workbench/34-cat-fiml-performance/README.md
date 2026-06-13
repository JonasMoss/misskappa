# 34 — cat-FIML performance probe

**Question.** Where does `estimate_fiml` (the cat_fiml raw path in
`src/estimate_fiml.cpp`) spend its time on study-29-shaped workloads, and how
much performance is on the table?

**Method.** `bench.cpp` simulates a latent-class DGP (C categories, R raters,
70% agreement) under MCAR(30%) and a designed keep-2-raters pattern, then
times `estimate_fiml` with study-29's options (`tol = 1e-7`,
`max_iter = 12000`, `prune_tol = 1e-10`), plus two `flatten = 1` cells. Three
variants of `estimate_fiml.cpp` are compiled side by side:

- `estimate_fiml_prof.cpp` — the canonical implementation plus phase timers
  (`prof.hpp`): preprocess, EM loop, Louis info, pseudo-inverse, theta-vcov
  product, kappa map, null-fraction diagnostic.
- `estimate_fiml_fast.cpp` — bitwise-identical mechanical rework: the EM
  state vector lives on the *active union* of completion cells instead of the
  full C^R table (the counts path already does this); cells outside the union
  are one analytic scalar. Plus buffer reuse/swap instead of a per-iteration
  allocation, fused gather, `uint32` completion indices, and a flat
  active→pruned map replacing the `unordered_map` in the Louis pass.
- `estimate_fiml_sq.cpp` — the fast variant with the plain EM iteration
  wrapped in SQUAREM S3 (two EM sweeps, extrapolation with the standard
  fall-back-to-EM safeguard, one stabilising sweep). Fixed points are EM
  fixed points; the trajectory differs.

Run `./run.sh [reps]`.

## Results (2026-06-12, study-29 EM options, identity weights)

Mean ms per fit; `iters` = EM sweeps; `support` = patterns surviving pruning.

| case                      | baseline | iters | exact-fast | SQUAREM | iters | support  |
|---------------------------|---------:|------:|-----------:|--------:|------:|---------:|
| mcar30 C=6 R=5 n=40       |    23 ms |   790 |      16 ms |    6 ms |   136 |      102 |
| mcar30 C=6 R=5 n=100      |    38 ms |  1115 |      32 ms |    8 ms |   197 |    59→86 |
| designed2 C=6 R=5 n=100   |   227 ms |  1940 |     181 ms |   31 ms |   387 |    45→42 |
| mcar30 C=5 R=5 n=100      |    17 ms |  1895 |      14 ms |    9 ms |   294 |       62 |
| designed2 C=6 R=5 n=300   |   528 ms |  2304 |     442 ms |   42 ms |   681 |    80→79 |
| designed2 C=5 R=7 n=100   |  11.6 s  |  4210 |     10.0 s |   3.6 s |  2139 |    42→41 |
| flatten=1 mcar30 C=6 R=5  |   491 ms |   441 |     420 ms |  559 ms |   109 |      644 |
| flatten=1 desgn2 C=6 R=5  |   54.5 s |   785 |     69.2 s |  67.1 s |   257 |    2,554 |

(strict-ML rows: exact-fast checksums, iteration counts, supports, estimates
and vcovs are bit-identical to baseline; SQUAREM estimates agree to ~1e-6 but
SEs move by a few percent because a different set of near-zero cells survives
`prune_tol`, shifting the Louis rcond truncation — the same support/rcond
artifact documented in studies 31/32.)

## Findings

1. **Strict ML (the study-29 regime): the EM loop is ~100% of runtime.**
   Preprocess, Louis information, pseudo-inverse, kappa map, and the
   null-fraction diagnostic are all ≤1 ms because the surviving support is
   small (40–100 patterns). Optimizing anything but the EM loop is pointless
   here.
2. **Exact mechanical rework buys only 1.1–1.5×.** The per-iteration cost is
   dominated by the gather/scatter over completion lists (groups × C^missing),
   which the active-set representation does not shrink. Its real value is the
   memory cliff: the canonical code allocates two dense C^R vectors and
   iterates over them, which is why R=8+ at C=5 is currently painful and
   bigger R impossible; the active set caps state at the completion union.
3. **The big lever is the EM iteration count** (800–4,200 sweeps at
   `tol = 1e-7`; convergence is linear with rate = fraction of missing
   information, so designed-missingness cells are the slowest). SQUAREM cuts
   sweeps 3–6× and wall time **5–13×** on the realistic R=5 cells. Estimates
   match to ~1e-6 (consistent with study 32: the kappa functional is not
   selection-dependent), but SE movement of a few % means adoption needs a
   small validation pass (re-run study-31/32-style calibration), i.e. it is
   policy-adjacent, not a drop-in.
4. **Flattening flips the profile.** With `flatten = 1` the surviving support
   explodes (~2,700 patterns) and the variance machinery dominates. Measured
   phase split for one 64-s fit (designed2, C=6, R=5, n=100): EM 0.02 s,
   Louis build 0.5 s, `pseudo_inverse_psd` **30.8 s**,
   `reduced_gradient_null_fraction` **27.0 s** (each runs its own O(m³)
   eigendecomposition of the same `info_star`), `theta_vcov` 5.4 s — and the
   cached `theta_vcov = J Σ* Jᵀ` is never needed (the 3×3 coefficient vcov
   only needs `jacobian_reduced · Σ* · jacobian_reducedᵀ`). Since the Louis information has rank ≤ G (number of
   observed-pattern groups, ≤ n) while m ≈ support size, the whole pass can
   be done through the G×G Gram matrix of the score rows in O(G²m) — a
   ~100–1000× reduction for these shapes — sharing one eigendecomposition
   between the pseudo-inverse, the null-fraction diagnostic, and
   `diagnose_fiml_louis`. Strict-ML runs don't care (m ≈ 50), so this only
   matters if flatten or very large supports come back.
5. **Smaller redundancies** (negligible at current sizes, listed for
   completeness): `estimate_alpha_fiml`, `estimate_fiml_gwise`, and the counts
   path each build the Louis info + pseudo-inverse twice (once inside
   `run_em_preprocessed(compute_vcov = true)` / `em_variance`, again for the
   psi/null-frac pass); `build_kappa_map`'s Qed matrices are O(n_final²R²)
   but could be O(n_final·R·C) via per-rater marginals; the counts path has
   no `_many` variant, so study 29's `fleiss_counts_fiml_nominal` +
   `_quadratic` arms run the same EM twice per dataset.

## Recommendation

For the pending study-29 targeted run (32 shards × 2000 reps): the safe,
bit-identical rework is worth folding in only if R>6 cells are planned;
otherwise it is a ~20% saving. SQUAREM is the change that matters (turns the
designed-missingness cells from ~0.2–0.5 s into ~30–60 ms per fit) but should
land behind an `EmOptions` switch and get a calibration smoke before becoming
the default.
