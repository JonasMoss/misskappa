# Probe 44 — finite-sample performance of vector-valued quadratic kappa

**Question.** Probe 43 locked the *math* of the vector block-compound-symmetry (block-CS)
estimator on complete data at large n (projection == MLE; point consistency at n=8k; d=1
reduction; IF/vcov contract). This probe is the *finite-sample performance* story — the
vector analogue of probe 40 — for the three array estimators on a subjects×raters×features
input with squared vector loss `(x − y)' W (x − y)`:

- `pairwise` — pairwise-available covariance moment (MCAR; Conger/Fleiss only);
- `nt_fiml` — saturated normal-theory FIML (+ Brennan-Prediger);
- `nt_fiml_cs` — block compound symmetry `Σ = J_R ⊗ Ψ + I_R ⊗ E` (+ Brennan-Prediger).

Headline: under MCAR the block-CS projection **preserves** the contractions Conger/Fleiss
read off — `tB = sum((J_R⊗W)·Σ)` (sum of all blocks), `tT = sum((I_R⊗W)·Σ)` (sum of
diagonal blocks), `tG` (means only) — so `nt_fiml_cs` should target the *same* Conger/Fleiss
as `nt_fiml` even when block-CS is false. The open cell is MAR×misspecified, where the
constrained EM fixed point might not reproduce that invariance.

## Design (`parity_and_coverage.R`)

Probe-40 idiom; reuses probe 43's `make_sigma` / `project_block_cs` / `vec_coef`. `R = 4`,
free non-equal rater means. Screen mode (run here): DGPs {`blockcs` (block-CS true),
`camps` (saturated-but-not-block-CS — rater diagonal blocks scaled by `s=(.7,1,1.3,1.6)` +
a coupled (1,2) rater pair)}, `n ∈ {120,300}`, `d=2`, `W=identity`, mech
{complete, mcar(30%, rater-level), mar(rater-1 anchor)}, `B=300`. Rater-level amputation
(whole d-vector observed or missing). Conger/Fleiss truth is exact analytic from the true
`(μ,Σ)`; BP truth is a large-n Monte-Carlo (its chance baseline is data-range-dependent).
A self-test fits `nt_fiml` on n=60k and checks `coef ≈ vec_coef` (guards the rater-major /
feature-minor stacking) — passed, max |err| = 1e-4. All cells `conv = 1.000` (no FIML
failures, no pairwise overlap-guard at n≥120). Runtime 43s. Targeted mode (t6 non-Gaussian,
`d∈{1,2,3}`, coupling W, larger n/B) coded but not run.

## Results — realized consistency map (Conger, n=300, bias / cov_wald)

| | **`blockcs` (CS true)** | **`camps` (CS misspecified)** |
|---|---|---|
| **complete** | all three identical: −0.0012 / 0.95 | all three identical: −0.0027 / 0.95 |
| **MCAR** | pairwise −0.0017/0.91 · fiml −0.0003/0.97 · **cs −0.0003/0.96** | pairwise −0.0050/0.93 · fiml −0.0023/0.95 · **cs −0.0027/0.94** |
| **MAR** | **pairwise −0.0223/0.88** · fiml −0.0030/0.94 · cs −0.0022/0.94 | **pairwise −0.0293/0.85** · fiml −0.0032/0.94 · **cs +0.0039/0.91** |

### Findings

1. **Constraint-invariance holds (the headline).** On **complete** data the three estimators
   give *numerically identical* Conger/Fleiss (|nt_fiml − nt_fiml_cs| = 0 to 1e-16; pairwise
   too) — these coefficients are functions of the constraint-preserved contractions, so the
   CS constraint costs **nothing** in point accuracy. Under **MCAR** `nt_fiml` and
   `nt_fiml_cs` stay near-identical (|diff| ~1e-3, shrinking with n) and **consistent even on
   `camps` where block-CS is false** (cs bias −0.0027, cov 0.94). The MCAR projection-
   invariance prediction is confirmed.

2. **No efficiency gain from the correct constraint** — contrary to the going-in conjecture.
   Because Conger/Fleiss read off constraint-invariant contractions, `nt_fiml_cs` and
   `nt_fiml` have the *same* sd/RMSE when block-CS is true (sd 0.0224 vs 0.0224 at
   blockcs/mcar/n300). The CS route's payoff is **not** efficiency for these coefficients —
   it is (a) needing only a *connected* (not complete) co-observation graph and (b) numerical
   stability (block-averaged Σ, no singular per-pattern inversions). Both FIML paths are more
   efficient than `pairwise` under MCAR (sd 0.022 vs 0.030), which uses only co-observed rows.

3. **`pairwise` is MAR-inconsistent** (as expected, it is an MCAR moment estimator): bias
   −0.022 (blockcs) / −0.029 (camps), coverage collapsing to 0.85–0.88 and *worsening* with n
   — the inconsistency signature. Both FIML methods stay MAR-consistent (|bias| < 0.004).

4. **The MAR×misspecified cell is essentially benign here.** `nt_fiml_cs` on `camps`/MAR has
   bias +0.0039 at n=300 (vs saturated −0.0032; the two FIML methods straddle zero, |diff|
   ~0.007 — the largest gap anywhere) and does **not** clearly grow with n (0.0033→0.0039).
   Coverage 0.91, driven by a slightly optimistic SE (se_ratio 0.91, shared with saturated).
   So at this misspecification strength we do **not** reproduce the strong MAR+misspec
   inconsistency probe 40 found for variance-component *alpha* (bias −0.018, cov 0.87). Two
   readings, not separated by this run: (i) the vector Conger/Fleiss ratio is more robust to
   covariance misspecification than alpha; (ii) the `camps`×anchor-MAR pairing is mild.
   **Follow-up:** a stronger / adversarial `camps` (and a MAR mechanism aimed at the
   contraction directions the CS model gets most wrong) in targeted mode.

5. **Brennan-Prediger SE is optimistic — the one weak spot.** BP se_ratio 0.29–0.53 and 0%
   coverage against the large-n MC truth. Two compounding causes: (a) **the analytic BP SE
   ignores the sampling variability of the observed-range chance baseline**
   `D_e^BP = Σ_a W_aa(hi_a−lo_a)²/6` (se_ratio ≪ 1 is truth-independent, so this is real, not
   a truth artifact); (b) the chance baseline drifts with n, so a fixed large-n truth is not
   the finite-n estimand. Conger/Fleiss are immune (range-free). BP needs either a
   range-aware delta-method SE or an explicit fixed working range before its CI can be trusted.

## Takeaways for the package / paper

- For Conger/Fleiss, `nt_fiml_cs` is a free win over saturated `nt_fiml`: identical accuracy
  under complete/MCAR, MAR-consistent, but works on connected-only graphs and is more stable.
  Frame its value as **robustness/identifiability, not efficiency**.
- `pairwise` (vector) is MCAR-only — do not use under suspected MAR.
- The vector **Brennan-Prediger SE is not yet trustworthy** (ignores chance-baseline sampling
  variation); flag before exposing BP CIs, or fix with a range-aware SE.

## Resolution of the BP weak spot (2026-06-15)

Finding 5's root cause was an interface gap, not a math gap: the scalar path lets the user
fix the chance baseline via `values` (`quadratic_bp_chance`), but the vector path
(`vector_bp_chance`, new in 1.2.0) was implemented without that hook and always used the
observed per-feature range. Now wired through: `kappa(arr, estimator="nt_fiml"/"nt_fiml_cs",
W=W, ranges=<features×2 lo/hi>)` fixes the per-feature width, so the BP chance denominator is
deterministic — the SE becomes exact and the estimand stops drifting with n. Verified:
observed-range BP drifts 0.605→0.703→0.758 over n=60→240→1000, fixed-range BP is stable
0.708→0.743→0.748 (sampling noise only). Passing the observed range reproduces the old
default exactly. A re-run of the BP rows with `ranges=` set (so the BP truth is a fixed
population value, not a large-n MC) is the clean way to assess BP coverage — deferred to the
targeted run.

## Stress study (`stress.R`, 2026-06-15)

Three targeted regimes (mclapply-parallel) chasing what the screen run left open.

### Regime 1 — adversarial MAR (this overturns the screen "benign" verdict)

> **PARTIALLY SUPERSEDED (2026-06-15) — see "Regime 1 REBUILD" below.** The
> `camps_sigma` knob drives diagonal rater blocks negative for `s > 1.67`, so it is
> PSD-floored into a near-singular truth at `s ≥ 2`. The CS bias curve at `s ≤ 1` is
> sound, but the "saturated fails to converge at `s ≥ 2`" finding is a **DGP-degeneracy
> artifact**, not a strong-covariance effect. The well-conditioned rebuild
> (`regime1_wc.R`) shows saturated `nt_fiml` converges (`conv = 1`) at ALL `s`.

Two-camp `camps` covariance with a tunable misspecification strength `s` (`s=0` is
exact block-CS; camps {1,2} and {3,4} get unequal rater variances + extra within-camp
coupling), plus a **camp-aimed MAR** (anchor camp A, drop camp B on the anchor value, so
the cross-camp blocks the CS model averages away are under-observed in a value-dependent way).

- **MCAR control — contraction-invariance confirmed at ANY strength.** `nt_fiml_cs` Conger
  bias stays `|·| ≤ 6e-4` and coverage `0.94–0.96` for every `s ∈ {0,0.5,1,2,4}` at n=1500.
  Under MCAR the block-CS projection is consistent no matter how false block-CS is — the
  screen-run prediction, now stress-tested to extreme misspecification.
- **MAR — `nt_fiml_cs` is genuinely inconsistent, severity grows with `s` (the headline).**
  Conger bias and coverage (cov_wald) under MAR:

  | s | bias n=400 | bias n=1500 | cov n=400 | cov n=1500 |
  |---|---|---|---|---|
  | 0   | −0.002 | −0.001 | 0.94 | 0.95 |
  | 0.5 | +0.008 | +0.011 | 0.92 | **0.78** |
  | 1   | +0.023 | +0.022 | 0.74 | **0.34** |
  | 2   | +0.030 | +0.031 | 0.60 | **0.12** |

  Bias persists (does not shrink with n) and coverage **collapses as n grows** — the
  asymptotic-inconsistency signature. `s=1` is the cleanest case (non-degenerate PD
  covariance, stable +0.022 bias, coverage 0.74→0.34). So the screen run's near-zero
  `camps` bias was an artifact of mild misspecification; the real MAR×misspec failure mode
  exists and matches probe 39/40's "benign→severe" map for the vector coefficients.
- **Saturated `nt_fiml` is the consistent fix — but turns fragile.** It stays unbiased
  (bias ~0, cov ~0.95) at `s ≤ 1`, but **fails to converge (NA)** at `s ≥ 2`: the strong
  camps covariance makes the saturated per-pattern blocks singular. (Caveat: at `s=4` the DGP
  itself is near-singular and PSD-floored, so read `s≥2` as "extreme".) Net: under MAR there
  is a real **consistency-vs-stability tradeoff** — CS is stable but biased off-model;
  saturated is consistent but numerically fragile exactly where the structure is strong.
  **[SUPERSEDED — the `s ≥ 2` non-convergence and hence the "consistency-vs-stability
  tradeoff" are artifacts of the PSD-floored degenerate truth. On the well-conditioned
  rebuild (below) saturated converges at ALL `s` and stays consistent + calibrated: under MAR
  on a non-degenerate truth, saturated `nt_fiml` is the all-around answer, no tradeoff.]**

### Regime 2 — robustness / conv-rate (the selling point the screen run couldn't show)

Clean block-CS DGP, small R + sparse/designed missingness; **convergence rate is the metric**.

| design | `nt_fiml` (saturated) | `pairwise` | `nt_fiml_cs` |
|---|---|---|---|
| `mcar70` (70% rater-missing) | collapses: conv **0.00–0.40** | sags at n=40 (0.77–0.93) | **0.95–1.00** |
| `random2` (2 raters/subject) | fails small-n/large-R (**0.01** at R=5,n=40) | 0.86 at R=5,n=40 | **1.00** |

`nt_fiml_cs` (connected-graph-only, block-averaged Σ) is the lone survivor where the
saturated per-pattern inversions go singular. This is the concrete payoff of "needs only a
connected co-observation graph," and it is invisible at the screen run's n≥120 (all conv=1).

### Regime 3 — fixed-range Brennan-Prediger calibration (the fix confirmed)

With `ranges=` set, the BP chance baseline is deterministic so BP has an analytic population
value. BP is now **well-calibrated everywhere** (both W, all mechanisms, n∈{120,300,1000}):
bias `≤ 1e-3`, `se_ratio 0.95–1.04`, `cov_wald 0.93–0.97`. Compare the screen run's
observed-range BP (`se_ratio 0.29–0.53`, coverage **0.00**). The `ranges=` hook fully fixes
the one BP weak spot — BP is trustworthy once the range is supplied.

### Stress takeaways for the paper

- **MCAR:** `nt_fiml_cs` is the recommendation — consistent at any misspecification,
  connected-graph-only, numerically robust, BP calibrated with `ranges=`.
- **MAR:** saturated `nt_fiml` is the recommendation — on a well-conditioned truth it is
  consistent, stable (`conv = 1`), and calibrated at ANY misspecification strength (see the
  rebuild below; the old "fails to converge / consistency-vs-stability tradeoff" was a
  degeneracy artifact). `nt_fiml_cs` is graded-inconsistent under MAR (bias +0.01→+0.05,
  coverage collapses with n) and should be used only when block-CS is approximately right OR
  the design is too sparse for saturated to fit (the regime-2 corner). Frame honestly
  (mirrors the structured-ICC story in probe 39).
- Artifacts: `stress_regime{1,2,3}.rds` + `stress_regime1_wc.rds`. Re-run:
  `Rscript .../stress.R [regimes…]`, `Rscript .../regime1_wc.R`.

## Regime 1 REBUILD — well-conditioned adversarial knob (`regime1_wc.R`, 2026-06-15)

Disentangles the two things the original `camps` knob confounded: (i) `nt_fiml_cs`'s MAR
bias from block-CS misspecification, and (ii) saturated `nt_fiml` "non-convergence" — which
at `s ≥ 2` was really *the truth is degenerate* (PSD-floored), not *missingness broke the EM*.

**New DGP — latent camp factor, PD at every `s`.** `Y_ij = μ_j + b_i + c_i^{camp(j)} + ε_ij`
with `b_i~N(0,Ψ)`, `c_i^{camp}~N(0,sG)`, `ε~N(0,E)` ⇒
`Σ = J_R⊗Ψ + I_R⊗E + Σ_camp (1_camp 1_camp')⊗(sG)`. Diagonal blocks `= Ψ + sG + E` (constant),
within-camp off-diag `= Ψ + sG`, cross-camp off-diag `= Ψ`. `s=0` is exact block-CS; `s>0`
makes it false in a graded way. **Smallest eigenvalue is PINNED at `λ_min(E)=0.4382` for all
`s`** (rcond 0.074→0.036 over `s∈{0,…,4}`), so every per-pattern `Σ_oo ⪰ I⊗E ≻ 0` — the
saturated EM's inversions are never population-singular. Knob is provably clean. `B=300`,
`n∈{400,1500}`, `R=4`, `d=2`, `W=I`, mech {mcar 30%, camp-aimed mar}.

**Finding A — clean graded CS inconsistency under MAR (Conger):**

| s | bias n=400 | bias n=1500 | cov n=400 | cov n=1500 |
|---|---|---|---|---|
| 0   | −0.001 | −0.001 | 0.95 | 0.96 |
| 0.5 | +0.012 | +0.010 | 0.90 | 0.83 |
| 1   | +0.018 | +0.018 | 0.86 | 0.61 |
| 2   | +0.032 | +0.033 | 0.70 | 0.19 |
| 4   | +0.049 | +0.045 | 0.49 | **0.07** |

Bias persists (n-stable), coverage collapses with n — asymptotic-inconsistency signature, now
monotone across the FULL `s` range with no degeneracy confound (the old run could only trust
`s ≤ 1`).

**Finding B — the disentangling headline. Saturated `nt_fiml` converges at ALL `s`.**
`conv = 1.000` for every `s`, both mechanisms, both n; MAR bias `|·| ≤ 0.0023`, coverage
0.92–0.96 throughout. So the original "saturated fails at `s ≥ 2`" was **entirely a
DGP-degeneracy artifact**. On a well-conditioned truth there is **no consistency-vs-stability
tradeoff under MAR**: saturated is consistent, stable, AND calibrated even at extreme
misspecification. Saturated's only genuine failure is regime-2 *sample* singularity
(sparse/small-n), independent of the population conditioning — that is where CS's stability
edge actually lives.

**Finding C — contraction-invariance re-confirmed to `s=4`.** MCAR `nt_fiml_cs` bias
`|·| ≤ 0.0028`, coverage 0.93–0.97 at every `s`: the block-CS projection is consistent no
matter how false block-CS is, now stress-tested on a non-degenerate truth.

**Consequence — ridge rescue (PASSOFF sim-mod 2) is moot, not built.** Its motivation was
"save saturated from strong-covariance non-convergence," which Finding B shows does not occur.
The residual sparse corner (regime 2) is owned by CS already (consistent + stable there); a
ridge would only help the narrow MAR+sparse+misspec triple-corner, at the SE-optimism cost the
scalar precedent already documented (`project_nt_fiml_degenerate_sigma`, study 33: rescue
recovers failures but coverage drops to ~0.35, smoked + rejected).

**Revised paper framing (split by mechanism, both coefficients):**
- **MCAR** → recommend `nt_fiml_cs`: consistent at any misspecification, connected-graph-only,
  numerically robust, no efficiency cost vs saturated (constraint-invariant contractions),
  BP-calibrated with `ranges=`.
- **MAR** → recommend saturated `nt_fiml`: consistent, stable, calibrated at any
  well-conditioned misspecification. `nt_fiml_cs` is graded-inconsistent (Finding A) — use it
  under MAR only if block-CS is ≈right or the design is too sparse for saturated to fit.

## Artifacts

`parity_and_coverage.R` (screen + targeted), `coverage_summary.{rds,csv}` (full summary),
`replicates.csv` (per-replicate). Re-run: `Rscript .../parity_and_coverage.R [--targeted]`.
