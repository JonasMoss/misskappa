# Alpha-missing simulation study setup

Design note for the paper-facing coefficient-alpha simulation. The active
backlog remains `dev/notes/todo.md`; this note records the actual grid so the
runner can be implemented once the normal-FIML alpha branch lands.

## Sources to borrow from

- Enders (2003): parallel thresholded-normal item batteries; `n = 100, 300,
  500`; 10 or 20 items; 3, 5, or 7 categories; half the items are targets for
  item-level missingness; target-item missing rates 15% and 30%; MCAR plus two
  MAR mechanisms. Useful for the missingness-rate and target-item framing, but
  too parallel-only for the main paper grid.
- Zhang and Yuan (2016): six-item one-factor calibrations. Tau-equivalent cell
  has standardized item variance, loading squared `0.60`, uniqueness `0.40`,
  and population alpha `0.90`. Non-tau-equivalent cell uses loading squared
  `0.20` for items 1-3 and `0.60` for items 4-6, with uniquenesses `0.80` and
  `0.40`; population alpha is about `0.778`. Their MAR mechanism keeps items 1
  and 4 complete and makes missingness in the remaining items depend on those
  observed anchors. This is the cleanest directly copied alpha setup.
- Savalei and Rhemtulla (2017): item-level missingness for composite models;
  `n = 200, 400, 600`; missing rates 5%, 15%, 30%; MCAR, MAR-linear, and
  MAR-nonlinear mechanisms using fully observed conditioning variables. Useful
  for the mechanism taxonomy and for target rates.
- van der Ark (2025): complete-data discrete/multinomial reliability SEs.
  Useful for the discrete scoring and coverage lens, but not a missing-data
  design.
- Existing misskappa experiments: use the anchor-MAR pattern from
  `experiments/studies/12-clean-mar-dgp/` and `experiments/probes/16-alpha-calibration-sweep/`
  instead of the old truth-dependent MAR mechanism.

## Measurement models

Use six to eight items. Six items keeps the categorical saturated FIML cell at
`5^6 = 15625` response patterns; eight items is continuous-only unless an
approximate categorical method is added.

Main continuous DGPs:

1. `tau6_essential`: six-item essential-tau one-factor model. Start from the
   Zhang-Yuan equal-loading calibration (`lambda^2 = 0.60`) but use mildly
   heterogeneous intercepts and residual variances, for example residual
   variances centered around `0.40`. Keep the exact Zhang-Yuan parallel cell
   (`lambda^2 = 0.60`, residual variance `0.40`, population alpha `0.90`) as a
   sanity/supplement row, not the only tau case.
2. `congeneric6_zy`: six-item non-tau one-factor model copied from Zhang and
   Yuan: loading squared `0.20` on items 1-3 and `0.60` on items 4-6, residual
   variances `0.80` and `0.40`. This gives a paper-backed congeneric
   comparator.
3. `congeneric8_gradient`: eight-item one-factor model with standardized item
   variances and loading sequence roughly `0.45, ..., 0.85`. This is the
   stronger congeneric stress case; normal FIML and pairwise only.
4. `twofactor8`: two correlated factors with four items per factor, factor
   correlation around `0.35`, and loading about `0.75`. This is explicitly
   non-congeneric/non-unidimensional. Alpha is still a covariance functional,
   but no reliability interpretation should be claimed.

Discrete/ordinal arm:

- Threshold `tau6_zy`, `congeneric6_zy`, and a six-item two-factor variant into
  five categories with symmetric thresholds. Report truth by quadrature for
  one-factor cells and by large complete-data Monte Carlo for the two-factor
  cell if two-dimensional quadrature is not worth implementing.
- Treat the earlier `paper5x6` bias as a known finite-sample limitation of
  saturated categorical FIML, not an implementation problem.

## Missingness mechanisms

Use complete data only as a sanity row, not as a headline mechanism.

1. `mcar15` / `mcar30`: Enders-style target-item MCAR. Choose half the items as
   targets and delete 15% or 30% of each target item independently. This makes
   the target rates comparable with the older alpha literature.
2. `mar_zy_light`: Zhang-Yuan six-item MAR. Items 1 and 4 are fully observed;
   missingness in items 5 and 6 depends on observed item 1, while missingness in
   items 2 and 3 depends on observed item 4. This yields light missingness
   (roughly 7% of cells) and is useful as a copied calibration check.
3. `mar_anchor15` / `mar_anchor30`: stronger clean MAR. Keep item 1 fully
   observed and let the other items' observation probabilities be logistic
   functions of item 1, tuned to the same 15% and 30% target rates. This is the
   main MAR cell because it is simple, ignorable from the observed item matrix,
   and already matches the kappa-paper clean-MAR style.
4. `mar_anchor_nonlinear30`: supplement/stress cell. Use a threshold or U-shaped
   function of item 1, following the Savalei-Rhemtulla linear/nonlinear split.

Avoid the old latent-truth MAR design; it is MNAR for the observed item matrix.

## Estimators

Paper-owned estimators:

- `pairwise`: pairwise-available covariance plug-in, with the new overlapping
  subsample SE once ported.
- `normal_fiml_sandwich`: `misskappa::alpha_continuous(se_type = "sandwich")`;
  main continuous MAR estimator.
- `normal_fiml_normal`: `misskappa::alpha_continuous(se_type = "normal")`;
  normal-theory comparator and diagnostic.
- `cat_fiml`: `misskappa::alpha(method = "fiml")`; discrete MAR estimator,
  capped to six five-category items in the main paper run.

Optional strawmen:

- `listwise`: useful because Enders and applied software discuss it.
- `pairwise_feldt_avg_n`: the average-pairwise-n Feldt interval as the
  incumbent software CI. Include only if cheap; this is a strawman, not a
  methods contribution.

Do not reproduce the full imputation zoo from Enders or Beland et al. unless a
reviewer asks. It dilutes the paper's inference question.

## Recommended grid

Pilot/smoke:

- DGPs: `tau6_essential`, `congeneric6_zy`.
- Mechanisms: `complete`, `mcar30`, `mar_anchor30`.
- `n = 200`; `B = 20`; continuous only, then ordinal with `B = 5`.

Main continuous run:

- DGPs: `tau6_essential`, `congeneric6_zy`, `congeneric8_gradient`,
  `twofactor8`.
- Distributions: normal for all cells; nonnormal only for `mcar30` and
  `mar_anchor30` as an SE robustness check.
- Mechanisms: `mcar15`, `mcar30`, `mar_zy_light`, `mar_anchor15`,
  `mar_anchor30`; `complete` as a sanity row.
- Sample sizes: `n = 200, 500, 1000`.
- Replications: start at `B = 500`; increase to `B = 1000` only for the final
  table cells if Monte Carlo SEs are still visible.

Main ordinal run:

- DGPs: `tau6_essential`, `congeneric6_zy`, six-item two-factor.
- Categories: five.
- Mechanisms: `mcar30`, `mar_anchor30`; `complete` as a sanity row.
- Sample sizes: `n = 250, 1000`.
- Replications: `B = 200` for the saturated categorical FIML cells. Larger
  `B` belongs on a scheduled long run, not an ordinary local pass.

Paper cut:

- Main text: continuous normal grid plus one ordinal table/figure.
- Supplement: `mar_zy_light`, nonnormal cells, listwise/Feldt strawmen,
  categorical timing/support diagnostics, and the known categorical FIML
  finite-sample bias caveat.

## Outputs

Each runner should write:

- `truth.csv`: DGP-level alpha truth, covariance summary, and truth method.
- `replicates.csv`: one row per replicate, method, DGP, mechanism, and `n`.
- `summary.csv`: bias, MC SD, mean SE, SE/SD ratio, Wald coverage, RMSE,
  failure rate, mean observed fraction, minimum pair count, and timing.
- `metadata.csv`: command, seed base, package versions, DGP grid, and method
  availability.

The report should lead with whether pairwise is calibrated under MCAR, biased
under MAR, and whether FIML buys enough efficiency/coverage to justify its cost.

## Implementation sequence

1. Wait for the normal-FIML alpha work to land cleanly in the package.
2. Add `experiments/studies/18-alpha-paper-simulation/` with a smoke path and the
   continuous runner first.
3. Add the ordinal/categorical arm after the continuous smoke is stable.
4. Run the capped paper pass, render the report, then migrate the paper-worthy
   CSVs/tables into the (private) alpha-missing paper repo once that scaffold
   is committed.
