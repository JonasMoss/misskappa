# Paper-C (quadratic / Conger's kappa) — framing note for a biostatistics submission

Working note seeded by the 2026-06-09 CCC exploration (study 28). Lives here next
to the study that motivated it; move to the private `papers/00_wip/quadratic/`
folder if it should travel with the manuscript. Naming throughout stays
"Conger's kappa" (no Cohen/Conger hedge, no "quadratically weighted" prefix in
body prose); Brennan–Prediger treated as before; Fleiss' kappa included.

## Target and one-line pitch

One-liner: *the quadratically weighted multirater agreement coefficient is Lin's
concordance correlation coefficient, generalized to R > 2 raters and repeated
measurements, and — the genuine novelty — estimated under missing data with
principled IPW/FIML and influence-function inference.*

**Venue — write agnostic, decide late** (the framing is venue-robust; only
theory-vs-applied emphasis and length shift). Shortlist by character:
- **Biometrics** — the *home lineage* of the CCC methods (Lin 1989, Carrasco &
  Jover 2003, Barnhart 2002 all here) + higher prestige than SiM. Aim here if the
  draft is theory-forward (identity-as-theorem, missing-data estimation theory,
  IF inference). The prestige + lineage play.
- **Statistics in Medicine** — bigger reach, faster, more variable quality (long
  tail of a high-volume journal); where King–Chinchilli + the repeated-measures
  CCC live. Safe fit / reach fallback.
- **SMMR / Biometrical Journal** — the agreement / repeated-measures / multilevel
  niche (Vanbelle's own home). Natural for the repeated-measures angle.
- **Biostatistics (Oxford)** — high methodological bar, selective; only if novelty
  is sharp. *Not* the same prestige as SiM — more selective.
- **Incentive structure:** check the Norwegian Kanalregister level before deciding
  — Biometrics / Biostatistics likely **level 2**, SiM probably level 1; if level-2
  credit matters this nudges toward Biometrics even at a higher bar / slower review.
- Net lean: **Biometrics if theory-forward and it ranks higher in the incentive
  system; SiM/SMMR as the reach/fit fallback.**

## Framing (the dual bridge)

- **Connect to CCC *and* to psychometrics.** The paper is the bridge: Conger's
  kappa with quadratic weights = Lin's CCC (van Oest & Moss). Speak to the
  method-comparison/biostat world (CCC) and the agreement-coefficient world
  (kappa) at once. The CCC connection is what unlocks the biostat venue; the
  psychometric grounding (chance-corrected agreement, Fréchet variances) is the
  home theory.
- **Repeated measurements are in from the start.** The vector / component-
  separable coefficient is the King–Chinchilli–Carrasco repeated-measures
  U-statistic CCC (study 28 verified: identical to `cccrm::cccUst` to 1e-16 for
  a diagonal weight; `Dmat` = the component-weight matrix W). Frame "occasions
  within subject" as vector components with a weight matrix — that single device
  covers replicates, time points, and clustered sites.

## Scope decisions (per discussion)

- **Use an ordinary PD weight only.** The weighted-Euclidean / Mahalanobis form
  `(x−y)ᵀ W (x−y) = ‖L(x−y)‖²`, W positive-definite. *Not* `cccUst`'s
  absolute-cross-term object `Σ D_st |x_s−y_s||x_t−y_t|`, which isn't induced by
  any metric and discards the sign of joint errors. They coincide on the
  diagonal; off-diagonal, ours is the correct generalization (and misskappa
  already rejects non-PD W — correct behavior). Off-diagonal W is meaningful here
  (CRACKLES couples left/right of each site) and essentially unused in the CCC
  literature, so it's a clean point of difference, not a corner case.
- **Coefficients:** Conger's kappa (lead), **Fleiss' kappa** (the pooled-marginal
  multirater form), and **Brennan–Prediger** (uniform baseline) — treated as in
  the prior paper.
- **Normal-theory FIML needs a more careful introduction than for Psychometrika.**
  EM / MAR / FIML are *familiar to biostatisticians* (Little & Rubin, Laird &
  Ware) — so do **not** pitch FIML as exotic; pitch it as "the standard
  missing-data toolkit, finally brought to agreement coefficients, which have sat
  on listwise deletion for 30 years." BUT the *specific construction* — EM over a
  saturated multivariate-normal covariance feeding a quadratic agreement
  functional, with robust delta-method / influence-function SEs — is unfamiliar
  *in this form*; spell it out (the EM steps, what's being estimated, the IF
  variance), don't assume it the way a Psychometrika draft could. Understatement
  plays better with this crowd than fanfare.

## Related-work imperative (the prerequisite before submitting)

A SiM reviewer knows the CCC line **cold**; under-citing it = desk-reject risk.
We *found* this connection mid-exploration, which means the current draft likely
under-engages it. Engage head-on, then position as **unification + what's new**:

- Cite/discuss: Lin 1989/2000 (CCC, TDI); **King & Chinchilli 2001** (generalized
  CCC, the kappa↔CCC `delta` dial); **King–Chinchilli–Carrasco 2007**
  (repeated-measures U-statistic CCC); **Carrasco & Jover 2003** (CCC via variance
  components); Carrasco 2013 (`cccrm`); Barnhart–Haber–Song 2002 (overall CCC),
  Barnhart–Haber–Lin 2007 (the scaled/unscaled taxonomy). DOIs in
  `dev/refs/concordance-refs.md`; all paywalled (institutional access).
- **Do not present "quadratic kappa = CCC" as the novelty** — the generalized CCC
  already spans categorical↔continuous (their `delta` dial proves it). The
  novelty is the four deltas below.

## What is genuinely new (lead with these)

1. **Missing data.** IPW (MCAR) and NT-FIML (MAR) for the weighted multirater
   coefficient. The entire CCC/method-comparison software stack is **listwise
   only** — this is the real contribution.
2. **Influence-function inference + equality tests** (Wald / joint), where the
   CCC packages reach SEs by U-statistic delta or bootstrap.
3. **The identity stated as a theorem with interpretations** (van Oest & Moss) —
   they had it as an estimator, not a result.
4. **R > 2 raters and categorical data** in one framework (Conger/Fleiss/BP).

## Worked example

- Use **`dat.vanbelle2019` (CRACKLES)** — the repeated-measures / vector example.
- Natural discussion point: contrast with **Vanbelle's own multilevel kappa**.
  She *collapses* the six sites under a within-cluster homogeneity assumption
  (her Table 2 rejects it at p<0.0001); we *keep* them as a vector. Study 28
  shows the two coincide **iff** homogeneity holds, differing otherwise by the
  between-site covariance of the raters' marginals in the chance term. Our
  per-site chance correction is the more honest one; and our vector SE already
  clusters at the patient level by construction (no separate multilevel
  correction needed).

## Simulations

- **Extend the quadratic simulation study (exp 23) to repeated measurements:**
  subjects × methods/raters × occasions, a PD weight W over occasions, MCAR and
  MAR amputation of occasions, comparing IPW / NT-FIML / listwise. This is the
  cell the CCC literature can't run (they're complete-data), so it's where the
  paper's contribution is demonstrated.
- Keep the existing two-rater / scalar coverage cells as the base case.

## Scope out / don't oversell (notes to self)

- **Don't build an "unbiased CCC."** Study 28 showed the denominator/location
  patch (and Carrasco's VC route) does ~nothing to the point estimate and isn't a
  proper correction anyway (removes ~20% of the O(1/n) bias; the real bias is the
  ratio nonlinearity, which a jackknife / 2nd-order IF removes ~fully). At most a
  one-paragraph remark; not a contribution.
- **Acknowledge but scope out the unscaled/interval family** — Bland–Altman LoA,
  TDI/CP, Deming/Passing–Bablok. That's the other half of method-comparison; we
  do the scaled/index side. One sentence so the reviewer sees we know it exists.
- **ICC:** a single CCC↔ICC pointer (Nickerson 1997), no ICC zoo.

## Tomorrow's first moves

1. Related-work pass on the CCC line (the gating task — read KCC 2001/2007,
   Carrasco & Jover 2003).
2. Re-draft the intro around the CCC↔psychometrics bridge + repeated measures.
3. Sketch the careful NT-FIML exposition for a biostat reader.
4. Spec the repeated-measures simulation cells (extend exp 23).
