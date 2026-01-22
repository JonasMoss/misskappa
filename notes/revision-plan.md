# Revision plan (Psychometrika target)

This is an internal working plan for smoothing the rough draft while keeping scope tight and length modest.

## Snapshot (current draft)

- PDF: `kappa-missing.pdf` is 6 pages and ~2800 words.
- Introduction length: ~326 words (from PDF text extraction).
- References in `kappa-missing.bib`: 7 entries.
- Obvious placeholders: Sections 4 (Simulations), 6 (R package), 7 (Concluding remarks) are currently empty in the PDF.
- One visible truncation/incomplete sentence in the Introduction: “We w”.

## Targets (keep it short + clean)

- Main paper length: aim ~10–12 pages (appendix can be separate/online if needed).
- Introduction: keep ~300–600 words (specialist journal; no padding).
- Citations: aim for “enough to signal competence” rather than exhaustive.
  - Minimum target: 15–25 total references, including:
    - (i) core kappa definitions (Cohen/Fleiss/Brennan–Prediger),
    - (ii) missing-data foundations (MCAR/MAR, IPW/Hájek),
    - (iii) prior work specifically on kappas with missing data / practical software (e.g. Gwet / irrCAC ecosystem),
    - (iv) inference references (U-statistics / delta method / asymptotic variance).
- Reader experience: an applied reader should be able to (a) choose an estimator, (b) compute a standard error, and (c) understand the assumptions without parsing every proof.

## Hard TODOs (must do)

- [ ] Fill in Section 4 (Simulations) with a minimal but convincing design:
  - [ ] 2–3 DGPs (at least one where available-case fails; at least one where it works).
  - [ ] Compare: available-case vs IPW (and any third estimator only if it earns its keep).
  - [ ] Output: 1 table + 1 figure max (unless Psychometrika format pushes otherwise).
  - [ ] Add reproducibility notes (seeds, n, R, C, missingness mechanism).
- [ ] Fill in Section 6 (R package misskappa):
  - [ ] 1-paragraph capability summary.
  - [ ] 10–20 line example (raw ratings + missing, show `coef()` + `vcov()` + CI).
  - [ ] Mention datasets included (e.g. Fleiss 1971 example already used in Section 5).
- [ ] Fill in Section 7 (Concluding remarks):
  - [ ] 1 short paragraph: main takeaway + practical guidance.
  - [ ] 1 short paragraph: what’s next (EM/ML efficiency story) without promising too much.
- [ ] Fix obvious draft roughness:
  - [ ] Complete/replace the “We w” sentence in the Introduction.
  - [ ] Tighten the abstract phrasing (it currently claims sims + package; align with what’s actually delivered).
- [ ] Readability pass on the inference section (2.2):
  - [ ] Add a “recipe paragraph” before influence-function algebra: what the user computes in practice.
  - [ ] Make notation local: every symbol used in (2.2) defined immediately nearby.
  - [ ] Ensure the reader understands: what is estimated vs what is target (κ vs κIPW target vs κ under MCAR/PMCAR).

## Soft TODOs (nice to have)

- [ ] Add 1–2 sentences up front clarifying the missingness assumption vocabulary (MCAR vs pairwise MCAR vs “missing uniformly at random”).
- [ ] Add a very small “narrative” motivation in the Introduction (one paragraph max): why available-case is tempting, when it breaks, why IPW helps.
- [ ] Notation cleanup:
  - [ ] Ensure consistent use of `X⋆` (full data) vs `X` (observed) vs `M` (missingness).
  - [ ] Ensure losses `l(·,·)` are introduced once and then referenced consistently.
- [ ] Citation pass: add missingness and IPW references beyond Tsiatis (only if genuinely needed).
- [ ] Decide how prominently to feature Gwet/irrCAC:
  - [ ] Keep it as a comparison baseline (short) unless it is central to the story.

## “Inference theory is readable” checklist

- [ ] Each theorem/proposition has a 1–2 sentence “what this means” immediately after.
- [ ] There is at least one worked mini-example that uses the final variance output (even if trivial).
- [ ] The appendix contains the heavier derivations; the main text keeps only the minimal formulas needed for implementation.
- [ ] Every estimator has: assumption summary + consistency statement + how SE is computed.

## Suggested sequence (low churn)

1. Fix “obvious roughness” (broken sentences, missing section stubs).
2. Draft the Simulations section around a minimal design (don’t overbuild).
3. Do the inference readability pass (add recipe + tighten notation).
4. Add citations in one dedicated pass near the end.
5. Final polish: “what is the contribution?” and “what is not covered?”

