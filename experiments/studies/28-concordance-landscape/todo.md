# Study 28 — follow-ups (CCC landscape / U-statistic connection)

Findings from the 2026-06-09 exploration. Ordered by value.

Status: **A, B, C done** (report rewritten, 4 new benchmark tables, re-rendered;
9 tables total). D in report. E = paywalled DOIs logged. F optional (paper frame).
Integration of this session's findings is complete; remaining items are external
(grab papers / write a paper), not study edits.

## A. Fold the U-statistic / `cccUst` findings into the report  [DONE]
The closest prior art to misskappa's quadratic/CCC path is the King–Chinchilli–
Carrasco U-statistic CCC (`cccrm::cccUst`), not anything on the kappa side. Add:

- **The `delta` dial** — `cccUst(delta=1)` = Lin CCC = misskappa quadratic Conger;
  `cccUst(delta=0)` = unweighted Cohen kappa. Verified to 5 dp vs misskappa / irr.
  One estimator spans the categorical↔continuous spectrum (= the van Oest–Moss
  bridge, built into software from the CCC side since King & Chinchilli 2001).
- **Repeated-measures mechanism** — each subject is a time-*profile* per method;
  `phi` gives WITHIN (matched) vs BETWEEN (chance) quadratic-form distances
  `v^T Dmat v`, `v = |X-Y|^delta`; CCC = 1 - U/V; subject = i.i.d. unit; the
  per-subject U-statistic projection IS the influence function.
- **Verified identity** — `cccUst(delta=1, Dmat=diag(w))` == misskappa vector
  quadratic kappa with `feature_weights=w`, to 1e-16. So **Dmat (time weights) =
  W (component weights)** and **time-profile = rating vector**: `cccUst`'s repeated
  CCC and misskappa's vector path are the same construction from opposite ends.
  Add as a new benchmark row + reframe the "repeated measurement" section.
- **Positioning paragraph** — KCC own the unified, distribution-free,
  repeated-measures CCC *estimator*; misskappa adds (i) the missing-data layer
  (IPW/FIML — `cccUst` is listwise only), (ii) influence-function inference +
  Wald/joint equality tests (they reach SEs by U-stat delta / bootstrap), (iii)
  the kappa=CCC *identity* (they had it as an estimator, not a theorem), (iv)
  R>2 raters and categorical via the nominal loss.

## B. Vanbelle cluster-collapse vs misskappa-vector note  [DONE]
- State the homogeneity assumption precisely: within a cluster (patient) all
  objects (sites) are identically distributed; her own Table 2 rejects it
  (p<0.0001), yet she pools anyway for small-n reasons.
- Verified result: misskappa vector == Vanbelle pooled **iff** homogeneity holds;
  otherwise they differ by `sum_j Cov_site(p_j(r1,.), p_j(r2,.))` in the chance
  term (misskappa chance per-site, Vanbelle pooled). 3-row sim confirms it
  (homogeneous diff ~7e-4 = MC noise; heterogeneous diff -0.046; real EXP
  0.519 vs 0.563, and 0.563 reproduces her Table 3 "All").

## C. Signed vs absolute cross-terms (methodological)  [DONE — in A's identity table]
- misskappa's full-W quadratic loss `(x-y)^T W (x-y)` is the weighted-Euclidean /
  Mahalanobis distance `||L(x-y)||^2` — the principled object. `cccUst` uses
  absolute cross-terms `sum_st D_st |x_s-y_s||x_t-y_t|`, not metric-induced, only a
  real distance when diagonal. Matters only off-diagonal — which is exactly the
  regime misskappa uses (exp-26 site-W couples left/right at 0.25) and the CCC
  literature never does. Document as a one-paragraph distinction (misskappa is the
  correct generalization off the diagonal).

## D. Unbiased-CCC future work (already in report, keep)
A real fix is a second-order influence-function / jackknife bias correction of the
ratio, NOT the denominator patch (removes only ~20% of the bias; jackknife ~90%+).
Open work: closed-form 2nd-order IF for the weighted multirater CCC so SE stays
consistent. Low priority (sub-percent above n~30).

## E. Grab the priority papers (institutional access)
Confirmed paywalled (Europe PMC), DOIs in `dev/refs/concordance-refs.md`:
King & Chinchilli 2001 (10.1002/sim.845); King–Chinchilli–Carrasco 2007
(10.1002/sim.2778); Carrasco et al. 2013 (10.1016/j.cmpb.2012.09.002).

## F. Paper frame (quadratic / Conger's kappa -> Statistics in Medicine)  [NOTE WRITTEN]
Full framing in `paper-framing.md`: dual CCC + psychometrics bridge; repeated
measurements via a PD weight W (not cccUst's absolute form); Conger/Fleiss/BP;
careful NT-FIML exposition for biostat; CRACKLES example + Vanbelle contrast;
extend exp-23 sims to repeated measures. Gating task: related-work pass on the
King-Chinchilli-Carrasco CCC line before submitting (a SiM reviewer knows it cold).
