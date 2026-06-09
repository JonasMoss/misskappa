# Experiment Index

The experiment tree is organized by lifecycle rather than topic. Use this file
as the source of truth for canonical paths and maintenance status.

## Lifecycle Buckets

- `studies/` contains reusable or paper-facing studies with some longevity
  promise.
- `workbench/` contains active work in progress.
- `probes/` contains diagnostics, stress tests, pilots, and validation checks.
- `archive/pre-redesign/` contains frozen records whose runners are not
  maintained against the current R API.

Always reference experiments by their canonical lifecycle path (the table
below). There are no root-level alias directories.

## Catalog

| Experiment | Canonical path | Status | Topic | Feeds | Current API | Results policy | Finding |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 01 coverage / IIF / Louis | `experiments/archive/pre-redesign/01-coverage-iif-louis/` | frozen | FIML/IPW coverage | Papers A/B | frozen | keep summaries | IPW coverage was near nominal under MCAR; Louis-SE behavior motivated the later spectrum diagnostic. |
| 02 rater model sensitivity | `experiments/archive/pre-redesign/02-rater-model-sensitivity/` | frozen | DGP sensitivity | Papers A/B | frozen | keep summaries | Paper simulation orderings depend on the rater model choice. |
| 03 AC vs IPW efficiency | `experiments/studies/03-ac-vs-ipw-efficiency/` | reusable study | IPW/available-case | Paper A | frozen | keep summaries | AC is lower-variance when consistent, but biased under non-exchangeable rater observation rates. |
| 04 counts sampling misspecification | `experiments/archive/pre-redesign/04-counts-sampling-misspec/` | frozen | counts FIML | Paper B | frozen | keep summaries | Aggregating heterogeneous raw ratings to counts before FIML can induce large bias. |
| 05 FIML sparsity scaling | `experiments/archive/pre-redesign/05-fiml-sparsity-scaling/` | frozen | saturated FIML scaling | Paper B | frozen | keep summaries | High-dimensional saturated FIML needs explicit finite-sample caveats. |
| 06 raw estimator scaling | `experiments/probes/06-raw-estimator-scaling/` | probe | engineering timing | library | frozen | keep summaries | Raw estimator timing record for engineering decisions. |
| 07 quadratic Edgeworth coverage | `experiments/archive/pre-redesign/07-quadratic-edgeworth-coverage/` | frozen | complete-data quadratic coverage | Paper C | frozen | keep summaries | Edgeworth corrections did not become the recommended main route. |
| 08 quadratic bootstrap literature | `experiments/archive/pre-redesign/08-quadratic-bootstrap-literature/` | frozen | quadratic bootstrap | Paper C | frozen | keep summaries | Bootstrap alternatives were explored for the quadratic coverage story. |
| 09 joint vcov pilot | `experiments/probes/09-joint-vcov-pilot/` | probe | joint inference | Papers A/C | frozen | keep summaries | Influence-function joint covariance supports Hausman-style and homogeneity tests. |
| 10 Louis spectrum | `experiments/probes/10-louis-spectrum/` | probe | FIML variance diagnostic | Paper B | frozen | keep summaries | Rank-truncating the Louis information inverse removes near-kernel variance blow-up. |
| 11 quadratic delta variance | `experiments/archive/pre-redesign/11-quadratic-delta-variance/` | frozen | quadratic variance | Paper C | frozen | keep summaries | Analytic variance route was explored as a no-bootstrap option. |
| 12 clean MAR DGP | `experiments/studies/12-clean-mar-dgp/` | reusable study | MAR DGP design | Paper B | frozen | keep summaries | Anchor-MAR and sequential MAR designs give clean non-FIML bias contrasts. |
| 12 quadratic rare disagreement | `experiments/archive/pre-redesign/12-quadratic-rare-disagreement/` | frozen | quadratic boundary behavior | Paper C | frozen | keep summaries | Rare-disagreement boundary diagnostics informed quadratic coverage recommendations. |
| 13 quadratic item hurdle | `experiments/archive/pre-redesign/13-quadratic-item-hurdle/` | frozen | item-level quadratic boundary | Paper C | frozen | keep summaries | Item-level hurdle diagnostics were useful but too blunt as final interval machinery. |
| 14 quadratic MCAR verification | `experiments/probes/14-quadratic-mcar-verification/` | probe | MCAR verification | Paper C | frozen | keep summaries | Pairwise quadratic behavior was checked under MCAR. |
| 15 alpha categorical smoke | `experiments/probes/15-alpha-categorical-smoke/` | probe | categorical alpha | alpha-missing | frozen | keep summaries | First mechanical check for categorical alpha FIML feasibility. |
| 16 alpha calibration sweep | `experiments/probes/16-alpha-calibration-sweep/` | probe | categorical alpha calibration | alpha-missing | frozen | keep summaries | Sparse high-dimensional categorical FIML shows finite-sample downward bias. |
| 17 alpha normal-FIML validation | `experiments/probes/17-alpha-fiml-normal-validation/` | probe | normal-FIML alpha validation | alpha-missing | frozen | keep summaries | Normal-FIML alpha path was validated against external covariance fits. |
| 18 alpha paper simulation | `experiments/studies/18-alpha-paper-simulation/` | paper-facing study | alpha simulation | alpha-missing | frozen | keep summaries | Scored-ordinal alpha simulation is fuller than the paper-local driver. |
| 19 alpha overnight regimes | `experiments/probes/19-alpha-overnight-regimes/` | probe | alpha regime search | alpha-missing | frozen | keep summaries | Overnight alpha regime exploration with local helper modules. |
| 20 alpha NT-FIML stress | `experiments/probes/20-alpha-nt-fiml-stress/` | probe | normal-FIML alpha stress | alpha-missing | frozen | keep summaries | Stress cells probe when normal-theory alpha FIML starts to fail. |
| 21 alpha KP2016 categorical | `experiments/probes/21-alpha-kp2016-categorical/` | probe | categorical alpha DGP | alpha-missing | frozen | keep summaries | Ordered-categorical KP2016-style cells document categorical-alpha limits. |
| 22 quadratic NT-FIML validation | `experiments/probes/22-quadratic-nt-fiml-validation/` | probe | quadratic NT-FIML validation | Paper C | frozen | keep summaries | NT-FIML versus external validation remains useful despite old SE framing. |
| 23 quadratic paper simulation | `experiments/studies/23-quadratic-paper-simulation/` | paper-facing study | quadratic simulation | Paper C | frozen | keep summaries | Main quadratic simulation runner, curated by the paper-local script. |
| 24 alpha equal cocron | `experiments/studies/24-alpha-equal-cocron/` | reusable study | alpha equality example | alpha-missing | current | keep summaries | Paired alpha equality test works on `cocron::knowledge` and MCAR amputations. |
| 25 kappa equal examples | `experiments/studies/25-kappa-equal-examples/` | reusable study | kappa equality examples | Paper C | current | keep summaries | Real-data examples demonstrate dependent and independent equality tests. |
| 26 CRACKLES vector kappa | `experiments/workbench/26-crackles-vector-kappa/` | active workbench | vector kappa | library / future paper | active WIP | keep summaries | Active CRACKLES pilot for component-separable vector kappa. |
