# Direct-on-counts exchangeable CS quadratic kappa

Status: **folded into the package** (2026-06-15). `r-package/R/kappa_counts_cs.R`
backs `kappa_counts(estimator = "nt_fiml")` directly on the counts; the
realizing-table detour (`.counts_exchangeable_scores()`) is removed. Probe 41
also prototyped `pairwise` and `listwise`, but neither shipped (pairwise
redundant with `fleiss_cuzick`; listwise is caller-side row dropping). Derivation
+ validation prototype: `experiments/workbench/41-counts-exchangeable-cs-direct`.

## Why

`kappa_counts(estimator = "nt_fiml", weight = "quadratic")` currently fits the
exchangeable common-mean compound-symmetry (CS) normal model
\(\Sigma = \psi J + e I\) by **expanding the counts into a realizing
subjects-by-`r_total` table** (`.counts_exchangeable_scores()`, a rotating
NA-padded placement) and then calling the raw-data quadratic normal FIML. That
works because under exchangeability the column placement is likelihood-invariant,
but it is a detour: it invents `r_total` columns, an EM over a padded table, and
the rotation hack only exists to keep columns from being all-NA.

The CS model has enough structure to fit **directly on the counts**, with no
table, no `r_total`, no EM.

## Sufficient statistics

A subject contributes only its observed sub-vector of size \(m_i\). Any
sub-vector of an exchangeable CS Gaussian is again CS Gaussian with the *same*
\((\mu,\psi,e)\), so the per-subject likelihood depends on the data only through

\[
m_i=\textstyle\sum_c n_{ic},\quad s_i=\sum_c n_{ic} v_c,\quad q_i=\sum_c n_{ic} v_c^2,
\]

read straight off the counts row \(n_{i\cdot}\) and the scores \(v\).

## Likelihood (closed form)

With \(\Sigma_m=\psi J_m+eI_m\): eigenvalues \(e+m\psi\) (once) and \(e\)
(multiplicity \(m-1\)), so

\[
\log|\Sigma_m| = \log(e+m\psi)+(m-1)\log e,\qquad
\Sigma_m^{-1}=\tfrac1e\!\left(I_m-\tfrac{\psi}{e+m\psi}J_m\right).
\]

Per subject, with \(A_i=q_i-2\mu s_i+m_i\mu^2\) and \(B_i=s_i-m_i\mu\),

\[
\ell_i=-\tfrac12\big[m_i\log 2\pi+\log(e+m_i\psi)+(m_i-1)\log e\big]
       -\tfrac1{2e}\Big(A_i-\tfrac{\psi}{e+m_i\psi}B_i^2\Big).
\]

Maximize \(\sum_i\ell_i\) over \((\mu,\psi,e)\) — a 3-parameter problem.
Profiling \(\mu\) gives the GLS mean
\(\hat\mu(\psi,e)=\frac{\sum_i s_i/(e+m_i\psi)}{\sum_i m_i/(e+m_i\psi)}\).
The minimal sufficient statistic is, per distinct row-total \(m\):
\((N_m,\ \sum s_i,\ \sum s_i^2,\ \sum q_i)\).

## Kappas (functions of \((\psi,e)\), `r_total`-free)

Observed quadratic disagreement of two raters of one subject is
\(D_o=\operatorname{Var}(Y_1-Y_2)=2e\). Hence

\[
\kappa_{\text{Fleiss}}=\frac{\psi}{\psi+e}\ (=\text{ICC}),\qquad
\kappa_{\text{BP}}=1-\frac{2e}{c_{\text{BP}}},\quad
c_{\text{BP}}=\tfrac{2}{C^2}\big(C\textstyle\sum v_c^2-(\sum v_c)^2\big).
\]

These match `.kqf_grad()` evaluated at the equal-mean CS point (\(t_3=0\),
\(t_2=R(\psi+e)\), \(t_1=R^2\psi+Re\)), confirming the existing coefficient code
and the `r_total`-invariance recorded in
`project_exchangeable_cat_fiml_exposure`.

## The estimators, one sufficient statistic

All estimate \((\mu,\psi,e)\) from the same per-subject stats and map to the
kappas; all carry robust (empirical sandwich) influence rows obeying the package
convention \(\mathrm{vcov}=n^{-2}\Psi^\top\Psi\).

**Shipped (`kappa_counts`):**

- **fiml** (`estimator = "nt_fiml"`) — maximize \(\sum_i\ell_i\). Casewise score
  \(g_i=\nabla_\theta\ell_i\) and total observed information by central
  differences; IF \(=\bar{\mathcal I}^{-1}g_i\), delta-method to \(\kappa\). The
  efficient, MAR-valid member.
- The **moment / available-case** member is the existing **`fleiss_cuzick`**
  (Fleiss-Cuzick 1979) — see below; it fills the "pairwise" slot, so no separate
  counts `pairwise` estimator is shipped.

Counts surface: `{fleiss_cuzick, cat_fiml, nt_fiml}`. **Listwise is not an
estimator option** — complete-case is `na.omit()`-style row dropping the caller
does before calling any estimator, not something the package bakes in.

**Documented but NOT shipped (redundant with `fleiss_cuzick`):**

- **pairwise** — the **weighted-pairwise (WP) intraclass-correlation estimator**
  (Pearson 1896; Karlin-Cameron-Williams 1981 for unbalanced sizes — *not*
  "available-case covariance", which is the saturated pairwise-deletion matrix;
  see `dev/refs/icc-stuff/refs.md`):
  \(\hat e=\frac{\sum_i(m_iq_i-s_i^2)}{2\sum_i\binom{m_i}{2}}\) from within-pair
  squared differences (\(\mu\)-free; this identity needs equal means — the
  unequal-mean generalisation is the GEE exchangeable-moment estimator and does
  not collapse to counts), \(\widehat{\psi+e}=\frac{\sum q_i}{\sum m_i}-\hat\mu^2\)
  from all tokens, \(\hat\psi=\widehat{\psi+e}-\hat e\). This is a *second*
  moment estimator in the same slot as `fleiss_cuzick`, differing only in
  subject weighting (pooled-pairs \(m_i(m_i-1)\) vs Fleiss-Cuzick's \(m_i-1\))
  and chance construction; identical on balanced complete data. Dropped from the
  package as redundant (probe 41 retains the prototype); if the pooled-pairs
  weighting is ever wanted, expose it as a weighting option under the
  Fleiss-Cuzick machinery (see `dev/notes/counts-weighting-fc-gls.md`).

## Verified (probe 41)

- **Exact**: direct fiml reproduces the table-expansion `nt_fiml` to
  \(\le 2.6\times10^{-7}\) (point *and* SE) across 30 random partial-count
  datasets; identical to 6 d.p. on `dat.fleiss1971`, where fiml = pairwise =
  listwise = `fleiss_cuzick` = `cat_fiml` (complete balanced data).
- **Robust SE**: se_ratio \(\approx1.0\), coverage \(\approx0.95\) for every
  estimator in its consistent regime.
- **Prototypes** (true-CS DGP, n=200): under MCAR fiml, the pairwise prototype,
  and the listwise prototype are all consistent (fiml ≈ pairwise efficient;
  listwise less efficient, slight small-\(n\) bias); under the benign anchor-MAR,
  fiml and pairwise stay consistent while **listwise is biased** (Fleiss bias
  \(-0.094\), coverage 0.74) — the textbook "listwise needs MCAR". Pairwise's MAR
  survival is mechanism-dependent (benign MAR); only fiml is the principled
  ignorable-missingness MLE, consistent with the fiml-paper thesis
  (`project_irrcac_se_check`).

## Landed / follow-up

Done: `kappa_counts(estimator = "nt_fiml")` now fits directly (no table); the
table helper is gone; the moment slot stays `"fleiss_cuzick"`. Shipped counts
surface: `{fleiss_cuzick, cat_fiml, nt_fiml}`. The standalone counts `"pairwise"`
was prototyped (probe 41) but **dropped as redundant** with `fleiss_cuzick` (same
moment slot); `"listwise"` was likewise dropped — complete-case is caller-side
`na.omit()`-style row dropping, not a baked-in estimator. testthat parity +
influence-identity tests added; `kappa_counts.Rd` regenerated.

The estimator story (each FIML paired with its moment foil):

| covariance | mean | FIML | moment foil |
|---|---|---|---|
| saturated | free \(\mu_j\) | raw `nt_fiml` (default) | raw `pairwise` (available-case cov.) |
| CS | equal \(\mu\) | counts `nt_fiml` | counts `fleiss_cuzick` |
| CS | free \(\mu_j\) | raw `nt_fiml_cs` | **— (open; would be `pairwise_cs`)** |

Follow-up (open): the bottom-right cell — the moment foil for the **free-mean CS**
(Conger) FIML — is the rater-centered pairwise-CS / GEE-exchangeable moment
estimator. Parked; needs rater identity (raw-data only), a `t3` mean-spread bias
correction, and is MCAR-only. Write-up in
`papers/00_wip/quadratic/dev/notes/cs-fixed-rater-pairwise.md`; see also
dev/notes/todo.md.
