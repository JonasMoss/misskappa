# Projected Satterthwaite calibration for robust scores

Status: research note. The immediate next step is a focused low-sample-size
probe in the sibling `magmaan` repository. This note records the mathematical
target and the order of work; the active checkboxes remain in
[`todo.md`](todo.md).

## Motivation

The Pearson-correlation interval project starts from an affine estimating
function, studentizes it with its empirical sandwich variance, and uses a
normal-reference Satterthwaite correction for the noise in that variance
estimate. Let \(\psi_i\) be the numerator influence function,
\(V=\operatorname{Var}(\psi_i)\), and \(h_i\) the influence function of
\(\widehat V\). Removing the part of the denominator fluctuation predictable
from the numerator gives

\[
 h_i^\perp
 = h_i-\frac{\operatorname{Cov}(\psi_i,h_i)}{V}\psi_i,
 \qquad
 \nu_{\mathrm{noise}}
 = \frac{2nV^2}{\operatorname{Var}(h_i^\perp)}.
\]

For Pearson correlation under the bivariate-normal working model,
\(\operatorname{Var}(h_i^\perp)=4V^2\), hence
\(\nu_{\mathrm{noise}}=n/2\). The important general object is the projected
denominator influence. The striking, coefficient-specific fact is that its
relative variance collapses to the constant \(4\).

The same construction applies to fully sandwich-studentized robust score tests
in SEM. This is not a proposal to further correct a normal-theory likelihood
ratio or score statistic whose nonnormal limit is a weighted chi-square
mixture. The statistic of interest is already asymptotically pivotal.

## Scalar scores

For a scalar unbiased score or estimating-function contribution \(\psi_i\),

\[
 T_n
 = \frac{\sqrt n\,\overline\psi}{\sqrt{\widehat V}},
 \qquad
 T_n^2 \overset{d}{\longrightarrow}\chi_1^2.
\]

If nuisance estimation does not affect \(\widehat V\) to first order, then
\(h_i=\psi_i^2-V\). Writing

\[
 \beta_2=\frac{E(\psi_i^4)}{V^2},
 \qquad
 \gamma_1=\frac{E(\psi_i^3)}{V^{3/2}},
\]

the projected relative variance is

\[
 \frac{\operatorname{Var}(h_i^\perp)}{V^2}
 =\beta_2-1-\gamma_1^2,
 \qquad
 \nu_{\mathrm{noise}}
 =\frac{2n}{\beta_2-1-\gamma_1^2}.
\]

This is the simplest statement of when the denominator is attractive. A
Gaussian linear score gives \(\nu\approx n\). A Gaussian quadratic score is
heavier tailed and generally gives a smaller fraction of \(n\). Estimated
nuisance quantities add terms to \(h_i\); Pearson correlation is an example
where those terms materially simplify the final answer.

Estimating \(\beta_2\) and \(\gamma_1\) directly is not the intended headline
method. For a quadratic score, their empirical estimates require information
through order eight in the raw observations. The working-model calculation
should reduce those moments to stable constants or traces of matrices built
from fitted second moments. An empirical moment or delete-one calculation is a
useful diagnostic and validation oracle, not necessarily a good low-\(n\)
procedure.

## Multivariate robust scores

Let \(\psi_i\in\mathbb R^q\) be the nuisance-efficient score contributions for
a rank-\(q\) null hypothesis, and let

\[
 B=E(\psi_i\psi_i^{\mathsf T}),\qquad
 Q_n=n\overline\psi^{\mathsf T}\widehat B^{-1}\overline\psi.
\]

When \(\widehat B\) consistently estimates the full score covariance,

\[
 Q_n\overset{d}{\longrightarrow}\chi_q^2
\]

under the null, including under fixed nonnormality. No eigenvalue-mixture
correction is needed for this statistic.

At low \(n\), the definition of \(\widehat B\) matters. The ordinary Hotelling
benchmark uses the centered sample covariance

\[
 \widehat B_c
 =\frac{1}{n-1}\sum_i
   (\psi_i-\overline\psi)(\psi_i-\overline\psi)^{\mathsf T}.
\]

Score sandwiches often use the raw outer product
\(\widehat B_r=n^{-1}\sum_i\psi_i\psi_i^{\mathsf T}\), because both estimate
\(B\) under the null and differ only by \(O_p(n^{-1})\). The difference is not
innocuous for the proposed low-\(n\) correction: \(\widehat B_r\) contains the
rank-one numerator term \(\overline\psi\overline\psi^{\mathsf T}\), so even
Gaussian score contributions do not give the ordinary Hotelling law. The
magmaan probe must audit which meat each score path uses and compare raw versus
centered meat before attributing any improvement to degrees of freedom.
Centering is a finite-sample meat correction, distinct from the reference-law
match.

Whiten the numerator with \(z_i=B^{-1/2}\psi_i\). Let \(H_i\) be the symmetric
matrix influence function of the normalized meat
\(B^{-1/2}\widehat B B^{-1/2}\). Project each element of \(H_i\) off the whole
numerator vector:

\[
 (H_i^\perp)_{ab}
 = (H_i)_{ab}
 -\operatorname{Cov}\{(H_i)_{ab},z_i\}^{\mathsf T}z_i.
\]

If the remaining meat fluctuation is approximated by a normalized Wishart
matrix, a natural Frobenius moment match is

\[
 \eta_{\mathrm{eff}}
 = \frac{nq(q+1)}{E\lVert H_i^\perp\rVert_F^2}.
\]

The corresponding approximate Hotelling reference is

\[
 \frac{\eta_{\mathrm{eff}}-q+1}
      {\eta_{\mathrm{eff}}q}Q_n
 \ \dot\sim\
 F_{q,\eta_{\mathrm{eff}}-q+1}.
\]

For \(q=1\), this reduces to
\(Q_n\dot\sim F_{1,\eta_{\mathrm{eff}}}
=t_{\eta_{\mathrm{eff}}}^2\) and reproduces the scalar projected
Satterthwaite formula. If the whitened score contributions themselves are
Gaussian and \(\widehat B=\widehat B_c\), their sample covariance is Wishart
and independent of their mean; the formula gives
\(\eta_{\mathrm{eff}}\approx n\), with the usual \(n-1\) finite-sample
covariance degrees of freedom.

For Gaussian covariance-structure models, the casewise scores are centered
quadratic forms rather than Gaussian vectors. Consequently:

- the numerator and empirical meat are generally dependent;
- the projection is not optional;
- \(E\lVert H_i^\perp\rVert_F^2\) involves Gaussian moments through order
  eight, but Isserlis reduction should express them as matrix traces;
- one universal fraction of \(n\) is unlikely, but a cheap
  model- and restriction-specific trace formula may still exist.

A single \(\eta_{\mathrm{eff}}\) cannot reproduce an arbitrary covariance
tensor for the elements of \(\widehat B\). The Hotelling approximation is most
credible when the projected meat fluctuation is close to Wishart in shape.
The probe should therefore report diagnostics of that approximation rather
than only its rejection rate.

## What counts as a useful simplification

The general projection always exists. The method retains the Pearson paper's
simplicity and performance advantages only when most of the following hold:

1. The tested score is scalar or low-dimensional.
2. The robust meat is estimable and well-conditioned at the target sample
   size.
3. Working-model moment reduction yields a constant or a small trace
   calculation, not an empirically estimated eighth-moment tensor.
4. Nuisance estimation can be included analytically without delete-one refits.
5. The effective degrees of freedom are constant over a confidence-set
   inversion, or the task is a test requiring only one null fit.
6. Bias correction of the meat is treated separately from its reference
   degrees of freedom.
7. Boundary and one-sided-tail behavior do not dominate the symmetric
   calibration.

The condition \(\eta_{\mathrm{eff}}>q-1\) is also substantive. Failure of this
condition, singular \(\widehat B\), or extreme meat conditioning should be
reported as evidence that the sample contains too little information for the
fully studentized rank-\(q\) test, not hidden by a fallback to a sharper
reference law.

## Coefficient priorities

The coefficient applications remain useful because they identify cases where
both numerator and denominator simplify.

### Partial correlation and one-predictor partial \(R^2\)

After residualizing two variables on \(k\) controls, the Pearson affine
estimating function applies to the residual pair. Under joint normality, the
residual cross-product matrix has the same Wishart geometry as an ordinary
correlation with effective sample size \(n-k\). The first conjecture to prove
is therefore

\[
 \nu_{\mathrm{partial}\ r}\approx\frac{n-k}{2}.
\]

This would preserve affine inversion and a constant critical value. A
single-predictor increment in \(R^2\) is the squared partial correlation and
comes with the same result.

### Coefficient alpha

For \(P\) complete items, coefficient alpha has the affine contribution

\[
 U_i(\alpha)
 =C_i-\left(1-\frac{\alpha}{c_P}\right)T_i,
 \qquad c_P=\frac{P}{P-1}.
\]

Under a Gaussian tau-equivalent working model, the whitened quadratic-form
matrix has eigenvalues proportional to

\[
 (\underbrace{1,\ldots,1}_{P-1},-(P-1)).
\]

A preliminary calculation, still requiring an independent symbolic audit,
gives

\[
 \frac{\operatorname{Var}(h_i^\perp)}{V^2}
 =6+\frac{4}{P(P-1)},
 \qquad
 \nu_\alpha
 =\frac{n}{3+2/\{P(P-1)\}}.
\]

Thus the reference approaches \(t_{n/3}\), independently of alpha, while a
general congeneric covariance should give a fitted trace formula. This is a
promising explanation for poor small-sample performance of affine alpha
inversion with a \(t_{n-1}\) cutoff. It cannot by itself repair an exact or
near-perfect-agreement boundary.

### Regression, ICC, and broader covariance ratios

For a scalar fixed-design regression contrast, the working-normal robust
variance is a quadratic residual form and gives the established trace rule

\[
 \nu=\frac{\{\operatorname{tr}(A)\}^2}{\operatorname{tr}(A^2)}.
\]

Bell--McCaffrey, Fay--Graubard, and the CR2/approximate-Hotelling literature
already occupy this application. It is a validation lineage rather than the
main novelty. Balanced Gaussian ICC and other ratios of covariance components
may again yield low-rank trace collapses, but exact normal-model \(F\)
procedures and existing small-sample corrections make them lower priority than
partial correlation and alpha.

## SEM target

The proposed SEM statistic is the fully sandwich-studentized nuisance-efficient
score. The same construction covers:

- a one-degree-of-freedom robust modification index;
- a joint score test for releasing several equality or fixed-parameter
  constraints;
- a nested-model score test requiring only the restricted fit;
- global exact fit, with the tested directions spanning the saturated
  residual-moment complement of the model tangent.

For a one-dimensional test, the proposed change is simply to replace the
ordinary \(\chi_1^2\) reference by \(F_{1,\eta_{\mathrm{eff}}}\). For
\(q>1\), use the approximate Hotelling reference above.

Score-based RMSEA, CFI, or other descriptive fit indices are not immediate
consequences. Pivotal null inference supplies an exact-fit or nested test.
Population misfit and noncentrality require a separate definition under fixed
alternatives, and baseline-model indices require two such definitions. Keep
that work deferred until the null test is demonstrably calibrated.

### Distinction from existing SEM corrections

Satorra--Bentler, mean-and-variance-adjusted, pEBA, and FMG procedures begin
with a normal-theory statistic whose nonnormal limit is a weighted sum of
\(\chi_1^2\) variables. They correct the mean, variance, or spectrum of that
nonpivotal limit.

The proposed statistic instead uses the full robust score covariance in its
quadratic form and has a \(\chi_q^2\) limit. Its Hotelling correction addresses
finite-sample uncertainty in estimating that covariance. MV-adjusted and FMG
statistics may be useful familiar comparators, but they are not derivations or
special cases of the proposed correction.

## Focused magmaan probe

The first probe should answer one question:

> Does projected working-normal Hotelling calibration improve the low-\(n\)
> null rejection rate of a fully sandwich-studentized SEM score without an
> unacceptable loss of availability or numerical stability?

Keep the initial design deliberately small.

### Model and data

- One correctly specified, complete-data, continuous one-factor model with
  five indicators, interior loadings, and comfortably positive residual
  variances.
- A normal control and one nonnormal family: Vale--Maurelli/Fleishman (VM),
  calibrated to the same population covariance.
- Begin with the moderate VM target of marginal skewness \(2\) and excess
  kurtosis \(7\). Add the severe \(3/21\) cell only if the moderate design is
  informative and numerically stable.
- Sample sizes \(n=30,50,100,200\). The smallest cells are diagnostics, not a
  promise that every rank-\(q\) meat is invertible.
- At least two true nulls: one scalar omitted residual covariance
  (\(q=1\)), and the model's global residual-score test
  (\(q\) equal to the model degrees of freedom). A small-rank joint nested
  constraint can be added after those two paths agree with their algebraic
  identities.

### Primary methods

1. Existing pivotal robust score with its current meat, referred to
   \(\chi_q^2\).
2. Centered-meat pivotal robust score, also referred to \(\chi_q^2\), as the
   finite-sample covariance baseline.
3. Projected working-normal Satterthwaite/Hotelling score using centered meat.
4. An unprojected working-normal Hotelling variant, solely to show whether the
   numerator--meat projection matters.
5. An empirical score-moment or delete-one effective-df diagnostic, not a
   recommended method.

Mean-and-variance-adjusted and FMG/pEBA statistics may be included as secondary
context once they can be generated from the same restricted fit and tested
directions. Do not let that comparison turn the probe into a method zoo.
Score-flip calibration is also deferred for the first pass; existing magmaan
experiments already document its cost and conditioning behavior.

### Mathematical and implementation checks

- Derive the influence function \(H_i\) of the robust score meat, including
  the effect of fitting the nuisance parameters under the null.
- Audit raw versus centered score meat and verify their exact algebraic
  relationship for \(q=1\).
- Verify that the scalar code reduces to
  \(2nV^2/\operatorname{Var}(h_i^\perp)\).
- Verify that Gaussian score contributions reproduce ordinary Hotelling
  scaling.
- Express the Gaussian covariance-score moments as deterministic matrix
  traces. Compare the analytic result with a delete-one or high-precision
  parametric calculation on a few fixed datasets.
- Compute both oracle-population and fitted-null
  \(\eta_{\mathrm{eff}}\) during development. Only the fitted-null version is
  eligible as a method.
- Keep meat bias correction separate. First determine whether reference-law
  calibration helps; add a BRL/HC-style correction only if observed meat bias
  remains a visible independent defect.

### Outcomes

The primary outcome is unconditional type-I rejection at \(5\%\), reported
cellwise with Monte Carlo uncertainty. Also retain:

- the distribution of \(\eta_{\mathrm{eff}}\) and
  \(\eta_{\mathrm{eff}}/q\);
- rejection by lower versus upper score direction for \(q=1\), where
  meaningful;
- singular-meat and \(\eta_{\mathrm{eff}}\le q-1\) rates;
- meat condition number and the distance of its projected fluctuation from the
  matched Wishart covariance shape;
- fit/test failures and runtime.

Do not add power until the normal and VM null cells show whether the correction
is calibrated. If the proposed method fails, the diagnostics should distinguish
meat bias, poor Wishart shape, numerator-linked fluctuation, and high-rank
conditioning.

## Existing local evidence to reuse

The sibling magmaan repository already contains four relevant experiment
lines:

- experiment 22 validates robust scalar modification indices and score tests
  against robust Wald and scaled-LRT behavior;
- experiment 44 finds that an empirical Kauermann--Carroll effective-df
  correction helps alpha and correlation intervals, while also showing that
  empirical fourth-moment degrees of freedom remain unstable under heavy
  tails at very small \(n\);
- experiment 65 shows that empirical whitening of global residual scores can
  become singular or badly conditioned when the model degrees of freedom
  approach the sample rank;
- experiment 68 provides the existing score-test comparator battery and VM
  designs, although its nonpivotal spectrum corrections answer a different
  asymptotic question from the fully sandwich-studentized score proposed here.

Reuse their model fixtures, VM calibration, score extraction, result columns,
and diagnostics where this does not create an experiment-to-experiment source
dependency. Shared mechanics belong in magmaan core or
`experiments/_support`; the new experiment must not source another experiment.

## Priorities

1. Audit the multivariate projection and Wishart moment match, including the
   exact scaling of the approximate Hotelling reference.
2. Audit raw versus centered score meat, then derive the fitted-nuisance
   contribution to the SEM meat influence and
   identify the smallest magmaan surface needed to expose it.
3. Create the focused normal-plus-VM low-\(n\) probe, starting with \(q=1\)
   and the five-indicator global score.
4. Decide from the probe whether the useful formula is a universal fraction,
   a fitted trace calculation, or too unstable to recommend.
5. Only after the SEM null result, return to the partial-correlation and alpha
   coefficient corollaries and consider score-based descriptive fit indices.
