# Quadratic Estimator Story

Notes on the specialized quadratic path (`estimate_quadratic`) and how to
describe it without overclaiming.

## What It Is

The quadratic estimator is not the IPW estimator with quadratic weights. It is
a closed-form moment estimator for the quadratically weighted Conger, Fleiss,
and Brennan-Prediger coefficients. It treats ratings as numeric scores and
estimates the saturated first and second moments:

- per-rater means;
- pairwise rater covariances;
- the asymptotic covariance of the derived summaries using third- and
  fourth-order moments;
- the final kappa covariance by the delta method.

This is close to the ICC / covariance-structure / limited-information SEM
story: quadratic kappa is a smooth function of the saturated mean-covariance
model.

## Relation To Available-Case And IPW

For complete data, the quadratic moment estimator and the usual U-statistic
available-case estimator with quadratic weights coincide numerically. Since
IPW weights are all one in complete data, IPW also coincides.

With missing data, the estimators are different finite-sample procedures:

- IPW is a Hajek-style pairwise loss-ratio estimator with inverse observation
  weights.
- The quadratic estimator is a plug-in moment estimator based on pairwise
  available means and covariances.

Under MCAR / PMCAR conditions that identify the same first and second moments,
they target the same quadratic coefficient asymptotically. They should not be
described as first-order equivalent in general, because the estimating
equations and influence functions differ once missingness enters.

## Efficiency Claim To Make

Safe version:

> For quadratic loss, the coefficient is a smooth function of saturated first
> and second moments. Under MCAR / PMCAR, the closed-form quadratic estimator is
> the delta-method plug-in based on pairwise-available moment equations, with a
> sandwich/Godambe covariance from the pairwise moment vector.

Possible limited-information version, if we want the SEM connection:

> It is efficient for the induced limited-information experiment that observes
> only the pairwise moment information.

This is a limited-information / Godambe claim, not a full observed-data
semiparametric efficiency claim.

## Claims To Avoid

Do not claim MAR efficiency for the current quadratic path. Pairwise-available
moments are generally not MAR consistent: whether a rater pair is observed can
depend on other observed ratings that are correlated with the pair. SEM's clean
MAR story belongs to observed-likelihood / FIML or to augmented estimating
equations, not to unaugmented pairwise deletion.

Also avoid saying simply that the quadratic estimator is "IPW but faster." It
has the same complete-data endpoint and the same MCAR estimand, but it is a
different estimator.

## Manuscript Framing

Good framing:

- complete data: exactly the ICC / covariance-moment special case;
- MCAR pairwise data: consistent, computationally clean, robust
  delta-method inference;
- pairwise-moment class: plausible limited-information / Godambe efficiency;
- full observed-data MCAR / MAR: not generally efficient, and MAR needs FIML or
  an augmented estimating-equation analogue.

The strongest efficiency statement should stay restricted to the special case
already in the manuscript: pairwise-complete, two-rater exchangeable data, where
there are no additional observed variables outside the pair that could improve
the estimator.
