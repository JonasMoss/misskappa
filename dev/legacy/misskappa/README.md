# misskappa

`misskappa` computes multi-rater agreement coefficients (Conger, Fleiss, Brennan–Prediger) with support for missing data.

## Quick start

Raw categorical ratings (subjects × raters):

```r
library(misskappa)

x <- matrix(c(1, 1, NA,
              2, 2, 2,
              1, 2, 1),
            nrow = 3, byrow = TRUE)

fit <- kappa_raw(x, method = "available", weight = "quadratic")
fit
coef(fit)
vcov(fit)
as.data.frame(fit)

# Or use the unified dispatcher:
kappa(x)
```

## Choosing a method

The `method` argument selects an estimator family:

- `method = "available"`: available-case estimator (uses observed pairs; simplest baseline).
- `method = "ipw"`: inverse-probability weighting adjustment (for settings where missingness depends on observed data via rater/subject patterns).
- `method = "ml"`: maximum likelihood via an EM algorithm (categorical raw/counts).
- `method = "quadratic"`: a separate moment-based “quadratic” code path (legacy/special-case).

Compatibility / comparison:

- `method = "gwet"`: supported for backwards compatibility and comparison against other implementations; not the primary recommended interface.

Notes:

- `kappa_counts()` expects count tables (subjects × categories) and does not retain per-rater missingness patterns, so methods like `"ipw"`/`"gwet"` are not currently available for counts.
- For continuous data, use `kappa_continuous()`; the weighting is implemented via continuous loss functions.

## Return value

All `kappa_*()` functions return an object of class `misskappa_estimate` with:

- `coef(x)` / `x$estimates`: point estimates
- `vcov(x)` / `x$vcov`: asymptotic variance-covariance matrix
- `as.data.frame(x)`: a convenient table (estimate + standard error)
