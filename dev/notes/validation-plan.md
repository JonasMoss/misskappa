# Validation plan

Two independent oracles cover the port:

## Primary: irrcacsmoke

`dev/irrcacsmoke` is the local numeric oracle fork for complete-data
agreement smoke tests. It is derived from `irrCAC`, but smoke-facing
functions return named numeric vectors (`estimate`, `variance`, `se`,
`pa`, `pe`) and use the asymptotic `Gamma / n` variance convention.

- `r-package/tests/testthat/test-parity-irrcacsmoke.R` skips unless
  `irrcacsmoke` is installed locally.
- Raw complete-data smokes compare Conger, Fleiss, and Brennan-Prediger
  estimates plus diagonal variances.
- Counts-format smokes compare Fleiss and Brennan-Prediger estimates plus
  diagonal variances on `dat.fleiss1971`.
- CI should not install oracle packages from CRAN; checked JSON fixtures can
  be added later if C++ golden tests need the same oracle.

## Secondary: dev/legacy/misskappa

The frozen legacy package is the secondary oracle while the port is in
flight, because it is the version we already trust empirically.

- `r-package/tests/testthat/test-parity-legacy.R` loads both the new and
  the legacy package and asserts numerical agreement on a curated set of
  inputs (matches `tests/fixtures/` where possible).
- The legacy package is unbuilt by default; this test installs it locally
  on first invocation and caches the install. It is skipped when the
  legacy install fails (e.g., on systems without Armadillo).

## Unit-level

- `tests/unit/loss_test.cpp` — symmetry, zero diagonal, closed-form known
  values for each weighting (8 categorical, 5 continuous).
- `tests/unit/available_test.cpp` — R=2 / no-missingness reduces to Cohen
  by hand.
- `tests/unit/ipw_test.cpp` — IPW under uniform MCAR equals available-case
  (property).
- `tests/unit/fiml_test.cpp` — EM converges on a small problem with a
  known maximum likelihood root.
- `tests/unit/inference_test.cpp` — finite-difference Jacobian sanity
  checks on variance entries.

## R-level

- `r-package/tests/testthat/test-dispatcher.R` — kappa(x, method, weight)
  routes correctly.
- `r-package/tests/testthat/test-s3.R` — print / coef / vcov / confint /
  as.data.frame.
- `r-package/tests/testthat/test-sim.R` — `sim$mcar`, `sim$mar`, `sim$jsm`
  produce reasonable output and are reproducible under a seed.
- `r-package/tests/testthat/test-parity-irrcac.R` — agreement vs installed
  irrCAC.

## Tolerances

- Closed-form vs library: `1e-12`.
- Library vs irrCAC: `1e-9`.
- Library vs legacy: `1e-9`.
- EM-based vs library closed form (sample-size dependent): `1e-6`.

Tolerances live in `tests/golden/` test code and are documented per fixture
when a fixture needs a looser bound.
