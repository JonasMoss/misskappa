# 06-louis-spectrum

Spectral diagnostic for the FIML Louis observed-information variance.

The runner simulates one DGP-A data set from experiment 01, fits raw FIML,
and decomposes Conger's-kappa variance across eigenvectors of the reduced
Louis observed-information matrix. It writes both the untruncated contribution
and the contribution retained by `em_options$info_rcond`.

Run from this directory after reinstalling the R package:

```sh
Rscript run_diagnostic.R
```

The key number in `results/summary.md` is the share of untruncated variance
coming from eigenvalues with `lambda / lambda_max <= 1e-3`. A large share means
the old pseudo-inverse was dominated by finite-sample lifts of weak or
unidentified directions.
