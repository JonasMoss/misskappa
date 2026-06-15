# Probe 43 — vector-valued block-compound-symmetry (block-CS) covariance

**Question.** Add the compound-symmetry covariance to the *vector-valued*
quadratic FIML: each rater records a d-feature profile, the stacked covariance is
Rd x Rd, and the structured model is block compound symmetry
`Sigma = J_R (x) Psi + I_R (x) E` (Psi, E both d x d). Lock the math before the
C++ surgery.

## Claims validated (`verify.R`)

1. **Projection = constrained MLE.** For block-CS the constrained Gaussian MLE
   given the unconstrained covariance S is *block averaging*:
   `Psi = mean off-diagonal d x d block of S`, `E = mean diagonal block - Psi`.
   Same fixed rater-eigenbasis argument as the scalar case: the orthonormal
   rater transform (Helmert (x) I_d) simultaneously block-diagonalises every
   candidate Sigma into one `R*Psi + E` block (mean direction) and (R-1) copies
   of `E` (contrast directions), so the likelihood factorises and the MLE is the
   eigen-block average = block average. Numerically: block average vs a BFGS
   maximiser of the Gaussian log-likelihood agree to **~1e-9**, negll gap
   **9e-16**. (At d=1 this reduces exactly to `project_variance_components`.)

2. **Coefficient consistency.** Conger/Fleiss/BP from `(mu, block-CS Sigma)` via
   the `J / I / (I - J/R) (x) W` contractions match the population coefficients
   to **0.0018 at n=8k** (sampling error; consistent).

3. **d=1 reduction.** The vector contractions reduce *exactly* to the scalar
   `t1 = 1'Sigma1`, `t2 = trSigma`, `t3 = sum(mu_j - mubar)^2` forms (diff 0).

## Vector Brennan-Prediger chance

Under the independent-uniform chance model the off-diagonal W terms drop out, so
`D_e^BP = sum_a W_aa * (range_a^2 / 6)` (per-feature continuous-uniform baseline,
W-weighted) — the vector generalisation of the scalar continuous BP baseline.

## Shipped

C++ (`src/estimate_normal_fiml.cpp`): `project_block_variance_components` (block
averaging + PSD guards on E and R*Psi+E), `variance_component_jacobian` extended
with a block-size argument (theta = (mu, vech Psi, vech E)), `vector_quadratic_gradient`
extended to the 3-vector incl. BP, and `estimate_vector_quadratic_normal_fiml`
relaxed to accept the variance-components covariance (it sets `block_size =
features`). The EM driver is generic over the stacked column count, so the only
covariance-specific touch points are the M-step projection and the info-reduction
Jacobian. `EmOptions.block_size` defaults to 1 (scalar, byte-identical to the
prior behaviour).

R: `kappa_vector_quadratic(method = "nt_fiml_cs")` (internal). Parity gates:
vector d=1 == scalar `nt_fiml_cs` to **3e-16** (point) / **1e-18** (vcov);
fitted Sigma is exactly block compound symmetric; IF contract
`vcov == crossprod(psi)/n^2` exact. C++ unit tests + testthat added.
