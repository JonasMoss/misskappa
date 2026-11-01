#ifndef MISSKAPPA_ESTIMATE_HPP
#define MISSKAPPA_ESTIMATE_HPP

// Estimator entry points of the misskappa C++ library. The library is the
// computational engine of the misskappa R package: its interface follows the
// package's needs and may change between releases without a deprecation
// period. Pin a commit if you build on it directly.

#include "misskappa/loss.hpp"
#include "misskappa/result.hpp"
#include "misskappa/types.hpp"

#include <cstdint>
#include <vector>

namespace misskappa {

enum class NormalCovarianceModel {
  saturated,
  variance_components,
};

enum class NormalMeanModel {
  free,
  equal,
};

struct EmOptions {
  int max_iter = 10000;
  double tol = 1e-8;
  double prune_tol = 1e-9;   // patterns with theta < prune_tol are dropped.
  double start_alpha = 0.1;  // smoothing on initial theta; small positive.
  double info_rcond = 5e-5;  // Louis eigenvalues <= info_rcond * lambda_max are dropped.
  // Total Dirichlet pseudo-mass spread uniformly over the full C^R pattern
  // table (delta = flatten / C^R added to each cell's expected count in the
  // M-step). With incomplete co-observation of the saturated joint the MLE is
  // a face of the simplex; any flatten > 0 makes the posterior mode unique
  // (the analytic center of that face as flatten -> 0) at the cost of
  // shrinking the fitted table toward uniform with weight
  // flatten / (n + flatten). 0 = strict ML (legacy behaviour, point on the
  // face determined by the deterministic start).
  double flatten = 0.0;
  // Normal-theory FIML covariance model. `saturated` is the existing
  // unrestricted mean/covariance route. `variance_components` is the
  // fixed-rater compound-symmetry model Sigma = psi * J + e * I with
  // psi >= 0, e > 0, fitted by the constrained EM M-step.
  NormalCovarianceModel covariance = NormalCovarianceModel::saturated;
  // Normal-theory FIML mean model. `free` estimates one mean per rater/item.
  // `equal` estimates a common mean and is currently supported with the
  // variance-components covariance model, giving the exchangeable normal
  // model appropriate when rater identity has been aggregated away.
  NormalMeanModel mean = NormalMeanModel::free;
  // Feature-block size for the variance-components covariance model. 1 is the
  // scalar compound symmetry Sigma = psi J + e I. For vector-valued ratings
  // with d features per rater, block_size = d selects the block compound
  // symmetry Sigma = J_R (x) Psi + I_R (x) E (Psi, E both d x d), with the
  // stacked columns ordered rater-major (rater 0 feature 0..d-1, ...).
  int block_size = 1;
};

struct GwiseOptions {
  // If g <= 0, use all raters (g = R).
  int g = 0;
  // Chance-tuple cap. Categorical g-wise estimators enumerate finite category
  // support; continuous g-wise estimators enumerate direct n^g item tuples.
  std::int64_t max_chance_tuples = 5000000;
};

struct NormalFimlEstimation {
  Estimation fit;
  RealVec mu;
  RealMat sigma;
  int iterations = 0;
  bool converged = false;
};

// --- Categorical raw-rating estimators -------------------------------------
//
// `ratings` is n x R; entries are non-negative integer category codes in
// [0, C-1], or `na_code` (= -1) for missing. `weights` is the C x C
// AGREEMENT matrix (1 on diagonal, partial agreement off-diagonal),
// matching irrCAC's identity.weights / quadratic.weights convention and
// what the misskappa::loss factories return. The estimator computes its
// disagreement coefficients on L = 1 - weights internally.
//
// Estimates returned: (Conger, Fleiss, Brennan-Prediger). Cohen's kappa is
// the R=2 case of Fleiss / Conger. The available-case and Gwet reweightings
// live in misskappa::research (research.hpp). With missing ratings, fixed-rater raw
// coefficients require every rater to be observed and every rater pair to be
// jointly observed by at least one subject; otherwise the saturated pairwise
// functional is not identified and the estimators return Error::not_identified.
Result<Estimation> estimate_ipw      (IntMatView ratings, RealMatView weights);
Result<Estimation> estimate_fiml     (IntMatView ratings, RealMatView weights, EmOptions opts);

// --- Counts-format input ---------------------------------------------------
//
// `counts` is n x C of non-negative integers: counts(i, k) is the number of
// raters who assigned subject i to category k. Row sum r_i is the number of
// raters who rated subject i (need not be the same across subjects). Per-
// rater identity is not preserved, so IPW / Gwet (which need rater-specific
// observation rates) are not meaningful here.
//
// Count-data moment weighting. The public estimator uses fleiss_cuzick.
// unit_weighted (research::estimate_counts) is for parity studies against
// software that averages each subject's row agreement equally.
enum class CountWeighting {
  fleiss_cuzick,  // row disagreement weight r_i - 1; F&C unequal-judges convention.
  unit_weighted,  // row disagreement weight 1; irrCAC distribution/count convention.
};

// Estimates returned: (Fleiss, Brennan-Prediger). Conger requires identified
// raters and is not in scope for counts. `estimate_available_counts` is the
// Fleiss--Cuzick count moment estimator. Rows
// with r_i < 2 contain no within-subject pair information and receive zero
// observed-disagreement weight; at least one row with r_i >= 2 is required.
Result<Estimation> estimate_available_counts(IntMatView counts, RealMatView weights);

// FIML / EM counterpart of estimate_available_counts. Models the full
// composition (length-C count vector summing to `r_total`) under
// exchangeable iid raters, and uses uniformly-random hypergeometric
// subsampling to relate full counts to observed counts. Each subject's
// observed row sum may be < r_total; FIML weights subjects by how much
// their partial counts pin down the full composition. Returned
// coefficients: (Fleiss, Brennan-Prediger).
//
// Assumption: raters are exchangeable iid from a common distribution.
// Counts data cannot identify rater-specific distributions; for the
// non-exchangeable case use rater-identified data and estimate_fiml.
Result<Estimation> estimate_fiml_counts(
    IntMatView counts, RealMatView weights, int r_total, EmOptions opts);

// --- Continuous-rating estimators -------------------------------------------
//
// `ratings` is n x R of real values; entries that are not finite (NaN or
// +/-Inf) are treated as missing. `loss` is a continuous loss kernel from
// the misskappa::loss continuous factories (identity_loss, linear_loss,
// quadratic_loss, etc.); it returns disagreement in [0, 1] for any pair of
// rating values.
//
// Estimates returned: (Conger, Fleiss). No Brennan-Prediger row — chance
// disagreement against a uniform reference distribution is not meaningful
// without a finite category count. With missing ratings, every rater pair must
// be jointly observed by at least one subject.
// The IPW and Gwet reweightings live in misskappa::research.
Result<Estimation> estimate_available_continuous(RealMatView ratings, loss::ContinuousLoss loss);

// Normal-theory FIML estimators for smooth functions of the fitted mean and
// covariance. The default covariance model is saturated; EmOptions may select
// the fixed-rater compound-symmetry covariance model. Inputs are real-valued
// with non-finite entries treated as missing. Returned influence rows
// correspond to rows with at least one observed entry; all-empty rows are
// dropped, matching the R-level convention. Scalar quadratic FIML returns
// (Conger, Fleiss, Brennan-Prediger); pass an empty `values` vector to use the
// continuous uniform BP baseline on the observed score range.
Result<NormalFimlEstimation> estimate_quadratic_normal_fiml(
    RealMatView ratings, const RealVec& values, EmOptions opts);
// `ranges`, when non-empty, fixes the per-feature score width (hi - lo) used by
// the Brennan-Prediger chance baseline instead of estimating it from the
// observed data, the vector analogue of the scalar `values` hook. It must be
// length `features`; an empty vector keeps the observed-range behaviour.
Result<NormalFimlEstimation> estimate_vector_quadratic_normal_fiml(
    RealMatView ratings, int features, RealMatView W, EmOptions opts,
    const RealVec& ranges = RealVec());

// --- Component-separable vector-rating estimators ---------------------------
//
// `ratings` is n x (R * p), with columns ordered rater-major:
//   rater 0 feature 0..p-1, rater 1 feature 0..p-1, ...
// Non-finite entries are missing component ratings. `features` is p.
// `loss` carries diagonal feature weights and a component loss. The estimator
// targets component-separable vector disagreement by summing weighted
// component losses over observed rater-feature overlaps and applying the
// loss transform (identity for Hamming/L1/squared, sqrt for RMS).
//
// Estimates returned: (Conger, Fleiss). Categorical full-pattern FIML and
// non-diagonal feature weights are intentionally out of scope here.
Result<Estimation> estimate_pairwise_vector(
    RealMatView ratings, int features, loss::ComponentSeparableLoss loss);
Result<Estimation> estimate_ipw_vector(
    RealMatView ratings, int features, loss::ComponentSeparableLoss loss);

// --- Closed rectangular g-wise estimators -----------------------------------
//
// Categorical / continuous IPW and categorical FIML estimators for
// symmetric g-wise disagreement kernels. This is the Frechet / Hubert-style
// multirater-distance family: observed disagreement averages a g-argument
// distance over within-item rater combinations; chance disagreement uses
// Cohen-type (distinct rater combinations) and Fleiss-type (all ordered rater
// tuples) V-statistics.
//
// The complete-data-only versions live in misskappa::research; with no missing
// ratings the IPW estimators reduce to them. `estimate_ipw_gwise` accepts categorical ratings with `na_code`, and
// `estimate_ipw_gwise_continuous` accepts non-finite missing ratings, under
// the same MCAR IPW convention as the pairwise IPW estimators.
// `estimate_fiml_gwise` fits the saturated categorical full-response
// distribution by EM under ignorable missingness. Estimates returned:
// (Cohen, Fleiss). With missing ratings, every requested rater g-tuple must be
// jointly observed by at least one subject.
Result<Estimation> estimate_ipw_gwise(
    IntMatView ratings, loss::GwiseCategoricalDistance distance,
    GwiseOptions opts = {});
Result<Estimation> estimate_fiml_gwise(
    IntMatView ratings, loss::GwiseCategoricalDistance distance,
    EmOptions em_opts, GwiseOptions gwise_opts = {});
Result<Estimation> estimate_ipw_gwise_continuous(
    RealMatView ratings, loss::GwiseContinuousDistance distance,
    GwiseOptions opts = {});

}  // namespace misskappa

#endif  // MISSKAPPA_ESTIMATE_HPP
