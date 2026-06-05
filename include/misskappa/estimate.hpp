#ifndef MISSKAPPA_ESTIMATE_HPP
#define MISSKAPPA_ESTIMATE_HPP

#include "misskappa/loss.hpp"
#include "misskappa/result.hpp"
#include "misskappa/types.hpp"

#include <cstdint>

namespace misskappa {

struct EmOptions {
  int max_iter = 10000;
  double tol = 1e-8;
  double prune_tol = 1e-9;   // patterns with theta < prune_tol are dropped.
  double start_alpha = 0.1;  // smoothing on initial theta; small positive.
  double info_rcond = 5e-5;  // Louis eigenvalues <= info_rcond * lambda_max are dropped.
};

enum class QuadraticVcov {
  empirical,
  normal,
  elliptical
};

struct QuadraticOptions {
  QuadraticVcov vcov = QuadraticVcov::empirical;
  // Relative Mardia kurtosis beta. beta = 1 gives the normal covariance.
  double relative_kurtosis = 1.0;
};

struct GwiseOptions {
  // If g <= 0, use all raters (g = R).
  int g = 0;
  // Chance-tuple cap. Categorical g-wise estimators enumerate finite category
  // support; continuous g-wise estimators enumerate direct n^g item tuples.
  std::int64_t max_chance_tuples = 5000000;
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
// the R=2 case of Fleiss / Conger.
Result<Estimation> estimate_available(IntMatView ratings, RealMatView weights);
Result<Estimation> estimate_ipw      (IntMatView ratings, RealMatView weights);
Result<Estimation> estimate_fiml     (IntMatView ratings, RealMatView weights, EmOptions opts);
Result<Estimation> estimate_gwet     (IntMatView ratings, RealMatView weights);

// --- Closed-form quadratic estimator (raw, real-valued) --------------------
//
// Moment-based estimator for the quadratically-weighted Cohen / Fleiss /
// Conger / Brennan-Prediger family. Treats categorical ratings as numeric
// scores; missing entries are NaN. Closed form via per-rater means and
// covariances; variance uses the asymptotic covariance of those moments
// with the delta method.
//
// `ratings` is n x R real-valued; `values` is the length-C category-score
// vector used to define the quadratic loss. Returns (Conger, Fleiss,
// Brennan-Prediger).
Result<Estimation> estimate_quadratic(
    RealMatView ratings, const RealVec& values, QuadraticOptions opts = {});

// Counts-format counterpart of estimate_quadratic. `counts` is n x C, R is
// the total number of raters. Returns (Fleiss, Brennan-Prediger).
Result<Estimation> estimate_quadratic_counts(IntMatView counts, const RealVec& values, int R);

// --- Counts-format input ---------------------------------------------------
//
// `counts` is n x C of non-negative integers: counts(i, k) is the number of
// raters who assigned subject i to category k. Row sum r_i is the number of
// raters who rated subject i (need not be the same across subjects). Per-
// rater identity is not preserved, so IPW / Gwet (which need rater-specific
// observation rates) are not meaningful here.
//
// Estimates returned: (Fleiss, Brennan-Prediger). Conger requires identified
// raters and is not in scope for counts. FIML for counts (with the multi-
// variate hypergeometric weights needed when row sums vary) is a future
// addition; see dev/notes/todo.md.
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
// without a finite category count.
Result<Estimation> estimate_available_continuous(RealMatView ratings, loss::ContinuousLoss loss);
Result<Estimation> estimate_ipw_continuous      (RealMatView ratings, loss::ContinuousLoss loss);
Result<Estimation> estimate_gwet_continuous     (RealMatView ratings, loss::ContinuousLoss loss);

// --- Coefficient alpha for scored categorical item batteries ----------------
//
// `ratings` is n x R with categories in [0, C-1] or `na_code`; `values` is the
// length-C score vector used to turn categories into item scores. The available
// estimator uses pairwise covariance moments and is an MCAR baseline. The FIML
// estimator fits the saturated C^R full-response distribution by EM and maps it
// to alpha through the implied covariance matrix, so it is intended for small
// fixed-category batteries.
//
// Estimate returned: (alpha).
Result<Estimation> estimate_alpha_available(IntMatView ratings, const RealVec& values);
Result<Estimation> estimate_alpha_fiml(
    IntMatView ratings, const RealVec& values, EmOptions opts);

// --- Closed rectangular g-wise estimators -----------------------------------
//
// Complete-data estimator for symmetric g-wise disagreement kernels. This is
// the Frechet / Hubert-style multirater-distance family: observed
// disagreement averages a g-argument distance over within-item rater
// combinations; chance disagreement uses Cohen-type (distinct rater
// combinations) and Fleiss-type (all ordered rater tuples) V-statistics.
//
// Missing ratings are intentionally not supported here. Categorical entries
// must be non-negative codes in [0, C-1], where C is carried by `distance`.
// Continuous entries must be finite. Estimates returned: (Cohen, Fleiss).
Result<Estimation> estimate_gwise(
    IntMatView ratings, loss::GwiseCategoricalDistance distance,
    GwiseOptions opts = {});
Result<Estimation> estimate_gwise_continuous(
    RealMatView ratings, loss::GwiseContinuousDistance distance,
    GwiseOptions opts = {});

}  // namespace misskappa

#endif  // MISSKAPPA_ESTIMATE_HPP
