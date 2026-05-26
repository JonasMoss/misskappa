#ifndef MISSKAPPA_ESTIMATE_HPP
#define MISSKAPPA_ESTIMATE_HPP

#include "misskappa/loss.hpp"
#include "misskappa/result.hpp"
#include "misskappa/types.hpp"

namespace misskappa {

struct EmOptions {
  int max_iter = 10000;
  double tol = 1e-8;
  double prune_tol = 1e-9;   // patterns with theta < prune_tol are dropped.
  double start_alpha = 0.1;  // smoothing on initial theta; small positive.
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
Result<Estimation> estimate_quadratic(RealMatView ratings, const RealVec& values);

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

}  // namespace misskappa

#endif  // MISSKAPPA_ESTIMATE_HPP
