#ifndef MISSKAPPA_ESTIMATE_HPP
#define MISSKAPPA_ESTIMATE_HPP

#include "misskappa/result.hpp"
#include "misskappa/types.hpp"

namespace misskappa {

struct EmOptions {
  int max_iter = 10000;
  double tol = 1e-8;
  double prune_tol = 1e-9;   // patterns with theta < prune_tol are dropped.
  double start_alpha = 0.1;  // smoothing on initial theta; small positive.
};

// Categorical raw-rating estimators. `ratings` is n x R; entries are
// non-negative integer category codes in [0, C-1], or `na_code` (= -1) for
// missing. `weights` is the C x C AGREEMENT matrix (1 on diagonal, partial
// agreement off-diagonal), matching irrCAC's identity.weights /
// quadratic.weights convention and what the misskappa::loss factories
// return. The estimator computes its disagreement coefficients on
// L = 1 - weights internally.
//
// Estimates returned: (Conger, Fleiss, Brennan-Prediger). Cohen's kappa is
// the R=2 case of Fleiss / Conger.
Result<Estimation> estimate_available(IntMatView ratings, RealMatView weights);
Result<Estimation> estimate_ipw      (IntMatView ratings, RealMatView weights);
Result<Estimation> estimate_fiml     (IntMatView ratings, RealMatView weights, EmOptions opts);
Result<Estimation> estimate_gwet     (IntMatView ratings, RealMatView weights);

}  // namespace misskappa

#endif  // MISSKAPPA_ESTIMATE_HPP
