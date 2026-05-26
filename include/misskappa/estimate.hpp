#ifndef MISSKAPPA_ESTIMATE_HPP
#define MISSKAPPA_ESTIMATE_HPP

#include "misskappa/result.hpp"
#include "misskappa/types.hpp"

namespace misskappa {

struct EmOptions {
  int max_iter = 500;
  double tol = 1e-9;
};

// Categorical raw-rating estimators. `ratings` is n x R; entries are
// non-negative integer category codes or `na_code` (= -1) for missing.
// `weights` is the c x c loss matrix returned by misskappa::loss factories.
// Each estimator returns Fleiss/Conger and Brennan-Prediger by default;
// Cohen falls out of the R=2 case of Fleiss/Conger.
Result<Estimation> estimate_available(IntMatView ratings, RealMatView weights);
Result<Estimation> estimate_ipw      (IntMatView ratings, RealMatView weights);
Result<Estimation> estimate_fiml     (IntMatView ratings, RealMatView weights, EmOptions opts);
Result<Estimation> estimate_gwet     (IntMatView ratings, RealMatView weights);

}  // namespace misskappa

#endif  // MISSKAPPA_ESTIMATE_HPP
