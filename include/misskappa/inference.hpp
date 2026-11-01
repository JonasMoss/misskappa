#ifndef MISSKAPPA_INFERENCE_HPP
#define MISSKAPPA_INFERENCE_HPP

#include "misskappa/estimate.hpp"

namespace misskappa {

// Exact delete-one evaluations of a complete-data coefficient functional.
// `estimate` has length K and `leave_one_out` is n x K, with row i equal to
// the estimate after deleting subject i. This is deliberately separate from
// Estimation: first-order influence-function inference remains the default,
// while callers can compose these finite-sample summaries with any smooth
// transformation or critical-value rule they choose.
struct DeleteOneEstimation {
  RealVec estimate;
  RealMat leave_one_out;
};

// Complete-data kappa delete-one summaries. In addition to the coefficient
// evaluations, the result exposes the ingredients needed by boundary-aware
// intervals for coefficients of the form 1 - observed / chance:
//
// - item_disagreement(i) is subject i's mean pairwise disagreement;
// - chance is the full-sample chance disagreement for each coefficient;
// - chance_leave_one_out is its n x K delete-one matrix.
//
// `item_disagreement` is empty when the supplied loss is not bounded in
// [0, 1]. The coefficient and chance summaries are still returned.
struct KappaDeleteOneEstimation {
  DeleteOneEstimation coefficient;
  RealVec item_disagreement;
  RealVec chance;
  RealMat chance_leave_one_out;
};

// Raw categorical ratings. The agreement matrix may be any finite symmetric
// C x C matrix with diagonal one; outputs follow the estimator order
// (Conger, Fleiss, Brennan--Prediger). Ratings must be complete and n >= 3.
Result<KappaDeleteOneEstimation> delete_one_kappa_complete(
    IntMatView ratings, RealMatView weights);

// Continuous ratings with a caller-supplied POD loss. Outputs are (Conger,
// Fleiss). Ratings must be finite and n >= 3.
Result<KappaDeleteOneEstimation> delete_one_kappa_complete(
    RealMatView ratings, loss::ContinuousLoss loss);

// Count-format ratings. Outputs are (Fleiss, Brennan--Prediger), using the
// same count weighting as the count moment estimators. Delete-one coefficient and
// chance values are exact even when row totals vary. Boundary item moments
// are exposed only when all row totals are equal and at least two.
Result<KappaDeleteOneEstimation> delete_one_kappa_counts_complete(
    IntMatView counts, RealMatView weights,
    CountWeighting weighting = CountWeighting::fleiss_cuzick);

}  // namespace misskappa

#endif  // MISSKAPPA_INFERENCE_HPP
