#ifndef MISSKAPPA_SRC_DETAIL_ESTIMATE_RAW_HPP
#define MISSKAPPA_SRC_DETAIL_ESTIMATE_RAW_HPP

#include "misskappa/result.hpp"
#include "misskappa/types.hpp"

namespace misskappa::detail {

// Pair-reweighting mode for the raw-ratings estimator family. Selects how
// `compute_inverse_weights` constructs the per-rater and per-pair inverse
// probabilities pi_j^{-1} and pi_{jk}^{-1}.
enum class Reweighting {
  available,  // no reweighting; pi_j^{-1} = 1, pi_{jk}^{-1} = 1.
  ipw,        // pi_j^{-1} = n / sum_i M_ij, pi_{jk}^{-1} = n / sum_i M_ij M_ik.
  gwet,       // pi_j^{-1} = n / sum_i M_ij, pi_{jk}^{-1} = 1.
};

// All three reweighting modes share the same kernel / influence-function
// pipeline. Public entry points in estimate.hpp call this with their mode.
//
// ratings: n x R; non-negative integer category codes in [0, C-1], or
//          na_code (= -1) for missing entries.
// weights: C x C agreement weight matrix from misskappa::loss factories.
// Returns Estimation{ estimates = (Conger, Fleiss, Brennan-Prediger), vcov }.
Result<Estimation> estimate_raw(
    IntMatView ratings, RealMatView weights, Reweighting mode);

}  // namespace misskappa::detail

#endif  // MISSKAPPA_SRC_DETAIL_ESTIMATE_RAW_HPP
