#ifndef MISSKAPPA_SRC_DETAIL_ESTIMATE_RAW_HPP
#define MISSKAPPA_SRC_DETAIL_ESTIMATE_RAW_HPP

#include "detail_inverse_weights.hpp"
#include "misskappa/result.hpp"
#include "misskappa/types.hpp"

namespace misskappa::detail {

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
