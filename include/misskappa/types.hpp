#ifndef MISSKAPPA_TYPES_HPP
#define MISSKAPPA_TYPES_HPP

#include <Eigen/Core>

namespace misskappa {

// Sentinel for a missing categorical rating. Categories are non-negative
// integers; -1 signals a missing entry. Picked deliberately so a single
// signed integer matrix can hold ratings + missingness without an auxiliary
// mask, mirroring how R-side ratings arrive after NA -> -1 translation.
inline constexpr int na_code = -1;

using IntMat = Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>;
using RealMat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
using RealVec = Eigen::Matrix<double, Eigen::Dynamic, 1>;

using IntMatView = Eigen::Ref<const IntMat>;
using RealMatView = Eigen::Ref<const RealMat>;

// Plain aggregate returned by every estimator entry point. Estimates are
// coefficient values in a stable order; vcov is the corresponding asymptotic
// covariance matrix. `labels` is a packed string of fixed-length names so
// the binding layer can label the output without allocating.
struct Estimation {
  RealVec estimates;
  RealMat vcov;
};

}  // namespace misskappa

#endif  // MISSKAPPA_TYPES_HPP
