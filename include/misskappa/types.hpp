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
// covariance matrix. `psi` is the per-subject influence-function matrix
// (n x K) for estimators that expose it; empty (0 x 0) otherwise. When
// populated, vcov = (1 / n^2) * psi^T psi up to numerical noise, and the
// caller can stack `psi` columns from independent fits on the same data to
// build a joint vcov across estimators / weight schemes / rater pairs.
struct Estimation {
  RealVec estimates;
  RealMat vcov;
  RealMat psi;
};

// psi_i = J * phi_i  =>  psi = phi * J^T. Centralised so estimators that
// already compute the moment-level IF matrix `phi` (n x M) and the
// delta-method Jacobian J (K x M) get the kappa-level IF matrix
// (n x K) in one place.
inline RealMat build_psi_from_phi(const RealMat& phi, const RealMat& J) {
  return phi * J.transpose();
}

}  // namespace misskappa

#endif  // MISSKAPPA_TYPES_HPP
