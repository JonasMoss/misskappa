#ifndef MISSKAPPA_DETAIL_PSI_HPP
#define MISSKAPPA_DETAIL_PSI_HPP

#include "misskappa/types.hpp"

namespace misskappa::detail {

// psi_i = J * phi_i  =>  psi = phi * J^T. Centralised so estimators that
// already compute the moment-level IF matrix `phi` (n x M) and the
// delta-method Jacobian J (K x M) get the kappa-level IF matrix
// (n x K) in one place.
inline RealMat build_psi_from_phi(const RealMat& phi, const RealMat& J) {
  return phi * J.transpose();
}

}  // namespace misskappa::detail

#endif  // MISSKAPPA_DETAIL_PSI_HPP
