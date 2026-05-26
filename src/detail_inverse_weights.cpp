#include "detail_inverse_weights.hpp"

namespace misskappa::detail {

namespace {
constexpr double zero_tol = 1e-9;
}

Result<InverseWeights> compute_inverse_weights(
    const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>& rating_mask,
    int n, int R, Reweighting mode) {
  InverseWeights w;
  w.pi_j_inv.setOnes(R);
  w.pi_jk_inv.setOnes(R, R);

  if (mode == Reweighting::available) return w;

  RealMat mask_d = rating_mask.cast<double>();
  RealVec pi_j = mask_d.colwise().sum().transpose() / n;
  if (pi_j.minCoeff() < zero_tol) {
    return std::unexpected(Error::singular_weight);
  }
  w.pi_j_inv = pi_j.cwiseInverse();

  if (mode == Reweighting::gwet) return w;

  // mode == ipw: pi_{jk} = (1/n) sum_i M_ij M_ik, then invert entrywise.
  RealMat pi_jk = (mask_d.transpose() * mask_d) / static_cast<double>(n);
  for (int j = 0; j < R; ++j) {
    for (int k = 0; k < R; ++k) {
      w.pi_jk_inv(j, k) = (pi_jk(j, k) > zero_tol) ? 1.0 / pi_jk(j, k) : 0.0;
    }
  }
  return w;
}

}  // namespace misskappa::detail
