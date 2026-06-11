#ifndef MISSKAPPA_DETAIL_PSD_INVERSE_HPP
#define MISSKAPPA_DETAIL_PSD_INVERSE_HPP

#include "misskappa/result.hpp"
#include "misskappa/types.hpp"

#include <Eigen/Eigenvalues>
#include <cmath>
#include <limits>

namespace misskappa::detail {

inline RealMat pseudo_inverse_psd(RealMat m, double rcond) {
  if (m.rows() == 0 || m.cols() == 0) return RealMat::Zero(m.rows(), m.cols());
  m = 0.5 * (m + m.transpose());

  Eigen::SelfAdjointEigenSolver<RealMat> es(m);
  if (es.info() != Eigen::Success) {
    return RealMat::Constant(m.rows(), m.cols(), std::numeric_limits<double>::quiet_NaN());
  }

  const RealVec& evals = es.eigenvalues();
  const RealMat& evecs = es.eigenvectors();
  const double lambda_max = evals.maxCoeff();
  if (!(lambda_max > 0.0) || !std::isfinite(lambda_max)) {
    return RealMat::Zero(m.rows(), m.cols());
  }
  const double rc = (std::isfinite(rcond) && rcond > 0.0) ? rcond : 0.0;
  const double threshold = rc * lambda_max;

  RealVec inv = RealVec::Zero(evals.size());
  for (Eigen::Index k = 0; k < evals.size(); ++k) {
    if (evals(k) > threshold) inv(k) = 1.0 / evals(k);
  }
  return evecs * inv.asDiagonal() * evecs.transpose();
}

// Per-row fraction of the delta-method gradient lying in the truncated null
// space of a (symmetric) information matrix -- the directions pseudo_inverse_psd
// drops at the same rcond. This is a diagnostic, not a gate: when the
// coefficient is an estimable function of the identified margins the saturated
// nuisance being flat does not break the point estimate; a fraction away from
// zero flags that the sample information is rank-deficient along directions the
// coefficient touches (sparse support), making the point estimate
// selection-dependent and the truncated SE optimistic about those directions.
// Shared by the saturated multinomial (cat_fiml) and normal-theory (nt_fiml)
// FIML backends so both report the diagnostic against the same truncation.
inline Result<RealVec> reduced_gradient_null_fraction(
    const RealMat& info_star, const RealMat& jacobian_reduced, double rcond) {
  RealVec frac = RealVec::Zero(jacobian_reduced.rows());
  if (jacobian_reduced.cols() == 0 || jacobian_reduced.rows() == 0) return frac;
  if (info_star.rows() != jacobian_reduced.cols()
      || info_star.cols() != jacobian_reduced.cols()) {
    return misskappa::unexpected(Error::numerical_error);
  }

  Eigen::SelfAdjointEigenSolver<RealMat> es(0.5 * (info_star + info_star.transpose()));
  if (es.info() != Eigen::Success) {
    return misskappa::unexpected(Error::numerical_error);
  }

  const RealVec& evals = es.eigenvalues();
  const RealMat& evecs = es.eigenvectors();
  const double lambda_max = evals.maxCoeff();
  const double rc = (std::isfinite(rcond) && rcond > 0.0) ? rcond : 0.0;
  const double threshold =
      (std::isfinite(lambda_max) && lambda_max > 0.0) ? rc * lambda_max : 0.0;

  for (Eigen::Index row = 0; row < jacobian_reduced.rows(); ++row) {
    const Eigen::VectorXd grad = jacobian_reduced.row(row).transpose();
    const double norm_sq = grad.squaredNorm();
    if (!(norm_sq > 0.0) || !std::isfinite(norm_sq)) continue;
    double null_sq = 0.0;
    for (Eigen::Index k = 0; k < evals.size(); ++k) {
      if (evals(k) > threshold) continue;
      const double projection = evecs.col(k).dot(grad);
      if (!std::isfinite(projection)) {
        return misskappa::unexpected(Error::numerical_error);
      }
      null_sq += projection * projection;
    }
    frac(row) = std::sqrt(std::min(1.0, null_sq / norm_sq));
  }
  return frac;
}

}  // namespace misskappa::detail

#endif  // MISSKAPPA_DETAIL_PSD_INVERSE_HPP
