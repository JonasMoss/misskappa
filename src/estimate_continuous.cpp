// Continuous-ratings estimators (available-case, IPW, Gwet).
//
// Algorithm parallels src/estimate_raw.cpp: same available / IPW / Gwet
// reweighting structure, same V/U-statistic kernel pipeline, same delta-
// method variance. The differences are:
//   - inputs are real-valued, with NaN (or +/-Inf) marking missing entries
//     instead of an integer sentinel;
//   - the loss is evaluated by calling a continuous loss kernel rather than
//     indexing a finite C x C agreement matrix;
//   - there are only two coefficients in the output (Conger, Fleiss); the
//     Brennan-Prediger chance baseline is not meaningful without a finite
//     category count.

#include "detail_inverse_weights.hpp"
#include "detail_kernel_moments.hpp"
#include "misskappa/estimate.hpp"

#include <cmath>
#include <limits>

namespace misskappa {

namespace {

constexpr double zero_tol = 1e-9;

Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>
build_finite_mask(RealMatView ratings) {
  Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> m(ratings.rows(), ratings.cols());
  for (Eigen::Index i = 0; i < ratings.rows(); ++i) {
    for (Eigen::Index j = 0; j < ratings.cols(); ++j) {
      m(i, j) = std::isfinite(ratings(i, j)) ? 1 : 0;
    }
  }
  return m;
}

Result<Estimation> estimate_continuous(
    RealMatView ratings, loss::ContinuousLoss loss, detail::Reweighting mode) {
  const int n = static_cast<int>(ratings.rows());
  const int R = static_cast<int>(ratings.cols());
  if (n < 1) return std::unexpected(Error::invalid_argument);
  if (R < 2) return std::unexpected(Error::invalid_argument);

  const auto mask = build_finite_mask(ratings);
  bool any_observed = mask.sum() > 0;
  if (!any_observed) return std::unexpected(Error::invalid_argument);

  auto wres = detail::compute_inverse_weights(mask, n, R, mode);
  if (!wres) return std::unexpected(wres.error());
  const RealVec& pi_j_inv = wres->pi_j_inv;
  const RealMat& pi_jk_inv = wres->pi_jk_inv;

  // Helper to evaluate the loss kernel.
  auto L = [&](double a, double b) {
    return loss.compute(a, b, loss.min_val, loss.max_val);
  };

  // --- Per-subject U-statistic for observed disagreement ---
  RealVec h_dN = RealVec::Zero(n);
  RealVec h_dD = RealVec::Zero(n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < R - 1; ++j) {
      if (!mask(i, j)) continue;
      for (int k = j + 1; k < R; ++k) {
        if (!mask(i, k)) continue;
        const double w_jk = pi_jk_inv(j, k);
        h_dN(i) += L(ratings(i, j), ratings(i, k)) * w_jk;
        h_dD(i) += w_jk;
      }
    }
  }
  const double psi_dN_hat = h_dN.mean();
  const double psi_dD_hat = h_dD.mean();

  // --- V-statistic kernels (chance disagreement under independence) ---
  detail::KernelMoments kernel_CN(n);
  detail::KernelMoments kernel_CD(n);
  detail::KernelMoments kernel_FN(n);
  detail::KernelMoments kernel_FD(n);
  for (int i = 0; i < n; ++i) {
    for (int ip = 0; ip < n; ++ip) {
      double h_cn = 0, h_cd = 0, h_fn = 0, h_fd = 0;
      for (int j = 0; j < R; ++j) {
        if (!mask(i, j)) continue;
        const double w_j = pi_j_inv(j);
        const double a = ratings(i, j);
        for (int k = 0; k < R; ++k) {
          if (!mask(ip, k)) continue;
          const double w_jk = w_j * pi_j_inv(k);
          const double l = L(a, ratings(ip, k));
          // Fleiss: all (j, k).
          h_fn += l * w_jk;
          h_fd += w_jk;
          if (j < k) {  // Conger: j < k only.
            h_cn += l * w_jk;
            h_cd += w_jk;
          }
        }
      }
      kernel_CN.add(i, ip, h_cn);
      kernel_CD.add(i, ip, h_cd);
      kernel_FN.add(i, ip, h_fn);
      kernel_FD.add(i, ip, h_fd);
    }
  }
  const double psi_CN_hat = kernel_CN.mean(n);
  const double psi_CD_hat = kernel_CD.mean(n);
  const double psi_FN_hat = kernel_FN.mean(n);
  const double psi_FD_hat = kernel_FD.mean(n);

  // --- Point estimates (Conger, Fleiss) ---
  const double d_hat   = (psi_dD_hat > zero_tol) ? psi_dN_hat / psi_dD_hat : 0.0;
  const double d_C_hat = (psi_CD_hat > zero_tol) ? psi_CN_hat / psi_CD_hat : 0.0;
  const double d_F_hat = (psi_FD_hat > zero_tol) ? psi_FN_hat / psi_FD_hat : 0.0;

  RealVec estimates(2);
  estimates(0) = (d_C_hat > zero_tol)
                     ? 1.0 - d_hat / d_C_hat
                     : std::numeric_limits<double>::quiet_NaN();
  estimates(1) = (d_F_hat > zero_tol)
                     ? 1.0 - d_hat / d_F_hat
                     : std::numeric_limits<double>::quiet_NaN();

  // --- Influence functions ---
  RealVec phi_dN = h_dN.array() - psi_dN_hat;
  RealVec phi_dD = h_dD.array() - psi_dD_hat;

  RealVec phi_CN = kernel_CN.influence(psi_CN_hat, n);
  RealVec phi_CD = kernel_CD.influence(psi_CD_hat, n);
  RealVec phi_FN = kernel_FN.influence(psi_FN_hat, n);
  RealVec phi_FD = kernel_FD.influence(psi_FD_hat, n);

  RealMat phi_matrix(n, 6);
  phi_matrix.col(0) = phi_dN;
  phi_matrix.col(1) = phi_dD;
  phi_matrix.col(2) = phi_CN;
  phi_matrix.col(3) = phi_CD;
  phi_matrix.col(4) = phi_FN;
  phi_matrix.col(5) = phi_FD;
  RealMat Gamma_hat = (phi_matrix.transpose() * phi_matrix) / static_cast<double>(n);

  // --- Delta method: disagreement covariance Sigma_hat ---
  RealMat J_d = RealMat::Zero(3, 6);
  if (psi_dD_hat > zero_tol) {
    J_d(0, 0) = 1.0 / psi_dD_hat;
    J_d(0, 1) = -psi_dN_hat / (psi_dD_hat * psi_dD_hat);
  }
  if (psi_CD_hat > zero_tol) {
    J_d(1, 2) = 1.0 / psi_CD_hat;
    J_d(1, 3) = -psi_CN_hat / (psi_CD_hat * psi_CD_hat);
  }
  if (psi_FD_hat > zero_tol) {
    J_d(2, 4) = 1.0 / psi_FD_hat;
    J_d(2, 5) = -psi_FN_hat / (psi_FD_hat * psi_FD_hat);
  }
  RealMat Sigma_hat = J_d * Gamma_hat * J_d.transpose();

  // --- Delta method: kappa covariance (Conger, Fleiss) ---
  RealMat J_kappa = RealMat::Zero(2, 3);
  if (d_C_hat > zero_tol) {
    J_kappa(0, 0) = -1.0 / d_C_hat;
    J_kappa(0, 1) = d_hat / (d_C_hat * d_C_hat);
  }
  if (d_F_hat > zero_tol) {
    J_kappa(1, 0) = -1.0 / d_F_hat;
    J_kappa(1, 2) = d_hat / (d_F_hat * d_F_hat);
  }
  RealMat kappa_cov = (J_kappa * Sigma_hat * J_kappa.transpose()) / static_cast<double>(n);

  return Estimation{std::move(estimates), std::move(kappa_cov)};
}

}  // namespace

Result<Estimation> estimate_available_continuous(
    RealMatView ratings, loss::ContinuousLoss loss) {
  return estimate_continuous(ratings, loss, detail::Reweighting::available);
}

Result<Estimation> estimate_ipw_continuous(
    RealMatView ratings, loss::ContinuousLoss loss) {
  return estimate_continuous(ratings, loss, detail::Reweighting::ipw);
}

Result<Estimation> estimate_gwet_continuous(
    RealMatView ratings, loss::ContinuousLoss loss) {
  return estimate_continuous(ratings, loss, detail::Reweighting::gwet);
}

}  // namespace misskappa
