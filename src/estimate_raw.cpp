#include "detail_estimate_raw.hpp"
#include "detail_inverse_weights.hpp"
#include "misskappa/estimate.hpp"

#include <cmath>
#include <limits>

namespace misskappa {

namespace {

constexpr double zero_tol = 1e-9;

struct KernelMoments {
  RealVec row_sum;
  RealVec col_sum;
  double total = 0.0;

  explicit KernelMoments(int n)
      : row_sum(RealVec::Zero(n)), col_sum(RealVec::Zero(n)) {}

  void add(int row, int col, double value) {
    row_sum(row) += value;
    col_sum(col) += value;
    total += value;
  }

  double mean(int n) const {
    return total / (static_cast<double>(n) * n);
  }

  RealVec influence(double psi, int n) const {
    const double inv_n = 1.0 / static_cast<double>(n);
    return ((row_sum.array() * inv_n - psi)
            + (col_sum.array() * inv_n - psi)).matrix();
  }
};

// Build mask: 1 if observed, 0 if na_code. Returns also (n, R) for convenience.
Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>
build_mask(IntMatView ratings) {
  Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> m(ratings.rows(), ratings.cols());
  for (Eigen::Index i = 0; i < ratings.rows(); ++i) {
    for (Eigen::Index j = 0; j < ratings.cols(); ++j) {
      m(i, j) = (ratings(i, j) == na_code) ? 0 : 1;
    }
  }
  return m;
}

}  // namespace

namespace detail {

Result<Estimation> estimate_raw(
    IntMatView ratings, RealMatView weights, Reweighting mode) {
  const int n = static_cast<int>(ratings.rows());
  const int R = static_cast<int>(ratings.cols());
  if (n < 1) return std::unexpected(Error::invalid_argument);
  if (R < 2) return std::unexpected(Error::invalid_argument);

  // Validate category indices: must be na_code or in [0, C-1] where C is
  // the dimension of the weight matrix.
  const int C = static_cast<int>(weights.rows());
  if (weights.cols() != C) return std::unexpected(Error::dimension_mismatch);
  if (C < 1) return std::unexpected(Error::invalid_argument);
  bool any_observed = false;
  for (Eigen::Index i = 0; i < ratings.rows(); ++i) {
    for (Eigen::Index j = 0; j < ratings.cols(); ++j) {
      const int x = ratings(i, j);
      if (x == na_code) continue;
      any_observed = true;
      if (x < 0 || x >= C) return std::unexpected(Error::invalid_argument);
    }
  }
  if (!any_observed) return std::unexpected(Error::invalid_argument);

  const auto mask = build_mask(ratings);
  auto wres = compute_inverse_weights(mask, n, R, mode);
  if (!wres) return std::unexpected(wres.error());
  const RealVec& pi_j_inv = wres->pi_j_inv;
  const RealMat& pi_jk_inv = wres->pi_jk_inv;

  // Convention: `weights` is the AGREEMENT matrix (1 on diagonal), matching
  // irrCAC's identity.weights / quadratic.weights and the standard
  // weighted-kappa literature. The kappa formulas below run on the
  // disagreement matrix L = 1 - W.
  const RealMat L = RealMat::Constant(C, C, 1.0) - weights;

  // --- Per-subject U-statistic kernels (observed-disagreement) ---
  // h_dN_i = sum over j<k of M_ij M_ik * loss(x_ij, x_ik) * pi_jk_inv(j, k)
  // h_dD_i = sum over j<k of M_ij M_ik * pi_jk_inv(j, k)
  RealVec h_dN = RealVec::Zero(n);
  RealVec h_dD = RealVec::Zero(n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < R - 1; ++j) {
      if (!mask(i, j)) continue;
      for (int k = j + 1; k < R; ++k) {
        if (!mask(i, k)) continue;
        const double w_jk = pi_jk_inv(j, k);
        const int a = ratings(i, j);
        const int b = ratings(i, k);
        h_dN(i) += L(a, b) * w_jk;
        h_dD(i) += w_jk;
      }
    }
  }
  const double psi_dN_hat = h_dN.mean();
  const double psi_dD_hat = h_dD.mean();

  // --- V-statistic kernels (chance disagreement) ---
  // kernel_CN(i, i') = sum over j<k of M_ij M_{i'k} * loss(x_ij, x_{i'k})
  //                                                 * pi_j_inv(j) * pi_j_inv(k)
  // kernel_FN(i, i') = sum over all j, k same, but without the j<k restriction.
  KernelMoments kernel_CN(n);
  KernelMoments kernel_CD(n);
  KernelMoments kernel_FN(n);
  KernelMoments kernel_FD(n);
  for (int i = 0; i < n; ++i) {
    for (int ip = 0; ip < n; ++ip) {
      double h_cn = 0, h_cd = 0, h_fn = 0, h_fd = 0;
      for (int j = 0; j < R; ++j) {
        if (!mask(i, j)) continue;
        const double w_j = pi_j_inv(j);
        const int a = ratings(i, j);
        for (int k = 0; k < R; ++k) {
          if (!mask(ip, k)) continue;
          const double w_jk = w_j * pi_j_inv(k);
          const int b = ratings(ip, k);
          const double l = L(a, b);
          // Fleiss: all (j, k).
          h_fn += l * w_jk;
          h_fd += w_jk;
          // Conger: j < k only.
          if (j < k) {
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

  // --- Point estimates ---
  const double d_hat   = (psi_dD_hat > zero_tol) ? psi_dN_hat / psi_dD_hat : 0.0;
  const double d_C_hat = (psi_CD_hat > zero_tol) ? psi_CN_hat / psi_CD_hat : 0.0;
  const double d_F_hat = (psi_FD_hat > zero_tol) ? psi_FN_hat / psi_FD_hat : 0.0;
  const double d_BP    = L.sum() / (static_cast<double>(C) * C);

  RealVec estimates(3);
  estimates(0) = (d_C_hat > zero_tol)
                     ? 1.0 - d_hat / d_C_hat
                     : std::numeric_limits<double>::quiet_NaN();
  estimates(1) = (d_F_hat > zero_tol)
                     ? 1.0 - d_hat / d_F_hat
                     : std::numeric_limits<double>::quiet_NaN();
  estimates(2) = (d_BP > zero_tol)
                     ? 1.0 - d_hat / d_BP
                     : std::numeric_limits<double>::quiet_NaN();

  // --- Influence functions ---
  // U-statistic IFs: phi_i = h_i - psi.
  RealVec phi_dN = h_dN.array() - psi_dN_hat;
  RealVec phi_dD = h_dD.array() - psi_dD_hat;

  // V-statistic IFs: phi_i = (mean_j K(i, j) - psi) + (mean_j K(j, i) - psi).
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

  // --- Delta method 1: disagreement covariance Sigma_hat ---
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

  // --- Delta method 2: kappa covariance ---
  // d_BP is a constant (depends only on the weight matrix, not the data),
  // so its row in J_kappa is just -1 / d_BP times the d-coordinate.
  RealMat J_kappa = RealMat::Zero(3, 3);
  if (d_C_hat > zero_tol) {
    J_kappa(0, 0) = -1.0 / d_C_hat;
    J_kappa(0, 1) = d_hat / (d_C_hat * d_C_hat);
  }
  if (d_F_hat > zero_tol) {
    J_kappa(1, 0) = -1.0 / d_F_hat;
    J_kappa(1, 2) = d_hat / (d_F_hat * d_F_hat);
  }
  if (d_BP > zero_tol) {
    J_kappa(2, 0) = -1.0 / d_BP;
  }

  RealMat kappa_cov = (J_kappa * Sigma_hat * J_kappa.transpose()) / static_cast<double>(n);

  return Estimation{std::move(estimates), std::move(kappa_cov)};
}

}  // namespace detail

// --- Public entry points -----------------------------------------------------

Result<Estimation> estimate_available(IntMatView ratings, RealMatView weights) {
  return detail::estimate_raw(ratings, weights, detail::Reweighting::available);
}

Result<Estimation> estimate_ipw(IntMatView ratings, RealMatView weights) {
  return detail::estimate_raw(ratings, weights, detail::Reweighting::ipw);
}

Result<Estimation> estimate_gwet(IntMatView ratings, RealMatView weights) {
  return detail::estimate_raw(ratings, weights, detail::Reweighting::gwet);
}

}  // namespace misskappa
