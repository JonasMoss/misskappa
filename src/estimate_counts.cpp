// Counts-format available-case estimator.
//
// Input is n x C of non-negative integer counts: counts(i, k) is the number
// of raters who assigned subject i to category k. r_i = sum_k counts(i, k)
// is the number of raters for subject i (may vary across subjects).
//
// Same moment-based pipeline as estimate_raw / estimate_continuous, but
// with closed-form U/V-statistic kernels that exploit the count structure:
// the within-subject pair sum becomes 0.5 * (N_u^T L N_u - diag(L) . N_u)
// and the chance kernel between subjects is N_u^T L N_v.
//
// Returned coefficients: (Fleiss, Brennan-Prediger). No Conger (raters are
// not identified in this input format). No IPW / Gwet (per-rater
// observation rates are aggregated away by counts).

#include "misskappa/estimate.hpp"

#include <cmath>
#include <limits>

namespace misskappa {

namespace {

constexpr double zero_tol = 1e-9;

}  // namespace

Result<Estimation> estimate_available_counts(IntMatView counts, RealMatView weights) {
  const int n = static_cast<int>(counts.rows());
  const int C = static_cast<int>(counts.cols());
  if (n < 1) return std::unexpected(Error::invalid_argument);
  if (weights.rows() != C || weights.cols() != C) {
    return std::unexpected(Error::dimension_mismatch);
  }
  if (C < 1) return std::unexpected(Error::invalid_argument);
  // Validate non-negative integer counts.
  for (Eigen::Index i = 0; i < counts.rows(); ++i) {
    for (Eigen::Index k = 0; k < counts.cols(); ++k) {
      if (counts(i, k) < 0) return std::unexpected(Error::invalid_argument);
    }
  }

  // Agreement convention on input; work on disagreement L = 1 - W.
  const RealMat L = RealMat::Constant(C, C, 1.0) - weights;

  // counts as doubles for the linear algebra.
  RealMat N = counts.cast<double>();        // n x C
  RealVec r = N.rowwise().sum();            // n; r(i) = total raters on subject i.
  RealVec L_diag = L.diagonal();

  // --- Per-subject U-statistic for observed disagreement ---
  // h_dN(i) = 0.5 * (N_i^T L N_i - diag(L) . N_i)  ; number of within-subject
  //                                                 rater-pair disagreement
  // h_dD(i) = r(i) * (r(i) - 1) / 2                ; number of rater pairs
  RealMat NLN = N * L * N.transpose();       // n x n; we only need diagonal.
  RealVec h_dN(n);
  RealVec h_dD(n);
  for (int i = 0; i < n; ++i) {
    h_dN(i) = 0.5 * (NLN(i, i) - L_diag.dot(N.row(i).transpose()));
    h_dD(i) = r(i) * (r(i) - 1.0) / 2.0;
  }
  const double psi_dN_hat = h_dN.mean();
  const double psi_dD_hat = h_dD.mean();

  // --- V-statistic kernels (chance disagreement) ---
  // kernel_FN(i, ip) = N_i^T L N_ip   (already in NLN)
  // kernel_FD(i, ip) = r(i) * r(ip)
  const RealMat& kernel_FN = NLN;
  RealMat kernel_FD = r * r.transpose();
  const double psi_FN_hat = kernel_FN.sum() / (static_cast<double>(n) * n);
  const double psi_FD_hat = kernel_FD.sum() / (static_cast<double>(n) * n);

  // --- Point estimates ---
  const double d_hat   = (psi_dD_hat > zero_tol) ? psi_dN_hat / psi_dD_hat : 0.0;
  const double d_F_hat = (psi_FD_hat > zero_tol) ? psi_FN_hat / psi_FD_hat : 0.0;
  const double d_BP    = L.sum() / (static_cast<double>(C) * C);

  RealVec estimates(2);
  estimates(0) = (d_F_hat > zero_tol)
                     ? 1.0 - d_hat / d_F_hat
                     : std::numeric_limits<double>::quiet_NaN();
  estimates(1) = (d_BP > zero_tol)
                     ? 1.0 - d_hat / d_BP
                     : std::numeric_limits<double>::quiet_NaN();

  // --- Influence functions ---
  RealVec phi_dN = h_dN.array() - psi_dN_hat;
  RealVec phi_dD = h_dD.array() - psi_dD_hat;

  auto v_stat_if = [](const RealMat& K, double psi) {
    RealVec row_mean = K.rowwise().mean();
    RealVec col_mean = K.colwise().mean().transpose();
    return RealVec((row_mean.array() - psi) + (col_mean.array() - psi));
  };
  RealVec phi_FN = v_stat_if(kernel_FN, psi_FN_hat);
  RealVec phi_FD = v_stat_if(kernel_FD, psi_FD_hat);

  RealMat phi_matrix(n, 4);
  phi_matrix.col(0) = phi_dN;
  phi_matrix.col(1) = phi_dD;
  phi_matrix.col(2) = phi_FN;
  phi_matrix.col(3) = phi_FD;
  RealMat Gamma_hat = (phi_matrix.transpose() * phi_matrix) / static_cast<double>(n);

  // --- Delta method ---
  RealMat J_d = RealMat::Zero(2, 4);
  if (psi_dD_hat > zero_tol) {
    J_d(0, 0) = 1.0 / psi_dD_hat;
    J_d(0, 1) = -psi_dN_hat / (psi_dD_hat * psi_dD_hat);
  }
  if (psi_FD_hat > zero_tol) {
    J_d(1, 2) = 1.0 / psi_FD_hat;
    J_d(1, 3) = -psi_FN_hat / (psi_FD_hat * psi_FD_hat);
  }
  RealMat Sigma_hat = J_d * Gamma_hat * J_d.transpose();

  RealMat J_kappa = RealMat::Zero(2, 2);
  if (d_F_hat > zero_tol) {
    J_kappa(0, 0) = -1.0 / d_F_hat;
    J_kappa(0, 1) = d_hat / (d_F_hat * d_F_hat);
  }
  if (d_BP > zero_tol) {
    J_kappa(1, 0) = -1.0 / d_BP;
  }
  RealMat kappa_cov = (J_kappa * Sigma_hat * J_kappa.transpose()) / static_cast<double>(n);

  return Estimation{std::move(estimates), std::move(kappa_cov)};
}

}  // namespace misskappa
