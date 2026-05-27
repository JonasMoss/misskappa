// Closed-form quadratic-loss estimator for counts-format input.
//
// Returns (Fleiss, Brennan-Prediger). The raw-input counterpart lives in
// src/estimate_quadratic.cpp; the math here is independent (no shared
// helpers) because the counts shape replaces the per-rater means and
// pairwise covariances with per-subject rescaled score sums.

#include "misskappa/estimate.hpp"

#include <cmath>
#include <limits>
#include <vector>

namespace misskappa {

namespace {

constexpr double zero_tol = 1e-9;

}  // namespace

Result<Estimation> estimate_quadratic_counts(IntMatView counts, const RealVec& values, int R_total) {
  const int C = static_cast<int>(values.size());
  if (C < 1) return std::unexpected(Error::invalid_argument);
  if (counts.cols() != C) return std::unexpected(Error::dimension_mismatch);
  if (R_total < 2) return std::unexpected(Error::invalid_argument);

  const double c1 = (2.0 / (C * C)) * (C * values.squaredNorm() - std::pow(values.sum(), 2));

  // Effective rows: at least one observed rating.
  std::vector<int> eff_rows;
  for (Eigen::Index i = 0; i < counts.rows(); ++i) {
    double row_sum = 0;
    for (Eigen::Index k = 0; k < counts.cols(); ++k) row_sum += counts(i, k);
    if (row_sum > 0) eff_rows.push_back(static_cast<int>(i));
  }
  if (eff_rows.size() < 2) return std::unexpected(Error::invalid_argument);

  const int n_eff = static_cast<int>(eff_rows.size());
  RealMat X_eff(n_eff, C);
  for (int k = 0; k < n_eff; ++k) {
    for (int c = 0; c < C; ++c) X_eff(k, c) = static_cast<double>(counts(eff_rows[k], c));
  }
  RealVec r_i = X_eff.rowwise().sum();
  RealVec v_sq(C);
  for (int c = 0; c < C; ++c) v_sq(c) = values(c) * values(c);

  // Per-subject summary stats.
  // s_i = R * (1 / r_i) * (X_i . values)             (sum of rater scores, rescaled to R raters)
  // q_i = R * (1 / r_i) * (X_i . values^2)           (sum of squared rater scores, rescaled)
  RealVec s_i(n_eff), q_i(n_eff);
  for (int k = 0; k < n_eff; ++k) {
    s_i(k) = (R_total / r_i(k)) * X_eff.row(k).dot(values);
    q_i(k) = (R_total / r_i(k)) * X_eff.row(k).dot(v_sq);
  }

  // Z(n_eff, 3) with columns (s, s^2, q).
  RealMat Z(n_eff, 3);
  Z.col(0) = s_i;
  Z.col(1) = s_i.array().square();
  Z.col(2) = q_i;
  RealVec phi = Z.colwise().mean().transpose();
  // Sample covariance (1/n divisor; matches arma::cov(_, 1)).
  RealMat Zc = Z.rowwise() - phi.transpose();
  RealMat Psi_dist = (Zc.transpose() * Zc) / static_cast<double>(n_eff);

  // Fleiss
  double est_f = std::numeric_limits<double>::quiet_NaN();
  RealVec grad_f = RealVec::Zero(3);
  const double Df = phi(2) - phi(0) * phi(0) / R_total;
  if (Df > zero_tol) {
    const double Nf = phi(1) - phi(0) * phi(0);
    const double r_inv = 1.0 / (R_total - 1.0);
    grad_f(0) = r_inv * ((-2.0 * phi(0) / Df) + (Nf * (2.0 * phi(0) / R_total)) / (Df * Df));
    grad_f(1) = r_inv / Df;
    grad_f(2) = r_inv * (-Nf / (Df * Df));
    est_f = r_inv * (Nf / Df - 1.0);
  }

  // Brennan-Prediger
  double est_bp = std::numeric_limits<double>::quiet_NaN();
  RealVec grad_bp = RealVec::Zero(3);
  if (std::abs(c1) > zero_tol) {
    const double d_bp = (2.0 / (R_total - 1.0)) * (phi(2) - phi(1) / R_total);
    est_bp = 1.0 - d_bp / c1;
    const double m = -(2.0 / (c1 * (R_total - 1.0)));
    grad_bp(0) = 0.0;
    grad_bp(1) = m * (-1.0 / R_total);
    grad_bp(2) = m;
  }

  RealMat G(3, 2);
  G.col(0) = grad_f;
  G.col(1) = grad_bp;

  RealVec estimates(2);
  estimates(0) = est_f;
  estimates(1) = est_bp;
  RealMat vcov = (G.transpose() * Psi_dist * G) / static_cast<double>(n_eff);
  return Estimation{std::move(estimates), std::move(vcov), {}};
}

}  // namespace misskappa
