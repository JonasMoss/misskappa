// Closed-form quadratic-loss estimator (Conger, Fleiss, Brennan-Prediger)
// for raw real-valued ratings. The counts-format counterpart is in
// src/estimate_quadratic_counts.cpp.
//
// Treats categorical ratings as numeric scores; missing entries are NaN.
// Math ported from dev/legacy/misskappa/src/kappaqp.cpp onto Eigen +
// Result<T>. Per-rater means and covariances are computed pairwise from
// observed entries; the asymptotic covariance of those moments is
// constructed from third- and fourth-order moments, and the delta method
// gives the asymptotic covariance of the three kappa coefficients.

#include "misskappa/estimate.hpp"

#include <cmath>
#include <limits>
#include <vector>

namespace misskappa {

namespace {

constexpr double zero_tol = 1e-9;

// Internal moment-bundle output of CalculatePsi.
struct AcovResults {
  RealMat psi;       // 3 x 3 covariance of (t1, t2, t3) summary statistics
  RealVec mu_hat;    // length R per-rater means
  RealMat sigma_hat; // R x R per-rater-pair covariances
  int R = 0;
};

// Build mask, M(i, j) = 1 if x(i, j) is finite, else 0.
Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>
build_finite_mask(RealMatView x) {
  Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> m(x.rows(), x.cols());
  for (Eigen::Index i = 0; i < x.rows(); ++i) {
    for (Eigen::Index j = 0; j < x.cols(); ++j) {
      m(i, j) = std::isfinite(x(i, j)) ? 1 : 0;
    }
  }
  return m;
}

// Heavy-lifter: from raw real-valued ratings + finiteness mask, compute the
// 3x3 covariance of (t1 = sum(sigma_hat), t2 = trace(sigma_hat),
// t3 = sum((mu - mean(mu))^2)) using 2nd / 3rd / 4th order moments.
Result<AcovResults> calculate_psi(
    const RealMat& x, const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>& M,
    int n_eff) {
  const int R = static_cast<int>(x.cols());
  AcovResults out;
  out.R = R;

  // Zero out non-finite entries in a working copy for later linear algebra.
  RealMat Xf = x;
  for (Eigen::Index i = 0; i < Xf.rows(); ++i) {
    for (Eigen::Index j = 0; j < Xf.cols(); ++j) {
      if (M(i, j) == 0) Xf(i, j) = 0.0;
    }
  }

  RealVec count = M.cast<double>().colwise().sum().transpose();  // length R
  out.mu_hat.resize(R);
  for (int j = 0; j < R; ++j) {
    out.mu_hat(j) = (count(j) > 0) ? Xf.col(j).sum() / count(j)
                                   : std::numeric_limits<double>::quiet_NaN();
  }

  // Y_hat = X with the rater-mean subtracted, NaNs zeroed out.
  RealMat Y_hat = x;
  for (Eigen::Index i = 0; i < Y_hat.rows(); ++i) {
    for (Eigen::Index j = 0; j < Y_hat.cols(); ++j) {
      Y_hat(i, j) = M(i, j) ? (x(i, j) - out.mu_hat(j)) : 0.0;
    }
  }

  RealVec p1 = count / static_cast<double>(n_eff);
  RealMat p2 = (M.cast<double>().transpose() * M.cast<double>()) / static_cast<double>(n_eff);

  // sigma_hat(j, k): pairwise covariance of raters j, k over rows where both observed.
  out.sigma_hat = RealMat::Zero(R, R);
  for (int j = 0; j < R; ++j) {
    for (int k = j; k < R; ++k) {
      // Find row indices where both raters observed.
      std::vector<int> idx;
      idx.reserve(n_eff);
      for (int i = 0; i < n_eff; ++i) {
        if (M(i, j) == 1 && M(i, k) == 1) idx.push_back(i);
      }
      if (!idx.empty()) {
        double sum_j = 0, sum_k = 0;
        for (int i : idx) { sum_j += x(i, j); sum_k += x(i, k); }
        const double mu_j_p = sum_j / idx.size();
        const double mu_k_p = sum_k / idx.size();
        double sum_prod = 0;
        for (int i : idx) {
          sum_prod += (x(i, j) - mu_j_p) * (x(i, k) - mu_k_p);
        }
        out.sigma_hat(j, k) = sum_prod / idx.size();
      }
      out.sigma_hat(k, j) = out.sigma_hat(j, k);
    }
  }

  // 3rd / 4th order moments and matching observation probabilities.
  // mu3[i][j][k] = E[Y_hat(:, i) * Y_hat(:, j) * Y_hat(:, k)] over rows where
  //               all three are observed.
  // mu4[i][j][k][l] = E[Y * Y * Y * Y] similarly.
  // p3, p4 are the corresponding observed-fractions.
  std::vector<double> p3(R * R * R, 0.0);
  std::vector<double> mu3(R * R * R, 0.0);
  std::vector<double> p4(R * R * R * R, 0.0);
  std::vector<double> mu4(R * R * R * R, 0.0);
  auto idx3 = [R](int i, int j, int k) { return ((i * R) + j) * R + k; };
  auto idx4 = [R](int i, int j, int k, int l) { return (((i * R) + j) * R + k) * R + l; };

  for (int i = 0; i < R; ++i) {
    for (int j = 0; j < R; ++j) {
      for (int k = 0; k < R; ++k) {
        std::vector<int> rows3;
        rows3.reserve(n_eff);
        for (int row = 0; row < n_eff; ++row) {
          if (M(row, i) && M(row, j) && M(row, k)) rows3.push_back(row);
        }
        p3[idx3(i, j, k)] = static_cast<double>(rows3.size()) / n_eff;
        if (!rows3.empty()) {
          double s = 0;
          for (int row : rows3) s += Y_hat(row, i) * Y_hat(row, j) * Y_hat(row, k);
          mu3[idx3(i, j, k)] = s / rows3.size();
        }
        for (int l = 0; l < R; ++l) {
          std::vector<int> rows4;
          rows4.reserve(n_eff);
          for (int row : rows3) {
            if (M(row, l)) rows4.push_back(row);
          }
          p4[idx4(i, j, k, l)] = static_cast<double>(rows4.size()) / n_eff;
          if (!rows4.empty()) {
            double s = 0;
            for (int row : rows4) {
              s += Y_hat(row, i) * Y_hat(row, j) * Y_hat(row, k) * Y_hat(row, l);
            }
            mu4[idx4(i, j, k, l)] = s / rows4.size();
          }
        }
      }
    }
  }

  // Assemble Psi (3 x 3 covariance of (t1, t2, t3)).
  out.psi = RealMat::Zero(3, 3);
  RealVec v_dmu = 2.0 * (out.mu_hat.array() - out.mu_hat.mean()).matrix();

  for (int i = 0; i < R; ++i) {
    for (int j = 0; j < R; ++j) {
      for (int k = 0; k < R; ++k) {
        for (int l = 0; l < R; ++l) {
          if (p2(i, j) > zero_tol && p2(k, l) > zero_tol) {
            const double gamma = mu4[idx4(i, j, k, l)] - out.sigma_hat(i, j) * out.sigma_hat(k, l);
            const double term = (p4[idx4(i, j, k, l)] / (p2(i, j) * p2(k, l))) * gamma;
            out.psi(0, 0) += term;
            if (k == l) out.psi(0, 1) += term;
            if (i == j && k == l) out.psi(1, 1) += term;
          }
        }
      }
    }
  }
  out.psi(1, 0) = out.psi(0, 1);

  for (int i = 0; i < R; ++i) {
    for (int j = 0; j < R; ++j) {
      if (p1(i) > zero_tol && p1(j) > zero_tol) {
        out.psi(2, 2) += v_dmu(i) * v_dmu(j) * (p2(i, j) / (p1(i) * p1(j))) * out.sigma_hat(i, j);
      }
    }
  }

  for (int i = 0; i < R; ++i) {
    for (int j = 0; j < R; ++j) {
      for (int k = 0; k < R; ++k) {
        if (p1(i) > zero_tol && p2(j, k) > zero_tol) {
          const double omega = (p3[idx3(i, j, k)] / (p1(i) * p2(j, k))) * mu3[idx3(i, j, k)];
          out.psi(0, 2) += v_dmu(i) * omega;
          if (j == k) out.psi(1, 2) += v_dmu(i) * omega;
        }
      }
    }
  }
  out.psi(2, 0) = out.psi(0, 2);
  out.psi(2, 1) = out.psi(1, 2);

  return out;
}

}  // namespace

Result<Estimation> estimate_quadratic(RealMatView ratings, const RealVec& values) {
  const int C = static_cast<int>(values.size());
  if (C < 1) return std::unexpected(Error::invalid_argument);
  if (ratings.rows() < 1 || ratings.cols() < 2) {
    return std::unexpected(Error::invalid_argument);
  }

  // Brennan-Prediger constant: 2 / C^2 * (C * sum(v^2) - sum(v)^2).
  const double c1 = (2.0 / (C * C)) * (C * values.squaredNorm() - std::pow(values.sum(), 2));

  // Effective rows: at least 2 ratings observed.
  const auto M_full = build_finite_mask(ratings);
  std::vector<int> eff_rows;
  for (Eigen::Index i = 0; i < ratings.rows(); ++i) {
    int obs = 0;
    for (Eigen::Index j = 0; j < ratings.cols(); ++j) obs += M_full(i, j);
    if (obs > 1) eff_rows.push_back(static_cast<int>(i));
  }
  if (eff_rows.size() < 2) return std::unexpected(Error::invalid_argument);

  RealMat x_eff(eff_rows.size(), ratings.cols());
  Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> M_eff(eff_rows.size(), ratings.cols());
  for (std::size_t k = 0; k < eff_rows.size(); ++k) {
    x_eff.row(k) = ratings.row(eff_rows[k]);
    M_eff.row(k) = M_full.row(eff_rows[k]);
  }
  const int n_eff = static_cast<int>(eff_rows.size());
  const int R = static_cast<int>(ratings.cols());

  auto acov = calculate_psi(x_eff, M_eff, n_eff);
  if (!acov) return std::unexpected(acov.error());

  // Summary stats t1 = sum(Sigma), t2 = tr(Sigma), t3 = sum((mu - mean(mu))^2).
  const double t1 = acov->sigma_hat.sum();
  const double t2 = acov->sigma_hat.diagonal().sum();
  RealVec mu_dev = acov->mu_hat.array() - acov->mu_hat.mean();
  const double t3 = mu_dev.squaredNorm();

  // Numerators / denominators for Fleiss and Conger.
  const double NF = t1 - t2 - t3;
  const double DF = (R - 1.0) * (t2 + t3);
  const double NC = t1 - t2;
  const double DC = (R - 1.0) * t2 + static_cast<double>(R) * t3;

  const double fleiss_est = (std::abs(DF) > zero_tol)
                                ? NF / DF
                                : std::numeric_limits<double>::quiet_NaN();
  const double conger_est = (std::abs(DC) > zero_tol)
                                ? NC / DC
                                : std::numeric_limits<double>::quiet_NaN();
  const double d_obs = (2.0 / (R - 1.0)) * (t2 + t3 - t1 / R);
  const double bp_est = (std::abs(c1) > zero_tol)
                            ? 1.0 - d_obs / c1
                            : std::numeric_limits<double>::quiet_NaN();

  // Gradients of (Fleiss, Conger, BP) wrt (t1, t2, t3).
  RealVec gF = RealVec::Zero(3);
  RealVec gC = RealVec::Zero(3);
  RealVec gBP = RealVec::Zero(3);
  if (std::abs(DF) > zero_tol) {
    gF(0) = 1.0 / DF;
    gF(1) = (-DF - NF * (R - 1.0)) / (DF * DF);
    gF(2) = (-DF - NF * (R - 1.0)) / (DF * DF);
  }
  if (std::abs(DC) > zero_tol) {
    gC(0) = 1.0 / DC;
    gC(1) = (-DC - NC * (R - 1.0)) / (DC * DC);
    gC(2) = (-NC * static_cast<double>(R)) / (DC * DC);
  }
  if (std::abs(c1) > zero_tol) {
    const double m = -(2.0 / (c1 * (R - 1.0)));
    gBP(0) = m * (-1.0 / R);
    gBP(1) = m;
    gBP(2) = m;
  }

  // scaled_acov = G^T * psi * G, in (Fleiss, Conger, BP) order, then divide
  // by n_eff. Reorder rows/cols to (Conger, Fleiss, BP) so the output matches
  // the misskappa raw-categorical convention.
  RealMat G(3, 3);
  G.col(0) = gF;
  G.col(1) = gC;
  G.col(2) = gBP;
  const RealMat scaled = G.transpose() * acov->psi * G;  // (Fleiss, Conger, BP) ordering
  RealMat vcov_reordered(3, 3);
  const int perm[3] = {1, 0, 2};  // map (Conger, Fleiss, BP) -> (Fleiss, Conger, BP)
  for (int r = 0; r < 3; ++r) {
    for (int c = 0; c < 3; ++c) {
      vcov_reordered(r, c) = scaled(perm[r], perm[c]);
    }
  }
  const RealMat vcov = vcov_reordered / static_cast<double>(n_eff);

  RealVec estimates(3);
  estimates(0) = conger_est;
  estimates(1) = fleiss_est;
  estimates(2) = bp_est;
  return Estimation{std::move(estimates), std::move(vcov), {}};
}

}  // namespace misskappa
