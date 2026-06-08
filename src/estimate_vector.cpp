#include "detail_kernel_moments.hpp"
#include "misskappa/estimate.hpp"

#include <cmath>
#include <limits>
#include <vector>

namespace misskappa {

namespace {

constexpr double zero_tol = 1e-9;

enum class VectorReweighting {
  pairwise,
  ipw,
};

struct VectorWeights {
  RealMat pi_inv;              // R x p
  std::vector<double> pair_inv; // (R * R) x p, row-major in (j, k, feature).
};

int col_index(int rater, int feature, int features) {
  return rater * features + feature;
}

int pair_index(int rater_a, int rater_b, int feature, int R, int features) {
  return (rater_a * R + rater_b) * features + feature;
}

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

bool valid_loss(const loss::ComponentSeparableLoss& loss, int features) {
  if (loss.feature_weights.size() != static_cast<Eigen::Index>(features)) return false;
  if (loss.compute == nullptr || loss.transform == nullptr ||
      loss.transform_derivative == nullptr) {
    return false;
  }
  double total = 0.0;
  for (int l = 0; l < features; ++l) {
    const double w = loss.feature_weights(l);
    if (!std::isfinite(w) || w < 0.0) return false;
    total += w;
  }
  return total > zero_tol;
}

Result<VectorWeights> compute_vector_weights(
    const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>& mask,
    int n, int R, int features,
    const loss::ComponentSeparableLoss& loss,
    VectorReweighting mode) {
  VectorWeights out;
  out.pi_inv = RealMat::Ones(R, features);
  out.pair_inv.assign(static_cast<std::size_t>(R * R * features), 1.0);
  if (mode == VectorReweighting::pairwise) return out;

  for (int j = 0; j < R; ++j) {
    for (int l = 0; l < features; ++l) {
      if (loss.feature_weights(l) <= zero_tol) continue;
      int count = 0;
      const int c = col_index(j, l, features);
      for (int i = 0; i < n; ++i) count += mask(i, c);
      if (count == 0) return misskappa::unexpected(Error::singular_weight);
      out.pi_inv(j, l) = static_cast<double>(n) / static_cast<double>(count);
    }
  }

  for (int j = 0; j < R; ++j) {
    for (int k = 0; k < R; ++k) {
      for (int l = 0; l < features; ++l) {
        if (loss.feature_weights(l) <= zero_tol) continue;
        int count = 0;
        const int cj = col_index(j, l, features);
        const int ck = col_index(k, l, features);
        for (int i = 0; i < n; ++i) {
          if (mask(i, cj) && mask(i, ck)) ++count;
        }
        const int idx = pair_index(j, k, l, R, features);
        out.pair_inv[static_cast<std::size_t>(idx)] =
            (count > 0) ? static_cast<double>(n) / static_cast<double>(count) : 0.0;
      }
    }
  }

  return out;
}

Result<Estimation> estimate_vector(
    RealMatView ratings, int features,
    loss::ComponentSeparableLoss loss,
    VectorReweighting mode) {
  const int n = static_cast<int>(ratings.rows());
  const int cols = static_cast<int>(ratings.cols());
  if (n < 1) return misskappa::unexpected(Error::invalid_argument);
  if (features < 1) return misskappa::unexpected(Error::invalid_argument);
  if (cols < 1 || cols % features != 0) {
    return misskappa::unexpected(Error::dimension_mismatch);
  }
  const int R = cols / features;
  if (R < 2) return misskappa::unexpected(Error::invalid_argument);
  if (!valid_loss(loss, features)) return misskappa::unexpected(Error::invalid_argument);

  const auto mask = build_finite_mask(ratings);
  auto wres = compute_vector_weights(mask, n, R, features, loss, mode);
  if (!wres) return misskappa::unexpected(wres.error());
  const VectorWeights& weights = *wres;

  RealVec h_dN = RealVec::Zero(n);
  RealVec h_dD = RealVec::Zero(n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < R - 1; ++j) {
      for (int k = j + 1; k < R; ++k) {
        for (int l = 0; l < features; ++l) {
          const double feature_weight = loss.feature_weights(l);
          if (feature_weight <= zero_tol) continue;
          const int cj = col_index(j, l, features);
          const int ck = col_index(k, l, features);
          if (!mask(i, cj) || !mask(i, ck)) continue;
          const double w = feature_weight *
              weights.pair_inv[static_cast<std::size_t>(pair_index(j, k, l, R, features))];
          h_dN(i) += w * loss.compute(ratings(i, cj), ratings(i, ck));
          h_dD(i) += w;
        }
      }
    }
  }
  const double psi_dN_hat = h_dN.mean();
  const double psi_dD_hat = h_dD.mean();
  if (psi_dD_hat <= zero_tol) return misskappa::unexpected(Error::invalid_argument);

  detail::KernelMoments kernel_CN(n);
  detail::KernelMoments kernel_CD(n);
  detail::KernelMoments kernel_FN(n);
  detail::KernelMoments kernel_FD(n);
  for (int i = 0; i < n; ++i) {
    for (int ip = 0; ip < n; ++ip) {
      double h_cn = 0.0;
      double h_cd = 0.0;
      double h_fn = 0.0;
      double h_fd = 0.0;
      for (int j = 0; j < R; ++j) {
        for (int k = 0; k < R; ++k) {
          for (int l = 0; l < features; ++l) {
            const double feature_weight = loss.feature_weights(l);
            if (feature_weight <= zero_tol) continue;
            const int cj = col_index(j, l, features);
            const int ck = col_index(k, l, features);
            if (!mask(i, cj) || !mask(ip, ck)) continue;
            const double w = feature_weight * weights.pi_inv(j, l) * weights.pi_inv(k, l);
            const double value = w * loss.compute(ratings(i, cj), ratings(ip, ck));
            h_fn += value;
            h_fd += w;
            if (j < k) {
              h_cn += value;
              h_cd += w;
            }
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
  if (psi_CD_hat <= zero_tol || psi_FD_hat <= zero_tol) {
    return misskappa::unexpected(Error::invalid_argument);
  }

  const double raw_d   = psi_dN_hat / psi_dD_hat;
  const double raw_C   = psi_CN_hat / psi_CD_hat;
  const double raw_F   = psi_FN_hat / psi_FD_hat;
  const double d_hat   = loss.transform(raw_d);
  const double d_C_hat = loss.transform(raw_C);
  const double d_F_hat = loss.transform(raw_F);

  RealVec estimates(2);
  estimates(0) = (d_C_hat > zero_tol)
                     ? 1.0 - d_hat / d_C_hat
                     : std::numeric_limits<double>::quiet_NaN();
  estimates(1) = (d_F_hat > zero_tol)
                     ? 1.0 - d_hat / d_F_hat
                     : std::numeric_limits<double>::quiet_NaN();

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

  RealMat J_d = RealMat::Zero(3, 6);
  auto fill_ratio_derivative = [&](int row, int n_col, int d_col,
                                   double num, double den, double raw) {
    if (den <= zero_tol) return;
    const double deriv = loss.transform_derivative(raw);
    if (!std::isfinite(deriv)) return;
    J_d(row, n_col) = deriv / den;
    J_d(row, d_col) = -deriv * num / (den * den);
  };
  fill_ratio_derivative(0, 0, 1, psi_dN_hat, psi_dD_hat, raw_d);
  fill_ratio_derivative(1, 2, 3, psi_CN_hat, psi_CD_hat, raw_C);
  fill_ratio_derivative(2, 4, 5, psi_FN_hat, psi_FD_hat, raw_F);

  RealMat J_kappa = RealMat::Zero(2, 3);
  if (d_C_hat > zero_tol) {
    J_kappa(0, 0) = -1.0 / d_C_hat;
    J_kappa(0, 1) = d_hat / (d_C_hat * d_C_hat);
  }
  if (d_F_hat > zero_tol) {
    J_kappa(1, 0) = -1.0 / d_F_hat;
    J_kappa(1, 2) = d_hat / (d_F_hat * d_F_hat);
  }

  const RealMat J_combined = J_kappa * J_d;
  RealMat kappa_cov =
      (J_combined * Gamma_hat * J_combined.transpose()) / static_cast<double>(n);
  RealMat psi_kappa = build_psi_from_phi(phi_matrix, J_combined);

  return Estimation{std::move(estimates), std::move(kappa_cov), std::move(psi_kappa)};
}

}  // namespace

Result<Estimation> estimate_pairwise_vector(
    RealMatView ratings, int features, loss::ComponentSeparableLoss loss) {
  return estimate_vector(ratings, features, std::move(loss), VectorReweighting::pairwise);
}

Result<Estimation> estimate_ipw_vector(
    RealMatView ratings, int features, loss::ComponentSeparableLoss loss) {
  return estimate_vector(ratings, features, std::move(loss), VectorReweighting::ipw);
}

}  // namespace misskappa
