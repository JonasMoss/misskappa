#include "misskappa/estimate.hpp"

#include <cmath>
#include <limits>
#include <vector>

namespace misskappa {

namespace {

constexpr double zero_tol = 1e-9;

struct PairwiseCovariance {
  RealMat sigma;
  RealMat mean_j;
  RealMat mean_k;
  RealMat p2;
};

Result<void> validate_alpha_inputs(IntMatView ratings, const RealVec& values) {
  const int n = static_cast<int>(ratings.rows());
  const int R = static_cast<int>(ratings.cols());
  const int C = static_cast<int>(values.size());
  if (n < 1 || R < 2 || C < 1) return std::unexpected(Error::invalid_argument);
  for (Eigen::Index k = 0; k < values.size(); ++k) {
    if (!std::isfinite(values(k))) return std::unexpected(Error::invalid_argument);
  }
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
  return {};
}

Result<PairwiseCovariance> pairwise_covariance(IntMatView ratings, const RealVec& values) {
  const int n = static_cast<int>(ratings.rows());
  const int R = static_cast<int>(ratings.cols());
  PairwiseCovariance out;
  out.sigma = RealMat::Zero(R, R);
  out.mean_j = RealMat::Zero(R, R);
  out.mean_k = RealMat::Zero(R, R);
  out.p2 = RealMat::Zero(R, R);

  for (int j = 0; j < R; ++j) {
    for (int k = j; k < R; ++k) {
      int count = 0;
      double sum_j = 0.0;
      double sum_k = 0.0;
      for (int i = 0; i < n; ++i) {
        const int xj = ratings(i, j);
        const int xk = ratings(i, k);
        if (xj == na_code || xk == na_code) continue;
        sum_j += values(xj);
        sum_k += values(xk);
        ++count;
      }
      if (count <= 0) return std::unexpected(Error::invalid_argument);
      const double muj = sum_j / static_cast<double>(count);
      const double muk = sum_k / static_cast<double>(count);
      double ss = 0.0;
      for (int i = 0; i < n; ++i) {
        const int xj = ratings(i, j);
        const int xk = ratings(i, k);
        if (xj == na_code || xk == na_code) continue;
        ss += (values(xj) - muj) * (values(xk) - muk);
      }
      const double s = ss / static_cast<double>(count);
      const double p = static_cast<double>(count) / static_cast<double>(n);
      out.sigma(j, k) = s;
      out.sigma(k, j) = s;
      out.mean_j(j, k) = muj;
      out.mean_k(j, k) = muk;
      out.mean_j(k, j) = muk;
      out.mean_k(k, j) = muj;
      out.p2(j, k) = p;
      out.p2(k, j) = p;
    }
  }
  return out;
}

double coefficient_alpha(double t1, double t2, int R) {
  if (std::abs(t1) <= zero_tol) return std::numeric_limits<double>::quiet_NaN();
  const double factor = static_cast<double>(R) / static_cast<double>(R - 1);
  return factor * (1.0 - t2 / t1);
}

RealMat alpha_gradient(double t1, double t2, int R) {
  RealMat g = RealMat::Zero(1, 2);
  if (std::abs(t1) <= zero_tol) return g;
  const double factor = static_cast<double>(R) / static_cast<double>(R - 1);
  g(0, 0) = factor * t2 / (t1 * t1);
  g(0, 1) = -factor / t1;
  return g;
}

}  // namespace

Result<Estimation> estimate_alpha_available(IntMatView ratings, const RealVec& values) {
  auto valid = validate_alpha_inputs(ratings, values);
  if (!valid) return std::unexpected(valid.error());

  const int n = static_cast<int>(ratings.rows());
  const int R = static_cast<int>(ratings.cols());
  auto cov = pairwise_covariance(ratings, values);
  if (!cov) return std::unexpected(cov.error());

  const double t1 = cov->sigma.sum();
  const double t2 = cov->sigma.diagonal().sum();

  RealVec estimates(1);
  estimates(0) = coefficient_alpha(t1, t2, R);

  RealMat phi = RealMat::Zero(n, 2);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < R; ++j) {
      for (int k = j; k < R; ++k) {
        const int xj = ratings(i, j);
        const int xk = ratings(i, k);
        if (xj == na_code || xk == na_code) continue;
        const double centered_product =
            (values(xj) - cov->mean_j(j, k)) * (values(xk) - cov->mean_k(j, k));
        const double term =
            (centered_product - cov->sigma(j, k)) / cov->p2(j, k);
        phi(i, 0) += (j == k) ? term : 2.0 * term;
        if (j == k) phi(i, 1) += term;
      }
    }
  }
  const RealMat gamma = (phi.transpose() * phi) / static_cast<double>(n);
  const RealMat J = alpha_gradient(t1, t2, R);
  RealMat vcov = (J * gamma * J.transpose()) / static_cast<double>(n);
  RealMat psi = build_psi_from_phi(phi, J);

  return Estimation{std::move(estimates), std::move(vcov), std::move(psi)};
}

}  // namespace misskappa
