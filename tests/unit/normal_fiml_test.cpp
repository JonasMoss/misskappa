#include "doctest.h"

#include "misskappa/estimate.hpp"

#include <cmath>
#include <limits>

using misskappa::EmOptions;
using misskappa::RealMat;
using misskappa::RealVec;
namespace ms = misskappa;

namespace {

constexpr double tol = 1e-8;
const double na_d = std::numeric_limits<double>::quiet_NaN();

struct Moments {
  RealVec mu;
  RealMat sigma;
};

Moments complete_moments(const RealMat& x) {
  Moments out;
  out.mu = x.colwise().mean().transpose();
  out.sigma = RealMat::Zero(x.cols(), x.cols());
  for (Eigen::Index i = 0; i < x.rows(); ++i) {
    const RealVec centered = x.row(i).transpose() - out.mu;
    out.sigma.noalias() += centered * centered.transpose();
  }
  out.sigma /= static_cast<double>(x.rows());
  return out;
}

double alpha_from_sigma(const RealMat& sigma) {
  const int p = static_cast<int>(sigma.cols());
  return (static_cast<double>(p) / static_cast<double>(p - 1)) *
         (1.0 - sigma.diagonal().sum() / sigma.sum());
}

RealVec scalar_quadratic_from_moments(const RealVec& mu, const RealMat& sigma) {
  const int R = static_cast<int>(mu.size());
  const double t1 = sigma.sum();
  const double t2 = sigma.diagonal().sum();
  const double t3 = (mu.array() - mu.mean()).square().sum();
  RealVec out(2);
  out(0) = (t1 - t2) / (static_cast<double>(R - 1) * t2 + static_cast<double>(R) * t3);
  out(1) = (t1 - t2 - t3) / (static_cast<double>(R - 1) * (t2 + t3));
  return out;
}

int vcol(int rater, int feature, int features) {
  return rater * features + feature;
}

RealVec vector_quadratic_from_moments(
    const RealVec& mu, const RealMat& sigma, int R, int features, const RealMat& W) {
  double tB = 0.0;
  double tT = 0.0;
  double tG = 0.0;
  for (int r = 0; r < R; ++r) {
    for (int s = 0; s < R; ++s) {
      const double p_rs = (r == s ? 1.0 : 0.0) - 1.0 / static_cast<double>(R);
      for (int a = 0; a < features; ++a) {
        for (int b = 0; b < features; ++b) {
          const int ca = vcol(r, a, features);
          const int cb = vcol(s, b, features);
          tB += W(a, b) * sigma(ca, cb);
          if (r == s) tT += W(a, b) * sigma(ca, cb);
          tG += p_rs * W(a, b) * mu(ca) * mu(cb);
        }
      }
    }
  }
  RealVec out(2);
  out(0) = (tB - tT) / (static_cast<double>(R - 1) * tT + static_cast<double>(R) * tG);
  out(1) = (tB - tT - tG) / (static_cast<double>(R - 1) * (tT + tG));
  return out;
}

void check_influence_contract(const ms::NormalFimlEstimation& e, int n, int k) {
  REQUIRE(e.fit.psi.rows() == n);
  REQUIRE(e.fit.psi.cols() == k);
  const RealMat via_psi =
      (e.fit.psi.transpose() * e.fit.psi) / std::pow(static_cast<double>(n), 2);
  CHECK((via_psi - e.fit.vcov).cwiseAbs().maxCoeff() < 1e-10);
}

RealMat complete_items() {
  RealMat x(8, 4);
  x << 1.0, 1.2, 0.9, 1.1,
       2.0, 2.2, 1.8, 2.1,
       3.0, 2.8, 3.2, 3.1,
       2.0, 2.4, 2.1, 2.3,
       1.0, 1.5, 1.2, 1.3,
       3.0, 2.9, 2.7, 2.8,
       2.0, 2.0, 2.1, 2.2,
       4.0, 3.8, 4.2, 4.1;
  return x;
}

}  // namespace

TEST_CASE("normal FIML alpha: complete data matches MLE covariance alpha") {
  const RealMat x = complete_items();
  const Moments m = complete_moments(x);

  auto r = ms::estimate_alpha_normal_fiml(x, EmOptions{});
  REQUIRE(r.has_value());
  CHECK(r->converged);
  CHECK(std::abs(r->fit.estimates(0) - alpha_from_sigma(m.sigma)) < tol);
  CHECK((r->mu - m.mu).cwiseAbs().maxCoeff() < tol);
  CHECK((r->sigma - m.sigma).cwiseAbs().maxCoeff() < tol);
  check_influence_contract(*r, static_cast<int>(x.rows()), 1);
}

TEST_CASE("normal FIML quadratic kappa: complete data matches moment formula") {
  const RealMat x = complete_items();
  const Moments m = complete_moments(x);

  auto r = ms::estimate_quadratic_normal_fiml(x, EmOptions{});
  REQUIRE(r.has_value());
  CHECK(r->converged);
  CHECK((r->fit.estimates - scalar_quadratic_from_moments(m.mu, m.sigma))
            .cwiseAbs()
            .maxCoeff() < tol);
  check_influence_contract(*r, static_cast<int>(x.rows()), 2);
}

TEST_CASE("normal FIML quadratic kappa: missing data returns finite IF covariance") {
  RealMat x = complete_items();
  x(1, 1) = na_d;
  x(3, 0) = na_d;
  x(5, 2) = na_d;

  auto r = ms::estimate_quadratic_normal_fiml(x, EmOptions{});
  REQUIRE(r.has_value());
  CHECK(r->fit.estimates.allFinite());
  CHECK(r->fit.vcov.allFinite());
  check_influence_contract(*r, static_cast<int>(x.rows()), 2);
}

TEST_CASE("normal FIML vector quadratic: complete data matches vector moment formula") {
  RealMat x(7, 6);
  x << 0.0, 1.0, 0.2, 1.1, 0.4, 1.2,
       1.0, 0.0, 1.1, 0.2, 1.2, 0.5,
       0.5, 0.4, 0.8, 0.6, 0.7, 0.9,
       1.4, 1.1, 1.6, 1.3, 1.8, 1.4,
       0.3, 0.2, 0.4, 0.3, 0.6, 0.4,
       1.2, 1.5, 1.1, 1.7, 1.4, 1.8,
       0.8, 0.7, 0.9, 1.0, 1.0, 1.1;
  RealMat W(2, 2);
  W << 1.0, 0.25,
       0.25, 1.5;
  const Moments m = complete_moments(x);

  auto r = ms::estimate_vector_quadratic_normal_fiml(x, 2, W, EmOptions{});
  REQUIRE(r.has_value());
  CHECK(r->converged);
  CHECK((r->fit.estimates - vector_quadratic_from_moments(m.mu, m.sigma, 3, 2, W))
            .cwiseAbs()
            .maxCoeff() < tol);
  check_influence_contract(*r, static_cast<int>(x.rows()), 2);
}
