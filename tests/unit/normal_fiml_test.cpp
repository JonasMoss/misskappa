#include "doctest.h"

#include "misskappa/estimate.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>

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

Moments complete_equal_mean_moments(const RealMat& x) {
  Moments out;
  out.mu = RealVec::Constant(x.cols(), x.sum() / static_cast<double>(x.size()));
  out.sigma = (x.transpose() * x) / static_cast<double>(x.rows()) -
              out.mu * out.mu.transpose();
  return out;
}

double alpha_from_sigma(const RealMat& sigma) {
  const int p = static_cast<int>(sigma.cols());
  return (static_cast<double>(p) / static_cast<double>(p - 1)) *
         (1.0 - sigma.diagonal().sum() / sigma.sum());
}

double continuous_bp_chance(const RealMat& x) {
  double lo = std::numeric_limits<double>::infinity();
  double hi = -std::numeric_limits<double>::infinity();
  for (Eigen::Index i = 0; i < x.rows(); ++i) {
    for (Eigen::Index j = 0; j < x.cols(); ++j) {
      if (!std::isfinite(x(i, j))) continue;
      lo = std::min(lo, x(i, j));
      hi = std::max(hi, x(i, j));
    }
  }
  return (hi - lo) * (hi - lo) / 6.0;
}

RealVec scalar_quadratic_from_moments(
    const RealVec& mu, const RealMat& sigma, double bp_chance) {
  const int R = static_cast<int>(mu.size());
  const double t1 = sigma.sum();
  const double t2 = sigma.diagonal().sum();
  const double t3 = (mu.array() - mu.mean()).square().sum();
  const double d_obs =
      (2.0 / static_cast<double>(R - 1)) * (t2 + t3 - t1 / static_cast<double>(R));
  RealVec out(3);
  out(0) = (t1 - t2) / (static_cast<double>(R - 1) * t2 + static_cast<double>(R) * t3);
  out(1) = (t1 - t2 - t3) / (static_cast<double>(R - 1) * (t2 + t3));
  out(2) = 1.0 - d_obs / bp_chance;
  return out;
}

RealMat variance_component_projection(const RealMat& sigma) {
  const int R = static_cast<int>(sigma.rows());
  const double trace = sigma.diagonal().sum();
  const double total = sigma.sum();
  const double contrast = (trace - total / static_cast<double>(R)) /
                          static_cast<double>(R - 1);
  const double mean = total / static_cast<double>(R);
  double e = 0.0;
  double psi = 0.0;
  if (mean >= contrast) {
    e = contrast;
    psi = (mean - contrast) / static_cast<double>(R);
  } else {
    e = trace / static_cast<double>(R);
    psi = 0.0;
  }
  RealMat out = RealMat::Constant(R, R, psi);
  out.diagonal().array() += e;
  return out;
}

int vcol(int rater, int feature, int features) {
  return rater * features + feature;
}

// Brennan-Prediger uniform chance for vector ratings: sum_a W_aa * range_a^2/6
// over the observed score range of feature a across raters.
double vector_bp_chance_ref(const RealMat& x, int R, int features, const RealMat& W) {
  double chance = 0.0;
  for (int a = 0; a < features; ++a) {
    double lo = std::numeric_limits<double>::infinity();
    double hi = -std::numeric_limits<double>::infinity();
    for (int r = 0; r < R; ++r) {
      for (int i = 0; i < x.rows(); ++i) {
        const double v = x(i, vcol(r, a, features));
        if (!std::isfinite(v)) continue;
        lo = std::min(lo, v);
        hi = std::max(hi, v);
      }
    }
    chance += W(a, a) * (hi - lo) * (hi - lo) / 6.0;
  }
  return chance;
}

RealVec vector_quadratic_from_moments(
    const RealVec& mu, const RealMat& sigma, int R, int features, const RealMat& W,
    double bp_chance) {
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
  const double d_obs =
      (2.0 / static_cast<double>(R - 1)) * (tT + tG - tB / static_cast<double>(R));
  RealVec out(3);
  out(0) = (tB - tT) / (static_cast<double>(R - 1) * tT + static_cast<double>(R) * tG);
  out(1) = (tB - tT - tG) / (static_cast<double>(R - 1) * (tT + tG));
  out(2) = 1.0 - d_obs / bp_chance;
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

RealMat benchmark_items(int n, int p, double missing_prob, unsigned seed) {
  std::mt19937 gen(seed);
  std::normal_distribution<double> normal(0.0, 1.0);
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  RealMat x(n, p);
  for (int i = 0; i < n; ++i) {
    const double factor = normal(gen);
    int observed = 0;
    for (int j = 0; j < p; ++j) {
      x(i, j) = 1.0 + 0.8 * factor + 0.6 * normal(gen);
      if (uniform(gen) < missing_prob) {
        x(i, j) = na_d;
      } else {
        ++observed;
      }
    }
    if (observed == 0) x(i, 0) = 1.0 + normal(gen);
  }
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

TEST_CASE("normal FIML alpha: variance components projects complete covariance") {
  const RealMat x = complete_items();
  const Moments m = complete_moments(x);
  const RealMat sigma_vc = variance_component_projection(m.sigma);
  EmOptions opts{};
  opts.covariance = ms::NormalCovarianceModel::variance_components;

  auto r = ms::estimate_alpha_normal_fiml(x, opts);
  REQUIRE(r.has_value());
  CHECK(r->converged);
  CHECK((r->mu - m.mu).cwiseAbs().maxCoeff() < tol);
  CHECK((r->sigma - sigma_vc).cwiseAbs().maxCoeff() < tol);
  CHECK(std::abs(r->fit.estimates(0) - alpha_from_sigma(sigma_vc)) < tol);
  check_influence_contract(*r, static_cast<int>(x.rows()), 1);
}

TEST_CASE("normal FIML alpha: variance components allows connected patterns") {
  RealMat x(8, 3);
  x << 1.0, 1.1, na_d,
       2.0, 2.1, na_d,
       3.0, 3.2, na_d,
       na_d, 1.0, 1.2,
       na_d, 2.0, 2.3,
       na_d, 3.1, 3.4,
       1.5, na_d, na_d,
       na_d, na_d, 2.6;

  auto saturated = ms::estimate_alpha_normal_fiml(x, EmOptions{});
  REQUIRE(!saturated.has_value());
  CHECK(saturated.error() == ms::Error::not_identified);

  EmOptions opts{};
  opts.covariance = ms::NormalCovarianceModel::variance_components;
  auto structured = ms::estimate_alpha_normal_fiml(x, opts);
  REQUIRE(structured.has_value());
  CHECK(structured->fit.estimates.allFinite());
  CHECK(structured->fit.vcov.allFinite());
  check_influence_contract(*structured, static_cast<int>(x.rows()), 1);
}

TEST_CASE("normal FIML quadratic kappa: complete data matches moment formula") {
  const RealMat x = complete_items();
  const Moments m = complete_moments(x);
  const double bp_chance = continuous_bp_chance(x);

  auto r = ms::estimate_quadratic_normal_fiml(x, RealVec(), EmOptions{});
  REQUIRE(r.has_value());
  CHECK(r->converged);
  CHECK((r->fit.estimates - scalar_quadratic_from_moments(m.mu, m.sigma, bp_chance))
            .cwiseAbs()
            .maxCoeff() < tol);
  check_influence_contract(*r, static_cast<int>(x.rows()), 3);
}

TEST_CASE("normal FIML quadratic kappa: complete data matches empirical BP") {
  const RealMat x = complete_items();
  RealVec values(5);
  values << 0.0, 1.0, 2.0, 3.0, 4.0;

  auto empirical = ms::estimate_quadratic(x, values);
  auto fiml = ms::estimate_quadratic_normal_fiml(x, values, EmOptions{});
  EmOptions opts{};
  opts.covariance = ms::NormalCovarianceModel::variance_components;
  auto structured = ms::estimate_quadratic_normal_fiml(x, values, opts);
  REQUIRE(empirical.has_value());
  REQUIRE(fiml.has_value());
  REQUIRE(structured.has_value());
  CHECK((fiml->fit.estimates - empirical->estimates).cwiseAbs().maxCoeff() < tol);
  CHECK((structured->fit.estimates - empirical->estimates).cwiseAbs().maxCoeff() < tol);
  check_influence_contract(*fiml, static_cast<int>(x.rows()), 3);
  check_influence_contract(*structured, static_cast<int>(x.rows()), 3);
}

TEST_CASE("normal FIML quadratic kappa: variance components projects complete covariance") {
  const RealMat x = complete_items();
  const Moments m = complete_moments(x);
  const RealMat sigma_vc = variance_component_projection(m.sigma);
  const double bp_chance = continuous_bp_chance(x);
  EmOptions opts{};
  opts.covariance = ms::NormalCovarianceModel::variance_components;

  auto r = ms::estimate_quadratic_normal_fiml(x, RealVec(), opts);
  REQUIRE(r.has_value());
  CHECK(r->converged);
  CHECK((r->mu - m.mu).cwiseAbs().maxCoeff() < tol);
  CHECK((r->sigma - sigma_vc).cwiseAbs().maxCoeff() < tol);
  CHECK((r->fit.estimates - scalar_quadratic_from_moments(m.mu, sigma_vc, bp_chance))
            .cwiseAbs()
            .maxCoeff() < tol);
  check_influence_contract(*r, static_cast<int>(x.rows()), 3);
}

TEST_CASE("normal FIML quadratic kappa: equal mean variance components are exchangeable") {
  const RealMat x = complete_items();
  const Moments m = complete_equal_mean_moments(x);
  const RealMat sigma_vc = variance_component_projection(m.sigma);
  const double bp_chance = continuous_bp_chance(x);
  EmOptions opts{};
  opts.covariance = ms::NormalCovarianceModel::variance_components;
  opts.mean = ms::NormalMeanModel::equal;

  auto r = ms::estimate_quadratic_normal_fiml(x, RealVec(), opts);
  REQUIRE(r.has_value());
  CHECK(r->converged);
  CHECK((r->mu - m.mu).cwiseAbs().maxCoeff() < tol);
  CHECK((r->sigma - sigma_vc).cwiseAbs().maxCoeff() < tol);
  CHECK((r->fit.estimates - scalar_quadratic_from_moments(m.mu, sigma_vc, bp_chance))
            .cwiseAbs()
            .maxCoeff() < tol);
  CHECK(std::abs(r->fit.estimates(0) - r->fit.estimates(1)) < tol);
  check_influence_contract(*r, static_cast<int>(x.rows()), 3);
}

TEST_CASE("normal FIML quadratic kappa: equal mean requires variance components") {
  const RealMat x = complete_items();
  EmOptions opts{};
  opts.mean = ms::NormalMeanModel::equal;

  auto r = ms::estimate_quadratic_normal_fiml(x, RealVec(), opts);
  REQUIRE(!r.has_value());
  CHECK(r.error() == ms::Error::invalid_argument);
}

TEST_CASE("normal FIML quadratic kappa: missing data returns finite IF covariance") {
  RealMat x = complete_items();
  x(1, 1) = na_d;
  x(3, 0) = na_d;
  x(5, 2) = na_d;

  auto r = ms::estimate_quadratic_normal_fiml(x, RealVec(), EmOptions{});
  REQUIRE(r.has_value());
  CHECK(r->fit.estimates.allFinite());
  CHECK(r->fit.vcov.allFinite());
  check_influence_contract(*r, static_cast<int>(x.rows()), 3);
}

TEST_CASE("normal FIML quadratic kappa: variance components allows connected patterns") {
  RealMat x(8, 3);
  x << 1.0, 1.1, na_d,
       2.0, 2.1, na_d,
       3.0, 3.2, na_d,
       na_d, 1.0, 1.2,
       na_d, 2.0, 2.3,
       na_d, 3.1, 3.4,
       1.5, na_d, na_d,
       na_d, na_d, 2.6;

  auto saturated = ms::estimate_quadratic_normal_fiml(x, RealVec(), EmOptions{});
  REQUIRE(!saturated.has_value());
  CHECK(saturated.error() == ms::Error::not_identified);

  EmOptions opts{};
  opts.covariance = ms::NormalCovarianceModel::variance_components;
  auto structured = ms::estimate_quadratic_normal_fiml(x, RealVec(), opts);
  REQUIRE(structured.has_value());
  CHECK(structured->fit.estimates.allFinite());
  CHECK(structured->fit.vcov.allFinite());
  check_influence_contract(*structured, static_cast<int>(x.rows()), 3);
}

TEST_CASE("normal FIML quadratic kappa: missing-data estimates and SEs are pinned") {
  const RealMat x = benchmark_items(500, 6, 0.2, 548);

  auto r = ms::estimate_quadratic_normal_fiml(x, RealVec(), EmOptions{});
  REQUIRE(r.has_value());
  CHECK(r->converged);
  CHECK(std::abs(r->fit.estimates(0) - 0.627565184290) < 1e-11);
  CHECK(std::abs(r->fit.estimates(1) - 0.627513328957) < 1e-11);
  CHECK(std::isfinite(r->fit.estimates(2)));
  CHECK(std::abs(std::sqrt(r->fit.vcov(0, 0)) - 0.018896194796) < 1e-8);
  CHECK(std::abs(std::sqrt(r->fit.vcov(1, 1)) - 0.018902343240) < 1e-8);
  check_influence_contract(*r, static_cast<int>(x.rows()), 3);
}

TEST_CASE("normal FIML quadratic kappa: singular covariance returns numerical error") {
  RealMat x(8, 3);
  x << 1.0, 1.0, 0.4,
       2.0, 2.0, 1.3,
       3.0, 3.0, 2.5,
       4.0, 4.0, 3.6,
       5.0, 5.0, 4.2,
       6.0, 6.0, 5.4,
       7.0, 7.0, 6.7,
       8.0, 8.0, 7.5;

  auto r = ms::estimate_quadratic_normal_fiml(x, RealVec(), EmOptions{});
  REQUIRE(!r.has_value());
  CHECK(r.error() == ms::Error::numerical_error);
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

  const double bp = vector_bp_chance_ref(x, 3, 2, W);
  auto r = ms::estimate_vector_quadratic_normal_fiml(x, 2, W, EmOptions{});
  REQUIRE(r.has_value());
  CHECK(r->converged);
  CHECK((r->fit.estimates - vector_quadratic_from_moments(m.mu, m.sigma, 3, 2, W, bp))
            .cwiseAbs()
            .maxCoeff() < tol);
  check_influence_contract(*r, static_cast<int>(x.rows()), 3);
}

TEST_CASE("normal FIML vector quadratic: block compound symmetry is fitted") {
  RealMat x(7, 6);
  x << 0.0, 1.0, 0.2, 1.1, 0.4, 1.2,
       1.0, 0.0, 1.1, 0.2, 1.2, 0.5,
       0.5, 0.4, 0.8, 0.6, 0.7, 0.9,
       1.4, 1.1, 1.6, 1.3, 1.8, 1.4,
       0.3, 0.2, 0.4, 0.3, 0.6, 0.4,
       1.2, 1.5, 1.1, 1.7, 1.4, 1.8,
       0.8, 0.7, 0.9, 1.0, 1.0, 1.1;
  RealMat W = RealMat::Identity(2, 2);
  const int R = 3, d = 2;
  EmOptions opts{};
  opts.covariance = ms::NormalCovarianceModel::variance_components;

  auto sat = ms::estimate_vector_quadratic_normal_fiml(x, d, W, EmOptions{});
  auto r = ms::estimate_vector_quadratic_normal_fiml(x, d, W, opts);
  REQUIRE(r.has_value());
  CHECK(r->converged);
  REQUIRE(sat.has_value());

  // fitted Sigma is block compound symmetric: every off-diagonal d x d block is
  // a common Psi, every diagonal block is Psi + E.
  const RealMat& S = r->sigma;
  const RealMat Psi = S.block(0, d, d, d);
  const RealMat E = S.block(0, 0, d, d) - Psi;
  for (int a = 0; a < R; ++a) {
    for (int b = 0; b < R; ++b) {
      const RealMat ref = (a == b) ? (Psi + E).eval() : Psi;
      CHECK((S.block(a * d, b * d, d, d) - ref).cwiseAbs().maxCoeff() < 1e-9);
    }
  }
  // genuine constraint: differs from the saturated fit
  CHECK((r->sigma - sat->sigma).cwiseAbs().maxCoeff() > 1e-6);
  check_influence_contract(*r, static_cast<int>(x.rows()), 3);
}

TEST_CASE("normal FIML vector quadratic: block CS reduces to scalar CS at d=1") {
  // d = 1 (one feature) must reproduce the scalar quadratic CS fit exactly.
  RealMat x(8, 4);
  x << 1.0, 1.2, 0.9, 1.1,
       2.0, 2.2, 1.8, 2.1,
       3.0, 2.8, 3.2, 3.1,
       2.0, 2.4, 2.1, 2.3,
       1.0, 1.5, 1.2, 1.3,
       3.0, 2.9, 2.7, 2.8,
       2.5, 2.3, 2.6, 2.4,
       1.5, 1.7, 1.4, 1.6;
  RealMat W = RealMat::Identity(1, 1);
  RealVec values(0);
  EmOptions opts{};
  opts.covariance = ms::NormalCovarianceModel::variance_components;

  auto vec = ms::estimate_vector_quadratic_normal_fiml(x, 1, W, opts);
  auto scl = ms::estimate_quadratic_normal_fiml(x, values, opts);
  REQUIRE(vec.has_value());
  REQUIRE(scl.has_value());
  CHECK((vec->fit.estimates - scl->fit.estimates).cwiseAbs().maxCoeff() < 1e-9);
  CHECK((vec->fit.vcov - scl->fit.vcov).cwiseAbs().maxCoeff() < 1e-9);
  CHECK((vec->sigma - scl->sigma).cwiseAbs().maxCoeff() < 1e-9);
}
