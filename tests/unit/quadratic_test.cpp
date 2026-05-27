#include "doctest.h"

#include "misskappa/estimate.hpp"
#include "misskappa/loss.hpp"

#include <Eigen/Eigenvalues>
#include <cmath>
#include <limits>

using misskappa::IntMat;
using misskappa::RealMat;
using misskappa::RealVec;
namespace ms = misskappa;

namespace {

constexpr double tol = 1e-9;
const double na_d = std::numeric_limits<double>::quiet_NaN();

struct QuadraticMoments {
  RealVec mu;
  RealMat sigma;
  RealVec p1;
  RealMat p2;
  Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> mask;
};

RealMat twelve_subject_3rater_quadratic() {
  RealMat x(12, 3);
  x <<
    1, 2, 1,
    2, 2, 3,
    3, 3, 3,
    4, 4, 5,
    5, 5, 5,
    1, 1, 2,
    3, 3, 4,
    4, 5, 4,
    2, 3, 2,
    5, 4, 5,
    1, 2, na_d,
    3, na_d, 4;
  return x;
}

QuadraticMoments moment_bundle(const RealMat& x) {
  const int n = static_cast<int>(x.rows());
  const int R = static_cast<int>(x.cols());
  QuadraticMoments out;
  out.mask.resize(n, R);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < R; ++j) out.mask(i, j) = std::isfinite(x(i, j)) ? 1 : 0;
  }
  out.p1 = out.mask.cast<double>().colwise().sum().transpose() / static_cast<double>(n);
  out.p2 = (out.mask.cast<double>().transpose() * out.mask.cast<double>()) /
           static_cast<double>(n);
  out.mu = RealVec::Zero(R);
  for (int j = 0; j < R; ++j) {
    double sum = 0.0;
    int count = 0;
    for (int i = 0; i < n; ++i) {
      if (!out.mask(i, j)) continue;
      sum += x(i, j);
      ++count;
    }
    out.mu(j) = sum / count;
  }
  out.sigma = RealMat::Zero(R, R);
  for (int j = 0; j < R; ++j) {
    for (int k = j; k < R; ++k) {
      double sum_j = 0.0;
      double sum_k = 0.0;
      int count = 0;
      for (int i = 0; i < n; ++i) {
        if (!out.mask(i, j) || !out.mask(i, k)) continue;
        sum_j += x(i, j);
        sum_k += x(i, k);
        ++count;
      }
      const double mu_j = sum_j / count;
      const double mu_k = sum_k / count;
      double sum = 0.0;
      for (int i = 0; i < n; ++i) {
        if (!out.mask(i, j) || !out.mask(i, k)) continue;
        sum += (x(i, j) - mu_j) * (x(i, k) - mu_k);
      }
      out.sigma(j, k) = sum / count;
      out.sigma(k, j) = out.sigma(j, k);
    }
  }
  return out;
}

double p4_for_indices(const QuadraticMoments& m, int a, int b, int c, int d) {
  int count = 0;
  for (int row = 0; row < m.mask.rows(); ++row) {
    if (m.mask(row, a) && m.mask(row, b) && m.mask(row, c) && m.mask(row, d)) ++count;
  }
  return static_cast<double>(count) / static_cast<double>(m.mask.rows());
}

RealMat brute_elliptical_summary_psi(const QuadraticMoments& m, double beta) {
  const int R = static_cast<int>(m.mu.size());
  RealMat psi = RealMat::Zero(3, 3);
  for (int a = 0; a < R; ++a) {
    for (int b = 0; b < R; ++b) {
      for (int c = 0; c < R; ++c) {
        for (int d = 0; d < R; ++d) {
          const double cov_y =
              beta * (m.sigma(a, c) * m.sigma(b, d) +
                      m.sigma(a, d) * m.sigma(b, c)) +
              (beta - 1.0) * m.sigma(a, b) * m.sigma(c, d);
          const double miss =
              p4_for_indices(m, a, b, c, d) / (m.p2(a, b) * m.p2(c, d));
          const double gamma = miss * cov_y;
          psi(0, 0) += gamma;
          if (c == d) psi(0, 1) += gamma;
          if (a == b) psi(1, 0) += gamma;
          if (a == b && c == d) psi(1, 1) += gamma;
        }
      }
    }
  }
  RealMat gamma_mu = RealMat::Zero(R, R);
  for (int j = 0; j < R; ++j) {
    for (int k = 0; k < R; ++k) {
      gamma_mu(j, k) = m.p2(j, k) * m.sigma(j, k) / (m.p1(j) * m.p1(k));
    }
  }
  const RealVec hmu = (m.mu.array() - m.mu.mean()).matrix();
  psi(2, 2) = 4.0 * (hmu.transpose() * gamma_mu * hmu)(0);
  return psi;
}

RealMat vcov_from_summary_psi(const QuadraticMoments& m, const RealVec& values,
                              const RealMat& psi) {
  const int C = static_cast<int>(values.size());
  const int R = static_cast<int>(m.mu.size());
  const int n = static_cast<int>(m.mask.rows());
  const double c1 = (2.0 / (C * C)) *
                    (C * values.squaredNorm() - std::pow(values.sum(), 2));
  const double t1 = m.sigma.sum();
  const double t2 = m.sigma.diagonal().sum();
  const double t3 = (m.mu.array() - m.mu.mean()).matrix().squaredNorm();
  const double NF = t1 - t2 - t3;
  const double DF = (R - 1.0) * (t2 + t3);
  const double NC = t1 - t2;
  const double DC = (R - 1.0) * t2 + static_cast<double>(R) * t3;

  RealVec gF = RealVec::Zero(3);
  RealVec gC = RealVec::Zero(3);
  RealVec gBP = RealVec::Zero(3);
  gF(0) = 1.0 / DF;
  gF(1) = (-DF - NF * (R - 1.0)) / (DF * DF);
  gF(2) = gF(1);
  gC(0) = 1.0 / DC;
  gC(1) = (-DC - NC * (R - 1.0)) / (DC * DC);
  gC(2) = (-NC * static_cast<double>(R)) / (DC * DC);
  const double bp_scale = -(2.0 / (c1 * (R - 1.0)));
  gBP(0) = bp_scale * (-1.0 / R);
  gBP(1) = bp_scale;
  gBP(2) = bp_scale;

  RealMat G(3, 3);
  G.col(0) = gF;
  G.col(1) = gC;
  G.col(2) = gBP;
  const RealMat scaled = G.transpose() * psi * G;
  RealMat out(3, 3);
  const int perm[3] = {1, 0, 2};
  for (int r = 0; r < 3; ++r) {
    for (int c = 0; c < 3; ++c) out(r, c) = scaled(perm[r], perm[c]);
  }
  return out / static_cast<double>(n);
}

}  // namespace

TEST_CASE("estimate_quadratic raw: matches legacy on 12-subject 3-rater fixture") {
  RealMat x = twelve_subject_3rater_quadratic();
  RealVec v(5);
  v << 1.0, 2.0, 3.0, 4.0, 5.0;

  auto r = ms::estimate_quadratic(x, v);
  REQUIRE(r.has_value());
  REQUIRE(r->estimates.size() == 3);

  // Point estimates frozen against dev/legacy/misskappa
  // kappa_raw(method="quadratic", values=1:5).
  CHECK(std::abs(r->estimates(0) - 0.8374174174174175) < tol);  // Conger
  CHECK(std::abs(r->estimates(1) - 0.8344881466279811) < tol);  // Fleiss
  CHECK(std::abs(r->estimates(2) - 0.8508539944903584) < tol);  // BP

  CHECK(std::abs(r->vcov(0, 0) - 0.003337365350253323) < 1e-12);
  CHECK(std::abs(r->vcov(1, 1) - 0.0036216673879180795) < 1e-12);
  CHECK(std::abs(r->vcov(2, 2) - 0.0019460442258634203) < 1e-12);
}

TEST_CASE("estimate_quadratic raw: variance is symmetric and PSD") {
  RealMat x = twelve_subject_3rater_quadratic();
  RealVec v(5);
  v << 1.0, 2.0, 3.0, 4.0, 5.0;
  auto r = ms::estimate_quadratic(x, v);
  REQUIRE(r.has_value());
  const RealMat& V = r->vcov;
  CHECK((V - V.transpose()).cwiseAbs().maxCoeff() < 1e-10);
  Eigen::SelfAdjointEigenSolver<RealMat> es(V);
  REQUIRE(es.info() == Eigen::Success);
  CHECK(es.eigenvalues().minCoeff() > -1e-9);
}

TEST_CASE("estimate_quadratic raw: normal vcov matches missing-data Gamma contraction") {
  RealMat x(6, 3);
  x << 1, 2, 1,
       2, 2, 3,
       3, 3, 3,
       4, 4, 5,
       5, 5, 5,
       na_d, na_d, na_d;
  RealVec v(5);
  v << 1.0, 2.0, 3.0, 4.0, 5.0;

  ms::QuadraticOptions opts;
  opts.vcov = ms::QuadraticVcov::normal;
  auto r = ms::estimate_quadratic(x, v, opts);
  REQUIRE(r.has_value());

  const QuadraticMoments m = moment_bundle(x);
  const RealMat expected = vcov_from_summary_psi(m, v, brute_elliptical_summary_psi(m, 1.0));
  CHECK((r->vcov - expected).cwiseAbs().maxCoeff() < 1e-12);
}

TEST_CASE("estimate_quadratic raw: normal vcov matches complete-data Gamma contraction") {
  RealMat x(5, 3);
  x << 1, 2, 1,
       2, 2, 3,
       3, 3, 3,
       4, 4, 5,
       5, 5, 5;
  RealVec v(5);
  v << 1.0, 2.0, 3.0, 4.0, 5.0;

  ms::QuadraticOptions opts;
  opts.vcov = ms::QuadraticVcov::normal;
  auto r = ms::estimate_quadratic(x, v, opts);
  REQUIRE(r.has_value());

  const QuadraticMoments m = moment_bundle(x);
  const RealMat expected = vcov_from_summary_psi(m, v, brute_elliptical_summary_psi(m, 1.0));
  CHECK((r->vcov - expected).cwiseAbs().maxCoeff() < 1e-12);
}

TEST_CASE("estimate_quadratic raw: elliptical vcov matches fourth-moment contraction") {
  RealMat x(6, 3);
  x << 1, 2, 1,
       2, 2, 3,
       3, 3, 3,
       4, 4, 5,
       5, 5, 5,
       na_d, na_d, na_d;
  RealVec v(5);
  v << 1.0, 2.0, 3.0, 4.0, 5.0;

  ms::QuadraticOptions opts;
  opts.vcov = ms::QuadraticVcov::elliptical;
  opts.relative_kurtosis = 1.6;
  auto r = ms::estimate_quadratic(x, v, opts);
  REQUIRE(r.has_value());

  const QuadraticMoments m = moment_bundle(x);
  const RealMat expected = vcov_from_summary_psi(m, v, brute_elliptical_summary_psi(m, 1.6));
  CHECK((r->vcov - expected).cwiseAbs().maxCoeff() < 1e-12);
}

TEST_CASE("estimate_quadratic raw: perfect agreement -> kappa = 1") {
  RealMat x(5, 3);
  x << 1, 1, 1,
       2, 2, 2,
       3, 3, 3,
       4, 4, 4,
       5, 5, 5;
  RealVec v(5);
  v << 1.0, 2.0, 3.0, 4.0, 5.0;
  auto r = ms::estimate_quadratic(x, v);
  REQUIRE(r.has_value());
  CHECK(std::abs(r->estimates(0) - 1.0) < tol);
  CHECK(std::abs(r->estimates(1) - 1.0) < tol);
  CHECK(std::abs(r->estimates(2) - 1.0) < tol);
}

TEST_CASE("estimate_quadratic raw: singleton rows contribute to marginal moments") {
  RealMat x(6, 3);
  x << 1, 2, 3,
       2, 3, 4,
       3, 4, 5,
       5, na_d, na_d,
       na_d, 1, na_d,
       na_d, na_d, 5;
  RealVec v(5);
  v << 1.0, 2.0, 3.0, 4.0, 5.0;
  auto r = ms::estimate_quadratic(x, v);
  REQUIRE(r.has_value());

  CHECK(std::abs(r->estimates(0) - 0.29357798165137616) < tol);  // Conger
  CHECK(std::abs(r->estimates(1) - 0.18661971830985918) < tol);  // Fleiss
  CHECK(std::abs(r->estimates(2) - 0.19791666666666685) < tol);  // BP
  CHECK(std::abs(r->vcov(0, 0) - 0.014298806444653409) < 1e-12);
  CHECK(std::abs(r->vcov(1, 1) - 0.026515177274675391) < 1e-12);
  CHECK(std::abs(r->vcov(2, 2) - 0.053811909239969126) < 1e-12);
}

TEST_CASE("estimate_quadratic raw: too few overlapping pairs -> invalid_argument") {
  RealMat x(2, 3);
  x << 1, na_d, na_d,
       na_d, 2, na_d;
  RealVec v(5);
  v << 1.0, 2.0, 3.0, 4.0, 5.0;
  auto r = ms::estimate_quadratic(x, v);
  REQUIRE(!r.has_value());
  CHECK(r.error() == ms::Error::invalid_argument);
}

TEST_CASE("estimate_quadratic_counts: matches legacy on Fleiss 1971") {
  IntMat x(30, 5);
  x <<
    0, 0, 0, 6, 0,
    0, 3, 0, 0, 3,
    0, 1, 4, 0, 1,
    0, 0, 0, 0, 6,
    0, 3, 0, 3, 0,
    2, 0, 4, 0, 0,
    0, 0, 4, 0, 2,
    2, 0, 3, 1, 0,
    2, 0, 0, 4, 0,
    0, 0, 0, 0, 6,
    1, 0, 0, 5, 0,
    1, 1, 0, 4, 0,
    0, 3, 3, 0, 0,
    1, 0, 0, 5, 0,
    0, 2, 0, 3, 1,
    0, 0, 5, 0, 1,
    3, 0, 0, 1, 2,
    5, 1, 0, 0, 0,
    0, 2, 0, 4, 0,
    1, 0, 2, 0, 3,
    0, 0, 0, 0, 6,
    0, 1, 0, 5, 0,
    0, 2, 0, 1, 3,
    2, 0, 0, 4, 0,
    1, 0, 0, 4, 1,
    0, 5, 0, 1, 0,
    4, 0, 0, 0, 2,
    0, 2, 0, 4, 0,
    1, 0, 5, 0, 0,
    0, 0, 0, 0, 6;
  RealVec v(5);
  v << 1.0, 2.0, 3.0, 4.0, 5.0;

  auto r = ms::estimate_quadratic_counts(x, v, /*R_total=*/6);
  REQUIRE(r.has_value());
  REQUIRE(r->estimates.size() == 2);

  // Frozen against dev/legacy/misskappa kappa_counts(method="quadratic").
  CHECK(std::abs(r->estimates(0) - 0.2840722495894910) < tol);  // Fleiss
  CHECK(std::abs(r->estimates(1) - 0.3338888888888888) < tol);  // BP

  CHECK(std::abs(r->vcov(0, 0) - 0.011948832185177609) < 1e-11);
  CHECK(std::abs(r->vcov(1, 1) - 0.010378693415637774) < 1e-11);
  CHECK(std::abs(r->vcov(0, 1) - 0.008380603825142614) < 1e-11);
}

TEST_CASE("estimate_quadratic_counts: dimension mismatch -> dimension_mismatch") {
  IntMat x(3, 3);
  x << 2, 1, 1,
       3, 0, 1,
       0, 3, 1;
  RealVec v(5);  // mismatch
  v << 1.0, 2.0, 3.0, 4.0, 5.0;
  auto r = ms::estimate_quadratic_counts(x, v, 4);
  REQUIRE(!r.has_value());
  CHECK(r.error() == ms::Error::dimension_mismatch);
}
