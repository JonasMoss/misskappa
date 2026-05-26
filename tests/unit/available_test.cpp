#include "doctest.h"

#include "misskappa/estimate.hpp"
#include "misskappa/loss.hpp"

#include <Eigen/Eigenvalues>
#include <cmath>

using misskappa::Estimation;
using misskappa::IntMat;
using misskappa::na_code;
using misskappa::RealMat;
using misskappa::RealVec;
namespace ms = misskappa;

namespace {

constexpr double tol = 1e-10;

// Build a 2-rater categorical matrix from a flat list (row-major).
IntMat ratings_2(int n, std::initializer_list<int> entries) {
  IntMat x(n, 2);
  auto it = entries.begin();
  for (int i = 0; i < n; ++i) {
    x(i, 0) = *it++;
    x(i, 1) = *it++;
  }
  return x;
}

// Cohen's (unweighted) kappa by the standard textbook formula.
// 2 raters, 2 categories, no missingness.
double cohen_unweighted(const IntMat& x) {
  const int n = x.rows();
  // Counts of (rater1, rater2) cells.
  int n00 = 0, n01 = 0, n10 = 0, n11 = 0;
  for (int i = 0; i < n; ++i) {
    if (x(i, 0) == 0 && x(i, 1) == 0) ++n00;
    else if (x(i, 0) == 0 && x(i, 1) == 1) ++n01;
    else if (x(i, 0) == 1 && x(i, 1) == 0) ++n10;
    else if (x(i, 0) == 1 && x(i, 1) == 1) ++n11;
  }
  const double po = static_cast<double>(n00 + n11) / n;
  const double p_r1_0 = static_cast<double>(n00 + n01) / n;
  const double p_r1_1 = static_cast<double>(n10 + n11) / n;
  const double p_r2_0 = static_cast<double>(n00 + n10) / n;
  const double p_r2_1 = static_cast<double>(n01 + n11) / n;
  const double pe = p_r1_0 * p_r2_0 + p_r1_1 * p_r2_1;
  return (po - pe) / (1.0 - pe);
}

}  // namespace

TEST_CASE("estimate_available: R=2, no missing, identity weights reproduces Cohen") {
  // 10 subjects, 2 raters, 2 categories. Hand-built example.
  IntMat x = ratings_2(10, {
    0, 0,
    0, 0,
    0, 1,
    1, 0,
    1, 1,
    1, 1,
    1, 1,
    0, 0,
    1, 1,
    0, 1,
  });
  auto W = ms::loss::identity_weights(2);
  REQUIRE(W.has_value());

  auto r = ms::estimate_available(x, *W);
  REQUIRE(r.has_value());
  REQUIRE(r->estimates.size() == 3);

  const double k_textbook = cohen_unweighted(x);

  // Conger reduces to Cohen for R=2 with no missingness.
  CHECK(std::abs(r->estimates(0) - k_textbook) < 1e-9);

  // Frozen against the dev/legacy/misskappa build for the same input. Locks
  // in the agreement-matrix convention and the Fleiss / BP formulas.
  CHECK(std::abs(r->estimates(0) - 0.4) < 1e-12);                  // Conger
  CHECK(std::abs(r->estimates(1) - 13.0 / 33.0) < 1e-12);          // Fleiss = 0.3939393939...
  CHECK(std::abs(r->estimates(2) - 0.4) < 1e-12);                  // Brennan-Prediger

  REQUIRE(r->vcov.rows() == 3);
  REQUIRE(r->vcov.cols() == 3);
  CHECK(std::abs(r->vcov(0, 0) - 0.08064) < 1e-9);
  CHECK(std::abs(r->vcov(1, 1) - 0.08525581945733879) < 1e-9);
  CHECK(std::abs(r->vcov(2, 2) - 0.084) < 1e-9);
}

TEST_CASE("estimate_available: perfect agreement -> kappa = 1") {
  IntMat x = ratings_2(5, {
    0, 0,
    1, 1,
    0, 0,
    1, 1,
    0, 0,
  });
  auto W = ms::loss::identity_weights(2);
  auto r = ms::estimate_available(x, *W);
  REQUIRE(r.has_value());
  // Disagreement is 0, so kappa = 1.
  CHECK(std::abs(r->estimates(0) - 1.0) < tol);
  CHECK(std::abs(r->estimates(1) - 1.0) < tol);
  CHECK(std::abs(r->estimates(2) - 1.0) < tol);
}

TEST_CASE("estimate_available: variance is symmetric and PSD") {
  IntMat x = ratings_2(12, {
    0, 0, 0, 1, 1, 0, 1, 1,
    0, 1, 1, 0, 0, 0, 1, 1,
    1, 1, 0, 1, 1, 0, 0, 0,
  });
  auto W = ms::loss::identity_weights(2);
  auto r = ms::estimate_available(x, *W);
  REQUIRE(r.has_value());

  const RealMat& V = r->vcov;
  REQUIRE(V.rows() == 3);
  REQUIRE(V.cols() == 3);

  // Symmetric.
  CHECK((V - V.transpose()).cwiseAbs().maxCoeff() < 1e-12);

  // PSD: every eigenvalue >= 0 (up to a small numerical floor).
  Eigen::SelfAdjointEigenSolver<RealMat> es(V);
  REQUIRE(es.info() == Eigen::Success);
  CHECK(es.eigenvalues().minCoeff() > -1e-10);
}

TEST_CASE("estimate_available: missing data does not crash, returns finite") {
  IntMat x(6, 3);
  x << 0, 1, na_code,
       1, 1, 0,
       na_code, 0, 1,
       0, 0, 0,
       1, na_code, 1,
       1, 0, 0;
  auto W = ms::loss::identity_weights(2);
  auto r = ms::estimate_available(x, *W);
  REQUIRE(r.has_value());
  for (Eigen::Index i = 0; i < r->estimates.size(); ++i) {
    CHECK(std::isfinite(r->estimates(i)));
  }
}

TEST_CASE("estimate_available: too few raters -> invalid_argument") {
  IntMat x(3, 1);
  x << 0, 1, 0;
  auto W = ms::loss::identity_weights(2);
  auto r = ms::estimate_available(x, *W);
  REQUIRE(!r.has_value());
  CHECK(r.error() == ms::Error::invalid_argument);
}

TEST_CASE("estimate_available: out-of-range category code -> invalid_argument") {
  IntMat x = ratings_2(3, {0, 1, 1, 5, 0, 0});
  auto W = ms::loss::identity_weights(2);
  auto r = ms::estimate_available(x, *W);
  REQUIRE(!r.has_value());
  CHECK(r.error() == ms::Error::invalid_argument);
}

TEST_CASE("estimate_available: 3 categories with quadratic weights produces finite result") {
  IntMat x = ratings_2(8, {
    0, 0,
    1, 1,
    2, 2,
    0, 1,
    1, 2,
    0, 2,
    2, 1,
    1, 0,
  });
  RealVec v(3);
  v << 1.0, 2.0, 3.0;
  auto W = ms::loss::quadratic_weights(3, v);
  REQUIRE(W.has_value());
  auto r = ms::estimate_available(x, *W);
  REQUIRE(r.has_value());
  for (Eigen::Index i = 0; i < r->estimates.size(); ++i) {
    CHECK(std::isfinite(r->estimates(i)));
  }
}
