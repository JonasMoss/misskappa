#include "doctest.h"

#include "misskappa/estimate.hpp"
#include "misskappa/loss.hpp"

#include <Eigen/Eigenvalues>
#include <cmath>
#include <limits>

using misskappa::RealMat;
using misskappa::RealVec;
namespace ms = misskappa;

namespace {

constexpr double tol = 1e-9;
const double na_d = std::numeric_limits<double>::quiet_NaN();

RealMat ten_subject_3rater_continuous() {
  // Matches the R fixture used to capture reference values; range 0.9 .. 4.6.
  RealMat x(10, 3);
  x <<
    1.0, 1.2, 0.9,
    3.5, 3.4, 3.6,
    2.0, 2.1, 1.9,
    4.0, 4.2, 3.8,
    1.5, 1.4, 1.6,
    3.0, 2.9, 3.1,
    2.5, 2.6, 2.4,
    4.5, 4.4, 4.6,
    1.0, na_d, 1.1,
    3.5, 3.6, na_d;
  return x;
}

// Legacy convention: min/max are taken over the finite entries of x.
std::pair<double, double> finite_range(const RealMat& x) {
  double mn = std::numeric_limits<double>::infinity();
  double mx = -std::numeric_limits<double>::infinity();
  for (Eigen::Index i = 0; i < x.rows(); ++i) {
    for (Eigen::Index j = 0; j < x.cols(); ++j) {
      if (!std::isfinite(x(i, j))) continue;
      if (x(i, j) < mn) mn = x(i, j);
      if (x(i, j) > mx) mx = x(i, j);
    }
  }
  return {mn, mx};
}

}  // namespace

TEST_CASE("Continuous available-case matches legacy: quadratic loss, fixture") {
  RealMat x = ten_subject_3rater_continuous();
  auto [mn, mx] = finite_range(x);
  auto loss = ms::loss::quadratic_loss(mn, mx);
  REQUIRE(loss.has_value());

  auto r = ms::estimate_available_continuous(x, *loss);
  REQUIRE(r.has_value());
  REQUIRE(r->estimates.size() == 2);

  // Frozen against dev/legacy/misskappa kappa_continuous(method="available",
  // weight="quadratic").
  CHECK(std::abs(r->estimates(0) - 0.9893902893012867) < tol);  // Conger
  CHECK(std::abs(r->estimates(1) - 0.9893337228643995) < tol);  // Fleiss

  CHECK(std::abs(r->vcov(0, 0) - 9.630587853753491e-06) < 1e-12);
  CHECK(std::abs(r->vcov(1, 1) - 9.780939091210242e-06) < 1e-12);
  CHECK(std::abs(r->vcov(0, 1) - 9.702861043956929e-06) < 1e-12);
}

TEST_CASE("Continuous IPW matches legacy: quadratic loss") {
  RealMat x = ten_subject_3rater_continuous();
  auto [mn, mx] = finite_range(x);
  auto loss = ms::loss::quadratic_loss(mn, mx);

  auto r = ms::estimate_ipw_continuous(x, *loss);
  REQUIRE(r.has_value());
  CHECK(std::abs(r->estimates(0) - 0.9889592202554336) < tol);
  CHECK(std::abs(r->estimates(1) - 0.9888908233278307) < tol);
}

TEST_CASE("Continuous Gwet matches legacy: quadratic loss") {
  RealMat x = ten_subject_3rater_continuous();
  auto [mn, mx] = finite_range(x);
  auto loss = ms::loss::quadratic_loss(mn, mx);

  auto r = ms::estimate_gwet_continuous(x, *loss);
  REQUIRE(r.has_value());
  CHECK(std::abs(r->estimates(0) - 0.9893896826901533) < tol);
  CHECK(std::abs(r->estimates(1) - 0.9893239524499280) < tol);
}

TEST_CASE("Continuous available-case matches legacy: linear loss") {
  RealMat x = ten_subject_3rater_continuous();
  auto [mn, mx] = finite_range(x);
  auto loss = ms::loss::linear_loss(mn, mx);

  auto r = ms::estimate_available_continuous(x, *loss);
  REQUIRE(r.has_value());
  CHECK(std::abs(r->estimates(0) - 0.8865710560625815) < tol);
  CHECK(std::abs(r->estimates(1) - 0.8857159516625119) < tol);
}

TEST_CASE("Continuous available-case: identity loss handles negative values") {
  // Identity (binary 0/1) loss on continuous data: agreement only when
  // values are exactly equal. With the fixture's near-equal but distinct
  // values it shows the disagreement-dominated regime; matches legacy.
  RealMat x = ten_subject_3rater_continuous();
  auto loss = ms::loss::identity_loss();

  auto r = ms::estimate_available_continuous(x, *loss);
  REQUIRE(r.has_value());
  CHECK(std::abs(r->estimates(0) - (-0.0038461538461538)) < tol);
  CHECK(std::abs(r->estimates(1) - (-0.0453333333333334)) < tol);
}

TEST_CASE("Continuous: complete data, all three methods agree exactly") {
  // No NaN entries: pi_j = 1 for all raters, so available / IPW / Gwet
  // collapse to the same estimator (matches the categorical-raw property).
  RealMat x(8, 3);
  x <<
    1.0, 1.2, 0.9,
    3.5, 3.4, 3.6,
    2.0, 2.1, 1.9,
    4.0, 4.2, 3.8,
    1.5, 1.4, 1.6,
    3.0, 2.9, 3.1,
    2.5, 2.6, 2.4,
    4.5, 4.4, 4.6;
  auto [mn, mx] = finite_range(x);
  auto loss = ms::loss::quadratic_loss(mn, mx);

  auto r_av  = ms::estimate_available_continuous(x, *loss);
  auto r_ipw = ms::estimate_ipw_continuous(x, *loss);
  auto r_gw  = ms::estimate_gwet_continuous(x, *loss);
  REQUIRE(r_av.has_value());
  REQUIRE(r_ipw.has_value());
  REQUIRE(r_gw.has_value());
  CHECK((r_av->estimates - r_ipw->estimates).cwiseAbs().maxCoeff() < tol);
  CHECK((r_av->estimates - r_gw->estimates).cwiseAbs().maxCoeff() < tol);
  CHECK((r_av->vcov - r_ipw->vcov).cwiseAbs().maxCoeff() < tol);
}

TEST_CASE("Continuous: perfect agreement -> kappa = 1") {
  RealMat x(5, 2);
  x << 1.5, 1.5,
       2.7, 2.7,
       3.1, 3.1,
       0.5, 0.5,
       4.0, 4.0;
  auto loss = ms::loss::quadratic_loss(0.5, 4.0);

  auto r = ms::estimate_available_continuous(x, *loss);
  REQUIRE(r.has_value());
  CHECK(std::abs(r->estimates(0) - 1.0) < tol);
  CHECK(std::abs(r->estimates(1) - 1.0) < tol);
}

TEST_CASE("Continuous: variance is symmetric and PSD") {
  RealMat x = ten_subject_3rater_continuous();
  auto [mn, mx] = finite_range(x);
  auto loss = ms::loss::quadratic_loss(mn, mx);

  auto r = ms::estimate_available_continuous(x, *loss);
  REQUIRE(r.has_value());
  const RealMat& V = r->vcov;
  CHECK((V - V.transpose()).cwiseAbs().maxCoeff() < 1e-12);
  Eigen::SelfAdjointEigenSolver<RealMat> es(V);
  REQUIRE(es.info() == Eigen::Success);
  CHECK(es.eigenvalues().minCoeff() > -1e-10);

  REQUIRE(r->psi.rows() == x.rows());
  REQUIRE(r->psi.cols() == 2);
  const RealMat psi_vcov =
      (r->psi.transpose() * r->psi) / std::pow(static_cast<double>(x.rows()), 2);
  CHECK((psi_vcov - V).cwiseAbs().maxCoeff() < 1e-12);
}

TEST_CASE("Continuous: too few raters -> invalid_argument") {
  RealMat x(3, 1);
  x << 1.0, 2.0, 3.0;
  auto loss = ms::loss::quadratic_loss(1.0, 3.0);
  auto r = ms::estimate_available_continuous(x, *loss);
  REQUIRE(!r.has_value());
  CHECK(r.error() == ms::Error::invalid_argument);
}

TEST_CASE("Continuous IPW rejects a rater with zero observations") {
  RealMat x(4, 3);
  x << 1.0, 2.0, na_d,
       2.0, 1.5, na_d,
       3.0, 2.5, na_d,
       1.5, 2.0, na_d;
  auto loss = ms::loss::quadratic_loss(1.0, 3.0);
  auto r = ms::estimate_ipw_continuous(x, *loss);
  REQUIRE(!r.has_value());
  CHECK(r.error() == ms::Error::singular_weight);
}
