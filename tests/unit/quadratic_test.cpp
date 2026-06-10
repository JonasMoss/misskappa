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
  REQUIRE(r->psi.rows() == x.rows());
  REQUIRE(r->psi.cols() == 3);
  const RealMat psi_vcov =
      (r->psi.transpose() * r->psi) / std::pow(static_cast<double>(x.rows()), 2);
  CHECK((psi_vcov - V).cwiseAbs().maxCoeff() < 1e-12);
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

TEST_CASE("estimate_quadratic raw: no co-observed pair -> not_identified") {
  RealMat x(2, 3);
  x << 1, na_d, na_d,
       na_d, 2, na_d;
  RealVec v(5);
  v << 1.0, 2.0, 3.0, 4.0, 5.0;
  auto r = ms::estimate_quadratic(x, v);
  REQUIRE(!r.has_value());
  CHECK(r.error() == ms::Error::not_identified);
}
