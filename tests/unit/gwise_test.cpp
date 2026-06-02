#include "doctest.h"

#include "misskappa/estimate.hpp"
#include "misskappa/loss.hpp"

#include <Eigen/Eigenvalues>
#include <cmath>

using misskappa::IntMat;
using misskappa::RealMat;
namespace ms = misskappa;

namespace {

constexpr double tol = 1e-9;

IntMat frechet_fixture_int() {
  IntMat x(4, 5);
  x << 0, 0, 1, 0, 0,
       0, 1, 2, 1, 1,
       1, 0, 0, 0, 0,
       1, 2, 3, 3, 4;
  return x;
}

RealMat frechet_fixture_real() {
  return frechet_fixture_int().cast<double>() + RealMat::Constant(4, 5, 1.0);
}

}  // namespace

TEST_CASE("G-wise nominal Frechet: g=2 matches pairwise identity kappa") {
  IntMat x = frechet_fixture_int();
  auto distance = ms::loss::frechet_nominal_distance(5);
  auto weights = ms::loss::identity_weights(5);
  REQUIRE(distance.has_value());
  REQUIRE(weights.has_value());

  auto gwise = ms::estimate_gwise(x, *distance, ms::GwiseOptions{2});
  auto pairwise = ms::estimate_available(x, *weights);
  REQUIRE(gwise.has_value());
  REQUIRE(pairwise.has_value());

  CHECK(std::abs(gwise->estimates(0) - pairwise->estimates(0)) < tol);
  CHECK(std::abs(gwise->estimates(1) - pairwise->estimates(1)) < tol);
  CHECK((gwise->vcov - pairwise->vcov.block(0, 0, 2, 2)).cwiseAbs().maxCoeff() < tol);
}

TEST_CASE("G-wise nominal Frechet: all-rater mode disagreement fixture") {
  IntMat x = frechet_fixture_int();
  auto distance = ms::loss::frechet_nominal_distance(5);
  REQUIRE(distance.has_value());

  auto r = ms::estimate_gwise(x, *distance);
  REQUIRE(r.has_value());
  REQUIRE(r->estimates.size() == 2);
  CHECK(std::abs(r->estimates(0) - 0.23353293413173815) < tol);
  CHECK(std::abs(r->estimates(1) - 0.21735788407113500) < tol);
  CHECK(r->psi.rows() == x.rows());
  CHECK(r->psi.cols() == 2);
}

TEST_CASE("G-wise absolute Frechet: all-rater median disagreement fixture") {
  RealMat x = frechet_fixture_real();
  auto distance = ms::loss::frechet_absolute_distance();
  REQUIRE(distance.has_value());

  auto r = ms::estimate_gwise_continuous(x, *distance);
  REQUIRE(r.has_value());
  CHECK(std::abs(r->estimates(0) - 0.45877378435517835) < tol);
  CHECK(std::abs(r->estimates(1) - 0.44821878125323300) < tol);

  const RealMat psi_vcov = (r->psi.transpose() * r->psi) / std::pow(static_cast<double>(x.rows()), 2);
  CHECK((psi_vcov - r->vcov).cwiseAbs().maxCoeff() < 1e-12);
}

TEST_CASE("G-wise Hubert disagreement: perfect agreement returns one") {
  IntMat x(3, 3);
  x << 0, 0, 0,
       1, 1, 1,
       2, 2, 2;
  auto distance = ms::loss::hubert_categorical_distance(3);
  REQUIRE(distance.has_value());

  auto r = ms::estimate_gwise(x, *distance);
  REQUIRE(r.has_value());
  CHECK(std::abs(r->estimates(0) - 1.0) < tol);
  CHECK(std::abs(r->estimates(1) - 1.0) < tol);
}

TEST_CASE("G-wise estimator rejects missing and oversized direct jobs") {
  IntMat x = frechet_fixture_int();
  auto distance = ms::loss::frechet_nominal_distance(5);
  REQUIRE(distance.has_value());

  x(0, 0) = ms::na_code;
  auto missing = ms::estimate_gwise(x, *distance);
  REQUIRE(!missing.has_value());
  CHECK(missing.error() == ms::Error::invalid_argument);

  x = frechet_fixture_int();
  auto oversized = ms::estimate_gwise(x, *distance, ms::GwiseOptions{5, 100});
  REQUIRE(!oversized.has_value());
  CHECK(oversized.error() == ms::Error::not_supported);
}

TEST_CASE("G-wise covariance is symmetric and PSD") {
  RealMat x = frechet_fixture_real();
  auto distance = ms::loss::frechet_quadratic_distance();
  REQUIRE(distance.has_value());

  auto r = ms::estimate_gwise_continuous(x, *distance);
  REQUIRE(r.has_value());
  CHECK((r->vcov - r->vcov.transpose()).cwiseAbs().maxCoeff() < 1e-12);

  Eigen::SelfAdjointEigenSolver<RealMat> es(r->vcov);
  REQUIRE(es.info() == Eigen::Success);
  CHECK(es.eigenvalues().minCoeff() > -1e-10);
}
