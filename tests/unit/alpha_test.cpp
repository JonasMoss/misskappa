#include "doctest.h"

#include "misskappa/estimate.hpp"

#include <cmath>

using misskappa::EmOptions;
using misskappa::IntMat;
using misskappa::na_code;
using misskappa::RealMat;
using misskappa::RealVec;
namespace ms = misskappa;

namespace {

IntMat complete_items() {
  IntMat x(6, 3);
  x << 0, 0, 1,
       1, 1, 1,
       2, 1, 2,
       1, 2, 2,
       0, 1, 0,
       2, 2, 1;
  return x;
}

IntMat missing_items() {
  IntMat x(7, 3);
  x << 0, 0, 1,
       1, na_code, 1,
       2, 1, 2,
       na_code, 2, 2,
       0, 1, 0,
       2, 2, na_code,
       1, 1, 1;
  return x;
}

RealVec three_scores() {
  RealVec v(3);
  v << 1.0, 2.0, 3.0;
  return v;
}

double complete_alpha_by_covariance(const IntMat& x, const RealVec& values) {
  const int n = static_cast<int>(x.rows());
  const int R = static_cast<int>(x.cols());
  RealMat y(n, R);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < R; ++j) y(i, j) = values(x(i, j));
  }
  RealVec mu = y.colwise().mean().transpose();
  RealMat sigma = RealMat::Zero(R, R);
  for (int i = 0; i < n; ++i) {
    const RealVec centered = y.row(i).transpose() - mu;
    sigma.noalias() += centered * centered.transpose();
  }
  sigma /= static_cast<double>(n);
  const double t1 = sigma.sum();
  const double t2 = sigma.diagonal().sum();
  return (static_cast<double>(R) / static_cast<double>(R - 1)) * (1.0 - t2 / t1);
}

void check_influence_reconstructs_vcov(const ms::Estimation& e, int n) {
  REQUIRE(e.psi.rows() == n);
  REQUIRE(e.psi.cols() == e.estimates.size());
  const RealMat psi_vcov =
      (e.psi.transpose() * e.psi) / std::pow(static_cast<double>(n), 2);
  CHECK((psi_vcov - e.vcov).cwiseAbs().maxCoeff() < 1e-10);
}

}  // namespace

TEST_CASE("estimate_alpha_available: complete data matches covariance alpha") {
  const IntMat x = complete_items();
  const RealVec values = three_scores();
  auto r = ms::estimate_alpha_available(x, values);
  REQUIRE(r.has_value());

  CHECK(std::abs(r->estimates(0) - complete_alpha_by_covariance(x, values)) < 1e-12);
  CHECK(r->vcov.rows() == 1);
  CHECK(r->vcov.cols() == 1);
  CHECK(r->vcov(0, 0) >= -1e-12);
  check_influence_reconstructs_vcov(*r, static_cast<int>(x.rows()));
}

TEST_CASE("estimate_alpha_fiml: complete data matches available alpha") {
  const IntMat x = complete_items();
  const RealVec values = three_scores();
  auto av = ms::estimate_alpha_available(x, values);
  REQUIRE(av.has_value());
  auto ml = ms::estimate_alpha_fiml(x, values, EmOptions{});
  REQUIRE(ml.has_value());

  CHECK(std::abs(ml->estimates(0) - av->estimates(0)) < 1e-10);
  CHECK(ml->vcov(0, 0) >= -1e-12);
  check_influence_reconstructs_vcov(*ml, static_cast<int>(x.rows()));
}

TEST_CASE("estimate_alpha_fiml: missing categorical fixture is finite") {
  const IntMat x = missing_items();
  const RealVec values = three_scores();
  auto av = ms::estimate_alpha_available(x, values);
  REQUIRE(av.has_value());
  auto ml = ms::estimate_alpha_fiml(x, values, EmOptions{});
  REQUIRE(ml.has_value());

  CHECK(std::isfinite(av->estimates(0)));
  CHECK(std::isfinite(ml->estimates(0)));
  CHECK(av->vcov(0, 0) >= -1e-12);
  CHECK(ml->vcov(0, 0) >= -1e-12);
  check_influence_reconstructs_vcov(*av, static_cast<int>(x.rows()));
  check_influence_reconstructs_vcov(*ml, static_cast<int>(x.rows()));
}

TEST_CASE("estimate_alpha: invalid inputs are rejected") {
  const RealVec values = three_scores();

  IntMat too_few_items(3, 1);
  too_few_items << 0, 1, 2;
  auto few = ms::estimate_alpha_available(too_few_items, values);
  REQUIRE(!few.has_value());
  CHECK(few.error() == ms::Error::invalid_argument);

  IntMat bad = complete_items();
  bad(0, 0) = 5;
  auto out_of_range = ms::estimate_alpha_fiml(bad, values, EmOptions{});
  REQUIRE(!out_of_range.has_value());
  CHECK(out_of_range.error() == ms::Error::invalid_argument);
}
