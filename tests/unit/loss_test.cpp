#include "doctest.h"

#include "misskappa/loss.hpp"

#include <cmath>
#include <numbers>

using misskappa::Error;
using misskappa::RealMat;
using misskappa::RealVec;
namespace loss = misskappa::loss;

namespace {

constexpr double tol = 1e-12;

bool is_symmetric(const RealMat& W) {
  return (W - W.transpose()).cwiseAbs().maxCoeff() < tol;
}

bool diagonal_is_one(const RealMat& W) {
  for (Eigen::Index i = 0; i < W.rows(); ++i) {
    if (std::abs(W(i, i) - 1.0) > tol) return false;
  }
  return true;
}

bool entries_in_unit_interval(const RealMat& W) {
  return W.minCoeff() >= -tol && W.maxCoeff() <= 1.0 + tol;
}

RealVec values_1_to_n(int n) {
  RealVec v(n);
  for (int i = 0; i < n; ++i) v(i) = static_cast<double>(i + 1);
  return v;
}

}  // namespace

TEST_CASE("identity_weights is the identity matrix") {
  for (int c : {1, 2, 3, 5, 10}) {
    auto r = loss::identity_weights(c);
    REQUIRE(r.has_value());
    CHECK((*r - RealMat::Identity(c, c)).cwiseAbs().maxCoeff() < tol);
  }
}

TEST_CASE("identity_weights rejects non-positive c") {
  CHECK(!loss::identity_weights(0).has_value());
  CHECK(!loss::identity_weights(-3).has_value());
  CHECK(loss::identity_weights(0).error() == Error::invalid_argument);
}

TEST_CASE("linear_weights known values on v = [1, 2, 3]") {
  auto r = loss::linear_weights(3, values_1_to_n(3));
  REQUIRE(r.has_value());
  const RealMat& W = *r;
  CHECK(is_symmetric(W));
  CHECK(diagonal_is_one(W));
  CHECK(entries_in_unit_interval(W));
  CHECK(std::abs(W(0, 1) - 0.5) < tol);
  CHECK(std::abs(W(0, 2) - 0.0) < tol);
  CHECK(std::abs(W(1, 2) - 0.5) < tol);
}

TEST_CASE("linear_weights wrong length -> dimension_mismatch") {
  RealVec v(2);
  v << 1.0, 2.0;
  auto r = loss::linear_weights(3, v);
  REQUIRE(!r.has_value());
  CHECK(r.error() == Error::dimension_mismatch);
}

TEST_CASE("quadratic_weights known values on v = [1, 2, 3]") {
  auto r = loss::quadratic_weights(3, values_1_to_n(3));
  REQUIRE(r.has_value());
  const RealMat& W = *r;
  CHECK(is_symmetric(W));
  CHECK(diagonal_is_one(W));
  CHECK(entries_in_unit_interval(W));
  CHECK(std::abs(W(0, 1) - 0.75) < tol);
  CHECK(std::abs(W(0, 2) - 0.0) < tol);
  CHECK(std::abs(W(1, 2) - 0.75) < tol);
}

TEST_CASE("ordinal_weights known values for c = 3") {
  auto r = loss::ordinal_weights(3);
  REQUIRE(r.has_value());
  const RealMat& W = *r;
  CHECK(is_symmetric(W));
  CHECK(diagonal_is_one(W));
  CHECK(entries_in_unit_interval(W));
  // Raw nkl*(nkl-1)/2: (0, 1, 3; 1, 0, 1; 3, 1, 0) / 3 -> 1 - that.
  CHECK(std::abs(W(0, 1) - (1.0 - 1.0 / 3.0)) < tol);
  CHECK(std::abs(W(0, 2) - 0.0) < tol);
  CHECK(std::abs(W(1, 2) - (1.0 - 1.0 / 3.0)) < tol);
}

TEST_CASE("radical_weights known values on v = [1, 2, 3]") {
  auto r = loss::radical_weights(3, values_1_to_n(3));
  REQUIRE(r.has_value());
  const RealMat& W = *r;
  CHECK(is_symmetric(W));
  CHECK(diagonal_is_one(W));
  CHECK(entries_in_unit_interval(W));
  const double sqrt2 = std::sqrt(2.0);
  CHECK(std::abs(W(0, 1) - (1.0 - 1.0 / sqrt2)) < tol);
  CHECK(std::abs(W(0, 2) - 0.0) < tol);
  CHECK(std::abs(W(1, 2) - (1.0 - 1.0 / sqrt2)) < tol);
}

TEST_CASE("ratio_weights symmetric / diagonal / range on v = [1, 2, 3]") {
  auto r = loss::ratio_weights(3, values_1_to_n(3));
  REQUIRE(r.has_value());
  CHECK(is_symmetric(*r));
  CHECK(diagonal_is_one(*r));
  CHECK(entries_in_unit_interval(*r));
}

TEST_CASE("circular_weights symmetric / diagonal / range on v = [1, 2, 3]") {
  auto r = loss::circular_weights(3, values_1_to_n(3));
  REQUIRE(r.has_value());
  CHECK(is_symmetric(*r));
  CHECK(diagonal_is_one(*r));
  CHECK(entries_in_unit_interval(*r));
}

TEST_CASE("bipolar_weights symmetric / diagonal / range on v = [1, 2, 3]") {
  auto r = loss::bipolar_weights(3, values_1_to_n(3));
  REQUIRE(r.has_value());
  CHECK(is_symmetric(*r));
  CHECK(diagonal_is_one(*r));
  CHECK(entries_in_unit_interval(*r));
}

TEST_CASE("Degenerate value vector collapses to identity") {
  RealVec v(3);
  v << 2.0, 2.0, 2.0;
  for (auto fn : {
           +[](int c, const RealVec& v) { return loss::linear_weights(c, v); },
           +[](int c, const RealVec& v) { return loss::quadratic_weights(c, v); },
           +[](int c, const RealVec& v) { return loss::radical_weights(c, v); },
       }) {
    auto r = fn(3, v);
    REQUIRE(r.has_value());
    CHECK((*r - RealMat::Identity(3, 3)).cwiseAbs().maxCoeff() < tol);
  }
}

// --- Continuous loss kernels ---

TEST_CASE("identity_loss: 0 on equality, 1 on disagreement") {
  auto r = loss::identity_loss();
  REQUIRE(r.has_value());
  CHECK(r->compute(1.0, 1.0, r->min_val, r->max_val) == 0.0);
  CHECK(r->compute(1.0, 2.0, r->min_val, r->max_val) == 1.0);
}

TEST_CASE("linear_loss: known values for [0, 2]") {
  auto r = loss::linear_loss(0.0, 2.0);
  REQUIRE(r.has_value());
  CHECK(std::abs(r->compute(0.0, 0.0, r->min_val, r->max_val) - 0.0) < tol);
  CHECK(std::abs(r->compute(0.0, 1.0, r->min_val, r->max_val) - 0.5) < tol);
  CHECK(std::abs(r->compute(0.0, 2.0, r->min_val, r->max_val) - 1.0) < tol);
}

TEST_CASE("quadratic_loss: known values for [0, 2]") {
  auto r = loss::quadratic_loss(0.0, 2.0);
  REQUIRE(r.has_value());
  CHECK(std::abs(r->compute(0.0, 1.0, r->min_val, r->max_val) - 0.25) < tol);
  CHECK(std::abs(r->compute(0.0, 2.0, r->min_val, r->max_val) - 1.0) < tol);
}

TEST_CASE("Continuous: zero-range collapses to identity behaviour") {
  auto r = loss::linear_loss(1.0, 1.0);
  REQUIRE(r.has_value());
  CHECK(r->compute(1.0, 1.0, r->min_val, r->max_val) == 0.0);
  CHECK(r->compute(1.0, 2.0, r->min_val, r->max_val) == 1.0);
}
