#include "doctest.h"

#include "misskappa/loss.hpp"

#include <cmath>
#include <limits>

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

// --- Validation guards on the metric weight factories ---

TEST_CASE("Metric weight factories reject non-positive c") {
  const RealVec v = values_1_to_n(3);
  CHECK(loss::linear_weights(0, v).error() == Error::invalid_argument);
  CHECK(loss::quadratic_weights(-1, v).error() == Error::invalid_argument);
  CHECK(loss::ordinal_weights(0).error() == Error::invalid_argument);
  CHECK(loss::radical_weights(0, v).error() == Error::invalid_argument);
  CHECK(loss::ratio_weights(0, v).error() == Error::invalid_argument);
  CHECK(loss::circular_weights(0, v).error() == Error::invalid_argument);
  CHECK(loss::bipolar_weights(-2, v).error() == Error::invalid_argument);
}

TEST_CASE("Metric weight factories reject a mismatched value-vector length") {
  RealVec v2(2);
  v2 << 1.0, 2.0;
  CHECK(loss::quadratic_weights(3, v2).error() == Error::dimension_mismatch);
  CHECK(loss::radical_weights(3, v2).error() == Error::dimension_mismatch);
  CHECK(loss::ratio_weights(3, v2).error() == Error::dimension_mismatch);
  CHECK(loss::circular_weights(3, v2).error() == Error::dimension_mismatch);
  CHECK(loss::bipolar_weights(3, v2).error() == Error::dimension_mismatch);
}

// --- Degenerate value vectors collapse to the identity matrix ---

TEST_CASE("ratio / circular / bipolar collapse to identity on flat scores") {
  const RealMat I3 = RealMat::Identity(3, 3);
  RealVec flat(3);
  flat << 2.0, 2.0, 2.0;
  CHECK((*loss::ratio_weights(3, flat) - I3).cwiseAbs().maxCoeff() < tol);
  CHECK((*loss::circular_weights(3, flat) - I3).cwiseAbs().maxCoeff() < tol);
  CHECK((*loss::bipolar_weights(3, flat) - I3).cwiseAbs().maxCoeff() < tol);
}

TEST_CASE("ratio_weights collapses to identity when scores are antisymmetric") {
  const RealMat I3 = RealMat::Identity(3, 3);
  RealVec antisym(3);
  antisym << -1.0, 0.0, 1.0;  // v_max + v_min == 0
  CHECK((*loss::ratio_weights(3, antisym) - I3).cwiseAbs().maxCoeff() < tol);
}

TEST_CASE("ordinal_weights collapses to identity for a single category") {
  auto r = loss::ordinal_weights(1);
  REQUIRE(r.has_value());
  CHECK(std::abs((*r)(0, 0) - 1.0) < tol);
}

// --- Special-case branches inside the ratio / bipolar kernels ---

TEST_CASE("ratio_weights special-cases a zero-sum category pair") {
  RealVec v(3);
  v << -1.0, 1.0, 2.0;  // pair (0, 1) sums to zero
  auto r = loss::ratio_weights(3, v);
  REQUIRE(r.has_value());
  CHECK(std::abs((*r)(0, 1) - 1.0) < tol);
}

TEST_CASE("bipolar_weights handles a zero-denominator pair") {
  RealVec v(3);
  v << 1.0, 1.0, 3.0;  // categories 0 and 1 tie at the minimum -> den == 0
  auto r = loss::bipolar_weights(3, v);
  REQUIRE(r.has_value());
  CHECK(is_symmetric(*r));
  CHECK(std::abs((*r)(0, 1) - 1.0) < tol);  // skipped pair keeps raw weight 0
}

// --- Continuous radical / ratio loss kernels ---

TEST_CASE("radical_loss: known values for [0, 4]") {
  auto r = loss::radical_loss(0.0, 4.0);
  REQUIRE(r.has_value());
  CHECK(std::abs(r->compute(0.0, 0.0, r->min_val, r->max_val) - 0.0) < tol);
  CHECK(std::abs(r->compute(0.0, 1.0, r->min_val, r->max_val) - 0.5) < tol);
  CHECK(std::abs(r->compute(0.0, 4.0, r->min_val, r->max_val) - 1.0) < tol);
}

TEST_CASE("radical_loss: zero range collapses to identity") {
  auto r = loss::radical_loss(2.0, 2.0);
  REQUIRE(r.has_value());
  CHECK(r->compute(2.0, 2.0, r->min_val, r->max_val) == 0.0);
  CHECK(r->compute(2.0, 5.0, r->min_val, r->max_val) == 1.0);
}

TEST_CASE("ratio_loss: known values for [1, 3]") {
  auto r = loss::ratio_loss(1.0, 3.0);
  REQUIRE(r.has_value());
  CHECK(std::abs(r->compute(2.0, 2.0, r->min_val, r->max_val) - 0.0) < tol);
  CHECK(std::abs(r->compute(1.0, 3.0, r->min_val, r->max_val) - 1.0) < tol);
  CHECK(r->compute(0.0, 0.0, r->min_val, r->max_val) == 0.0);  // a + b == 0
}

TEST_CASE("ratio_loss: degenerate parameters collapse to identity") {
  auto r1 = loss::ratio_loss(-1.0, 1.0);  // v_max + v_min == 0
  REQUIRE(r1.has_value());
  CHECK(r1->compute(1.0, 1.0, r1->min_val, r1->max_val) == 0.0);
  CHECK(r1->compute(1.0, 2.0, r1->min_val, r1->max_val) == 1.0);

  auto r2 = loss::ratio_loss(2.0, 2.0);  // zero range -> den_sq == 0
  REQUIRE(r2.has_value());
  CHECK(r2->compute(2.0, 2.0, r2->min_val, r2->max_val) == 0.0);
  CHECK(r2->compute(2.0, 9.0, r2->min_val, r2->max_val) == 1.0);
}

TEST_CASE("quadratic_loss: zero range collapses to identity") {
  auto r = loss::quadratic_loss(1.0, 1.0);
  REQUIRE(r.has_value());
  CHECK(r->compute(1.0, 1.0, r->min_val, r->max_val) == 0.0);
  CHECK(r->compute(1.0, 2.0, r->min_val, r->max_val) == 1.0);
}

// --- Component-separable vector loss guards ---

TEST_CASE("make_vector_loss rejects degenerate feature weights") {
  RealVec empty(0);
  CHECK(loss::hamming_vector_loss(empty).error() == Error::invalid_argument);

  RealVec negative(2);
  negative << 1.0, -1.0;
  CHECK(loss::absolute_vector_loss(negative).error() == Error::invalid_argument);

  RealVec nan_weight(2);
  nan_weight << 1.0, std::numeric_limits<double>::quiet_NaN();
  CHECK(loss::squared_vector_loss(nan_weight).error() == Error::invalid_argument);

  RealVec all_zero(2);
  all_zero << 0.0, 0.0;
  CHECK(loss::rms_vector_loss(all_zero).error() == Error::invalid_argument);
}

TEST_CASE("rms_vector_loss transform handles non-positive arguments") {
  RealVec fw(2);
  fw << 1.0, 1.0;
  auto r = loss::rms_vector_loss(fw);
  REQUIRE(r.has_value());
  CHECK(r->transform(0.0) == 0.0);
  CHECK(r->transform(-3.0) == 0.0);
  CHECK(std::abs(r->transform(4.0) - 2.0) < tol);
  CHECK(r->transform_derivative(0.0) == 0.0);
  CHECK(std::abs(r->transform_derivative(4.0) - 0.25) < tol);
}
