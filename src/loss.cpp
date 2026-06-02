#include "misskappa/loss.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <vector>

namespace misskappa::loss {

namespace {

constexpr double zero_tol = 1e-9;

// Convenience: build a c x c agreement matrix from a pairwise functor.
template <typename F>
RealMat fill_pairwise(int c, F&& fn) {
  RealMat W(c, c);
  for (int i = 0; i < c; ++i) {
    for (int j = 0; j < c; ++j) {
      W(i, j) = fn(i, j);
    }
  }
  return W;
}

bool length_ok(const RealVec& v, int c) {
  return v.size() == static_cast<Eigen::Index>(c);
}

}  // namespace

// All categorical factories return AGREEMENT matrices (1 on the diagonal,
// values in [0, 1] off-diagonal). The estimator turns these into the kappa
// numerator / denominator. Continuous factories below return LOSS kernels
// (0 on equality, larger for greater disagreement) — the legacy convention,
// preserved deliberately so estimator code carrying both is unsurprised.

Result<RealMat> identity_weights(int c) {
  if (c <= 0) return std::unexpected(Error::invalid_argument);
  return RealMat(RealMat::Identity(c, c));
}

Result<RealMat> linear_weights(int c, const RealVec& v) {
  if (c <= 0) return std::unexpected(Error::invalid_argument);
  if (!length_ok(v, c)) return std::unexpected(Error::dimension_mismatch);
  const double range = std::abs(v.maxCoeff() - v.minCoeff());
  if (range < zero_tol) return identity_weights(c);
  return fill_pairwise(c, [&](int i, int j) {
    return 1.0 - std::abs(v(i) - v(j)) / range;
  });
}

Result<RealMat> quadratic_weights(int c, const RealVec& v) {
  if (c <= 0) return std::unexpected(Error::invalid_argument);
  if (!length_ok(v, c)) return std::unexpected(Error::dimension_mismatch);
  const double range_sq = std::pow(v.maxCoeff() - v.minCoeff(), 2);
  if (range_sq < zero_tol) return identity_weights(c);
  return fill_pairwise(c, [&](int i, int j) {
    return 1.0 - std::pow(v(i) - v(j), 2) / range_sq;
  });
}

Result<RealMat> ordinal_weights(int c) {
  if (c <= 0) return std::unexpected(Error::invalid_argument);
  RealMat W = fill_pairwise(c, [](int i, int j) {
    const double nkl = std::abs(i - j) + 1.0;
    return nkl * (nkl - 1.0) / 2.0;
  });
  const double max_w = W.maxCoeff();
  if (max_w < zero_tol) return identity_weights(c);
  return RealMat(RealMat::Constant(c, c, 1.0) - W / max_w);
}

Result<RealMat> radical_weights(int c, const RealVec& v) {
  if (c <= 0) return std::unexpected(Error::invalid_argument);
  if (!length_ok(v, c)) return std::unexpected(Error::dimension_mismatch);
  const double range_sqrt = std::sqrt(std::abs(v.maxCoeff() - v.minCoeff()));
  if (range_sqrt < zero_tol) return identity_weights(c);
  return fill_pairwise(c, [&](int i, int j) {
    return 1.0 - std::sqrt(std::abs(v(i) - v(j))) / range_sqrt;
  });
}

Result<RealMat> ratio_weights(int c, const RealVec& v) {
  if (c <= 0) return std::unexpected(Error::invalid_argument);
  if (!length_ok(v, c)) return std::unexpected(Error::dimension_mismatch);
  const double v_min = v.minCoeff();
  const double v_max = v.maxCoeff();
  if (std::abs(v_max + v_min) < zero_tol) return identity_weights(c);
  const double den_term = (v_max - v_min) / (v_max + v_min);
  const double den_sq = std::pow(den_term, 2);
  if (den_sq < zero_tol) return identity_weights(c);
  return fill_pairwise(c, [&](int i, int j) {
    if (std::abs(v(i) + v(j)) < zero_tol) return 1.0;
    const double num_term = (v(i) - v(j)) / (v(i) + v(j));
    return 1.0 - std::pow(num_term, 2) / den_sq;
  });
}

Result<RealMat> circular_weights(int c, const RealVec& v) {
  if (c <= 0) return std::unexpected(Error::invalid_argument);
  if (!length_ok(v, c)) return std::unexpected(Error::dimension_mismatch);
  const double U = v.maxCoeff() - v.minCoeff() + 1.0;
  if (U < zero_tol) return identity_weights(c);
  RealMat W = fill_pairwise(c, [&](int i, int j) {
    return std::pow(std::sin(std::numbers::pi * (v(i) - v(j)) / U), 2);
  });
  const double max_w = W.maxCoeff();
  if (max_w < zero_tol) return identity_weights(c);
  return RealMat(RealMat::Constant(c, c, 1.0) - W / max_w);
}

Result<RealMat> bipolar_weights(int c, const RealVec& v) {
  if (c <= 0) return std::unexpected(Error::invalid_argument);
  if (!length_ok(v, c)) return std::unexpected(Error::dimension_mismatch);
  const double v_min = v.minCoeff();
  const double v_max = v.maxCoeff();
  RealMat W = RealMat::Zero(c, c);
  for (int i = 0; i < c; ++i) {
    for (int j = 0; j < c; ++j) {
      if (i == j) continue;
      const double den = ((v(i) + v(j)) - 2.0 * v_min) * (2.0 * v_max - (v(i) + v(j)));
      if (std::abs(den) > zero_tol) {
        W(i, j) = std::pow(v(i) - v(j), 2) / den;
      }
    }
  }
  const double max_w = W.maxCoeff();
  if (max_w < zero_tol) return identity_weights(c);
  return RealMat(RealMat::Constant(c, c, 1.0) - W / max_w);
}

// --- Continuous loss kernels ---
//
// These take two real-valued ratings (a, b) and return a loss in [0, 1]:
// 0 when a equals b, larger for greater disagreement. The legacy
// convention is to scale by the range so the maximum loss is 1.

namespace continuous_kernels {

double identity(double a, double b, double /*min_val*/, double /*max_val*/) {
  return (std::abs(a - b) < zero_tol) ? 0.0 : 1.0;
}

double linear(double a, double b, double min_val, double max_val) {
  const double range = std::abs(max_val - min_val);
  return std::abs(a - b) / range;
}

double quadratic(double a, double b, double min_val, double max_val) {
  const double range_sq = std::pow(max_val - min_val, 2);
  return std::pow(a - b, 2) / range_sq;
}

double radical(double a, double b, double min_val, double max_val) {
  const double range_sqrt = std::sqrt(std::abs(max_val - min_val));
  return std::sqrt(std::abs(a - b)) / range_sqrt;
}

double ratio(double a, double b, double min_val, double max_val) {
  if (std::abs(a + b) < zero_tol) return 0.0;
  const double den_term = (max_val - min_val) / (max_val + min_val);
  const double den_sq = std::pow(den_term, 2);
  const double num_term = (a - b) / (a + b);
  return std::pow(num_term, 2) / den_sq;
}

}  // namespace continuous_kernels

Result<ContinuousLoss> identity_loss() {
  return ContinuousLoss{0.0, 0.0, &continuous_kernels::identity};
}

Result<ContinuousLoss> linear_loss(double min_val, double max_val) {
  if (std::abs(max_val - min_val) < zero_tol) return identity_loss();
  return ContinuousLoss{min_val, max_val, &continuous_kernels::linear};
}

Result<ContinuousLoss> quadratic_loss(double min_val, double max_val) {
  if (std::abs(max_val - min_val) < zero_tol) return identity_loss();
  return ContinuousLoss{min_val, max_val, &continuous_kernels::quadratic};
}

Result<ContinuousLoss> radical_loss(double min_val, double max_val) {
  if (std::abs(max_val - min_val) < zero_tol) return identity_loss();
  return ContinuousLoss{min_val, max_val, &continuous_kernels::radical};
}

Result<ContinuousLoss> ratio_loss(double min_val, double max_val) {
  if (std::abs(max_val + min_val) < zero_tol) return identity_loss();
  const double den_term = (max_val - min_val) / (max_val + min_val);
  if (std::pow(den_term, 2) < zero_tol) return identity_loss();
  return ContinuousLoss{min_val, max_val, &continuous_kernels::ratio};
}

namespace gwise_kernels {

double frechet_nominal(const int* values, int g, int C) {
  std::vector<int> counts(static_cast<std::size_t>(C), 0);
  int mode_count = 0;
  for (int i = 0; i < g; ++i) {
    const int count = ++counts[static_cast<std::size_t>(values[i])];
    if (count > mode_count) mode_count = count;
  }
  return static_cast<double>(g - mode_count) / static_cast<double>(g);
}

double hubert_categorical(const int* values, int g, int /*C*/) {
  for (int i = 1; i < g; ++i) {
    if (values[i] != values[0]) return 1.0;
  }
  return 0.0;
}

double frechet_absolute(const double* values, int g) {
  std::vector<double> work(values, values + g);
  std::sort(work.begin(), work.end());
  const double median = (g % 2 == 1)
                            ? work[static_cast<std::size_t>(g / 2)]
                            : 0.5 * (work[static_cast<std::size_t>(g / 2 - 1)]
                                     + work[static_cast<std::size_t>(g / 2)]);
  double acc = 0.0;
  for (int i = 0; i < g; ++i) acc += std::abs(values[i] - median);
  return acc / static_cast<double>(g);
}

double frechet_quadratic(const double* values, int g) {
  double mean = 0.0;
  for (int i = 0; i < g; ++i) mean += values[i];
  mean /= static_cast<double>(g);
  double acc = 0.0;
  for (int i = 0; i < g; ++i) acc += std::pow(values[i] - mean, 2);
  return acc / static_cast<double>(g);
}

double hubert_continuous(const double* values, int g) {
  for (int i = 1; i < g; ++i) {
    if (std::abs(values[i] - values[0]) >= zero_tol) return 1.0;
  }
  return 0.0;
}

}  // namespace gwise_kernels

Result<GwiseCategoricalDistance> frechet_nominal_distance(int C) {
  if (C <= 0) return std::unexpected(Error::invalid_argument);
  return GwiseCategoricalDistance{C, &gwise_kernels::frechet_nominal};
}

Result<GwiseCategoricalDistance> hubert_categorical_distance(int C) {
  if (C <= 0) return std::unexpected(Error::invalid_argument);
  return GwiseCategoricalDistance{C, &gwise_kernels::hubert_categorical};
}

Result<GwiseContinuousDistance> frechet_absolute_distance() {
  return GwiseContinuousDistance{&gwise_kernels::frechet_absolute};
}

Result<GwiseContinuousDistance> frechet_quadratic_distance() {
  return GwiseContinuousDistance{&gwise_kernels::frechet_quadratic};
}

Result<GwiseContinuousDistance> hubert_continuous_distance() {
  return GwiseContinuousDistance{&gwise_kernels::hubert_continuous};
}

}  // namespace misskappa::loss
