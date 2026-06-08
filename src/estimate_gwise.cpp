#include "misskappa/estimate.hpp"

#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

namespace misskappa {

namespace {

constexpr double zero_tol = 1e-9;

bool checked_power(int base, int exponent, std::int64_t limit, std::int64_t& out) {
  if (base < 0 || exponent < 0 || limit < 1) return false;
  std::int64_t value = 1;
  for (int i = 0; i < exponent; ++i) {
    if (base != 0 && value > limit / base) return false;
    value *= base;
  }
  out = value;
  return true;
}

void combinations_rec(
    int R, int g, int start, std::vector<int>& current,
    std::vector<std::vector<int>>& out) {
  if (static_cast<int>(current.size()) == g) {
    out.push_back(current);
    return;
  }
  const int need = g - static_cast<int>(current.size());
  for (int r = start; r <= R - need; ++r) {
    current.push_back(r);
    combinations_rec(R, g, r + 1, current, out);
    current.pop_back();
  }
}

std::vector<std::vector<int>> combinations(int R, int g) {
  std::vector<std::vector<int>> out;
  std::vector<int> current;
  combinations_rec(R, g, 0, current, out);
  return out;
}

void rater_tuples_rec(
    int R, int g, std::vector<int>& current,
    std::vector<std::vector<int>>& out) {
  if (static_cast<int>(current.size()) == g) {
    out.push_back(current);
    return;
  }
  for (int r = 0; r < R; ++r) {
    current.push_back(r);
    rater_tuples_rec(R, g, current, out);
    current.pop_back();
  }
}

std::vector<std::vector<int>> rater_tuples(int R, int g) {
  std::vector<std::vector<int>> out;
  std::vector<int> current;
  rater_tuples_rec(R, g, current, out);
  return out;
}

template <typename Visitor>
void visit_item_tuples_rec(
    int n, int g, int depth, std::vector<int>& current, Visitor& visitor) {
  if (depth == g) {
    visitor(current);
    return;
  }
  for (int i = 0; i < n; ++i) {
    current[static_cast<std::size_t>(depth)] = i;
    visit_item_tuples_rec(n, g, depth + 1, current, visitor);
  }
}

template <typename Visitor>
void visit_item_tuples(int n, int g, Visitor& visitor) {
  std::vector<int> current(static_cast<std::size_t>(g), 0);
  visit_item_tuples_rec(n, g, 0, current, visitor);
}

Estimation finish_estimation(
    RealVec&& d_values, double C_hat, double F_hat, int g,
    RealVec&& mu_C_projection, RealVec&& mu_F_projection) {
  const int n = static_cast<int>(d_values.size());
  const double D_hat = d_values.mean();

  RealVec estimates(2);
  estimates(0) = (C_hat > zero_tol)
                     ? 1.0 - D_hat / C_hat
                     : std::numeric_limits<double>::quiet_NaN();
  estimates(1) = (F_hat > zero_tol)
                     ? 1.0 - D_hat / F_hat
                     : std::numeric_limits<double>::quiet_NaN();

  RealMat phi(n, 3);
  phi.col(0) = d_values.array() - D_hat;
  phi.col(1) = mu_C_projection.array() - static_cast<double>(g) * C_hat;
  phi.col(2) = mu_F_projection.array() - static_cast<double>(g) * F_hat;

  RealMat Gamma_hat = (phi.transpose() * phi) / static_cast<double>(n);

  RealMat J = RealMat::Zero(2, 3);
  if (C_hat > zero_tol) {
    J(0, 0) = -1.0 / C_hat;
    J(0, 1) = D_hat / (C_hat * C_hat);
  }
  if (F_hat > zero_tol) {
    J(1, 0) = -1.0 / F_hat;
    J(1, 2) = D_hat / (F_hat * F_hat);
  }

  RealMat vcov = (J * Gamma_hat * J.transpose()) / static_cast<double>(n);
  RealMat psi = build_psi_from_phi(phi, J);
  return Estimation{std::move(estimates), std::move(vcov), std::move(psi)};
}

RealMat categorical_margins(IntMatView ratings, int C) {
  const int n = static_cast<int>(ratings.rows());
  const int R = static_cast<int>(ratings.cols());
  RealMat probs = RealMat::Zero(R, C);
  for (int i = 0; i < n; ++i) {
    for (int r = 0; r < R; ++r) {
      probs(r, ratings(i, r)) += 1.0;
    }
  }
  probs /= static_cast<double>(n);
  return probs;
}

template <typename ProbAtPos>
double categorical_expectation(
    loss::GwiseCategoricalDistance distance, int g, ProbAtPos&& prob_at_pos) {
  std::vector<int> values(static_cast<std::size_t>(g), 0);
  double acc = 0.0;
  auto visitor = [&](const std::vector<int>& cats) {
    double prob = 1.0;
    for (int pos = 0; pos < g; ++pos) {
      const int cat = cats[static_cast<std::size_t>(pos)];
      values[static_cast<std::size_t>(pos)] = cat;
      prob *= prob_at_pos(pos, cat);
    }
    acc += prob * distance.compute(values.data(), g, distance.C);
  };
  visit_item_tuples(distance.C, g, visitor);
  return acc;
}

template <typename ProbAtPos>
double categorical_expectation_fixed(
    loss::GwiseCategoricalDistance distance, int g, int fixed_pos,
    int fixed_value, ProbAtPos&& prob_at_pos) {
  std::vector<int> values(static_cast<std::size_t>(g), 0);
  values[static_cast<std::size_t>(fixed_pos)] = fixed_value;
  double acc = 0.0;
  auto visitor = [&](const std::vector<int>& cats) {
    int cursor = 0;
    double prob = 1.0;
    for (int pos = 0; pos < g; ++pos) {
      if (pos == fixed_pos) continue;
      const int cat = cats[static_cast<std::size_t>(cursor++)];
      values[static_cast<std::size_t>(pos)] = cat;
      prob *= prob_at_pos(pos, cat);
    }
    acc += prob * distance.compute(values.data(), g, distance.C);
  };
  visit_item_tuples(distance.C, g - 1, visitor);
  return acc;
}

}  // namespace

Result<Estimation> estimate_gwise(
    IntMatView ratings, loss::GwiseCategoricalDistance distance,
    GwiseOptions opts) {
  const int n = static_cast<int>(ratings.rows());
  const int R = static_cast<int>(ratings.cols());
  const int g = (opts.g <= 0) ? R : opts.g;
  if (n < 1 || R < 2 || g < 2 || g > R) return misskappa::unexpected(Error::invalid_argument);
  if (distance.C <= 0 || distance.compute == nullptr) {
    return misskappa::unexpected(Error::invalid_argument);
  }

  for (Eigen::Index i = 0; i < ratings.rows(); ++i) {
    for (Eigen::Index j = 0; j < ratings.cols(); ++j) {
      const int x = ratings(i, j);
      if (x < 0 || x >= distance.C) return misskappa::unexpected(Error::invalid_argument);
    }
  }

  std::int64_t category_tuples = 0;
  std::int64_t category_projection_tuples = 0;
  if (!checked_power(distance.C, g, opts.max_chance_tuples, category_tuples)
      || !checked_power(distance.C, g - 1, opts.max_chance_tuples, category_projection_tuples)) {
    return misskappa::unexpected(Error::not_supported);
  }

  const auto c_raters = combinations(R, g);
  if (c_raters.empty()) return misskappa::unexpected(Error::invalid_argument);

  std::vector<int> values(static_cast<std::size_t>(g), 0);
  RealVec d_values = RealVec::Zero(n);
  for (int i = 0; i < n; ++i) {
    double acc = 0.0;
    for (const auto& raters : c_raters) {
      for (int pos = 0; pos < g; ++pos) values[static_cast<std::size_t>(pos)] = ratings(i, raters[pos]);
      acc += distance.compute(values.data(), g, distance.C);
    }
    d_values(i) = acc / static_cast<double>(c_raters.size());
  }

  const RealMat probs = categorical_margins(ratings, distance.C);
  const RealVec pooled = probs.colwise().mean().transpose();

  double C_hat = 0.0;
  for (const auto& raters : c_raters) {
    C_hat += categorical_expectation(
        distance, g, [&](int pos, int cat) { return probs(raters[pos], cat); });
  }
  C_hat /= static_cast<double>(c_raters.size());

  const double F_hat = categorical_expectation(
      distance, g, [&](int /*pos*/, int cat) { return pooled(cat); });

  RealVec mu_C_projection = RealVec::Zero(n);
  RealVec mu_F_projection = RealVec::Zero(n);
  for (int i = 0; i < n; ++i) {
    double c_projection = 0.0;
    double f_projection = 0.0;
    for (int pos = 0; pos < g; ++pos) {
      double c_at_pos = 0.0;
      for (const auto& raters : c_raters) {
        c_at_pos += categorical_expectation_fixed(
            distance, g, pos, ratings(i, raters[pos]),
            [&](int other_pos, int cat) { return probs(raters[other_pos], cat); });
      }
      c_projection += c_at_pos / static_cast<double>(c_raters.size());

      double f_at_pos = 0.0;
      for (int r = 0; r < R; ++r) {
        f_at_pos += categorical_expectation_fixed(
            distance, g, pos, ratings(i, r),
            [&](int /*other_pos*/, int cat) { return pooled(cat); });
      }
      f_projection += f_at_pos / static_cast<double>(R);
    }
    mu_C_projection(i) = c_projection;
    mu_F_projection(i) = f_projection;
  }

  return finish_estimation(std::move(d_values), C_hat, F_hat, g,
                           std::move(mu_C_projection), std::move(mu_F_projection));
}

Result<Estimation> estimate_gwise_continuous(
    RealMatView ratings, loss::GwiseContinuousDistance distance,
    GwiseOptions opts) {
  const int n = static_cast<int>(ratings.rows());
  const int R = static_cast<int>(ratings.cols());
  const int g = (opts.g <= 0) ? R : opts.g;
  if (n < 1 || R < 2 || g < 2 || g > R) return misskappa::unexpected(Error::invalid_argument);
  if (distance.compute == nullptr) return misskappa::unexpected(Error::invalid_argument);

  for (Eigen::Index i = 0; i < ratings.rows(); ++i) {
    for (Eigen::Index j = 0; j < ratings.cols(); ++j) {
      if (!std::isfinite(ratings(i, j))) return misskappa::unexpected(Error::invalid_argument);
    }
  }

  std::int64_t n_g = 0;
  std::int64_t n_g_minus_1 = 0;
  if (!checked_power(n, g, opts.max_chance_tuples, n_g)
      || !checked_power(n, g - 1, opts.max_chance_tuples, n_g_minus_1)) {
    return misskappa::unexpected(Error::not_supported);
  }

  const auto c_raters = combinations(R, g);
  const auto f_raters = rater_tuples(R, g);
  if (c_raters.empty() || f_raters.empty()) return misskappa::unexpected(Error::invalid_argument);

  std::vector<double> values(static_cast<std::size_t>(g), 0.0);
  RealVec d_values = RealVec::Zero(n);
  for (int i = 0; i < n; ++i) {
    double acc = 0.0;
    for (const auto& raters : c_raters) {
      for (int pos = 0; pos < g; ++pos) values[static_cast<std::size_t>(pos)] = ratings(i, raters[pos]);
      acc += distance.compute(values.data(), g);
    }
    d_values(i) = acc / static_cast<double>(c_raters.size());
  }

  RealVec mu_C_sum = RealVec::Zero(n);
  RealVec mu_F_sum = RealVec::Zero(n);
  double C_total = 0.0;
  double F_total = 0.0;

  auto eval_C = [&](const std::vector<int>& items) {
    double acc = 0.0;
    for (const auto& raters : c_raters) {
      for (int pos = 0; pos < g; ++pos) {
        values[static_cast<std::size_t>(pos)] = ratings(items[pos], raters[pos]);
      }
      acc += distance.compute(values.data(), g);
    }
    return acc / static_cast<double>(c_raters.size());
  };

  auto eval_F = [&](const std::vector<int>& items) {
    double acc = 0.0;
    for (const auto& raters : f_raters) {
      for (int pos = 0; pos < g; ++pos) {
        values[static_cast<std::size_t>(pos)] = ratings(items[pos], raters[pos]);
      }
      acc += distance.compute(values.data(), g);
    }
    return acc / static_cast<double>(f_raters.size());
  };

  auto visitor = [&](const std::vector<int>& items) {
    const double c = eval_C(items);
    const double f = eval_F(items);
    C_total += c;
    F_total += f;
    for (int pos = 0; pos < g; ++pos) {
      mu_C_sum(items[pos]) += c;
      mu_F_sum(items[pos]) += f;
    }
  };
  visit_item_tuples(n, g, visitor);

  const double C_hat = C_total / static_cast<double>(n_g);
  const double F_hat = F_total / static_cast<double>(n_g);
  RealVec mu_C_projection = mu_C_sum / static_cast<double>(n_g_minus_1);
  RealVec mu_F_projection = mu_F_sum / static_cast<double>(n_g_minus_1);
  return finish_estimation(std::move(d_values), C_hat, F_hat, g,
                           std::move(mu_C_projection), std::move(mu_F_projection));
}

}  // namespace misskappa
