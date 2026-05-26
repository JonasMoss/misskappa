// FIML / EM estimator for categorical raw ratings under MAR.
//
// Algorithm preserved from the legacy emdiscrete::run_em and
// misskappa::kappaml::kappa pipeline in dev/legacy/misskappa/src/, ported
// onto Eigen + Result<T> with no exceptions. Only the "raw" categorical
// path is implemented; counts-format is out of Phase 1 scope.

#include "misskappa/estimate.hpp"

#include <Eigen/QR>
#include <algorithm>
#include <cmath>
#include <limits>
#include <cstddef>
#include <cstdint>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

namespace misskappa {

namespace {

constexpr double zero_tol = 1e-9;
constexpr double singular_tol = 1e-12;

// --- Internal EM types ---

struct EmInput {
  std::size_t n_total_patterns = 0;
  int c = 0;
  int R = 0;
  // Each "group" is one observed pattern (unique row of the input). Within a
  // group, completion_indices[group_offsets[g] .. +group_n_completions[g])
  // lists the ranks (in the C^R completion space) of every full pattern
  // compatible with the observed row.
  std::vector<std::uint32_t> group_n_subjects;
  std::vector<std::uint32_t> group_offsets;
  std::vector<std::uint32_t> group_n_completions;
  std::vector<std::uint64_t> completion_indices;
};

struct EmRunResult {
  Eigen::VectorXd theta_hat;    // pruned, normalised
  RealMat vcov;                 // covariance of theta_hat
  std::vector<std::uint64_t> pattern_indices;  // ranks of surviving patterns
  int c = 0;
  int R = 0;
  int iterations = 0;
  bool converged = false;
  std::size_t n_subjects = 0;
};

// --- Encoding / decoding of categorical patterns as base-c integers ---

std::uint64_t rank_tuple(const std::vector<int>& row, int c) {
  std::uint64_t rank = 0;
  std::uint64_t power = 1;
  for (int i = static_cast<int>(row.size()) - 1; i >= 0; --i) {
    rank += static_cast<std::uint64_t>(row[i]) * power;
    if (i > 0) power *= static_cast<std::uint64_t>(c);
  }
  return rank;
}

std::vector<int> unrank_tuple(std::uint64_t rank, int R, int c) {
  std::vector<int> rt(R, 0);
  std::uint64_t temp = rank;
  for (int i = R - 1; i >= 0; --i) {
    rt[i] = static_cast<int>(temp % static_cast<std::uint64_t>(c));
    temp /= static_cast<std::uint64_t>(c);
  }
  return rt;
}

void enumerate_completions(
    std::vector<int>& row, const std::vector<int>& missing_positions,
    std::size_t k, int c, std::vector<std::vector<int>>& out) {
  if (k == missing_positions.size()) {
    out.push_back(row);
    return;
  }
  const int pos = missing_positions[k];
  for (int cat = 0; cat < c; ++cat) {
    row[pos] = cat;
    enumerate_completions(row, missing_positions, k + 1, c, out);
  }
}

// --- Preprocessing: build EM input from raw n x R matrix ---

Result<EmInput> preprocess_raw(IntMatView ratings, int c) {
  EmInput in;
  in.R = static_cast<int>(ratings.cols());
  in.c = c;

  // Validate ratings: must be na_code or in [0, c-1].
  for (Eigen::Index i = 0; i < ratings.rows(); ++i) {
    for (Eigen::Index j = 0; j < ratings.cols(); ++j) {
      const int v = ratings(i, j);
      if (v == na_code) continue;
      if (v < 0 || v >= c) return std::unexpected(Error::invalid_argument);
    }
  }

  // Bucket subjects by observed pattern (as a string key for stable ordering).
  std::map<std::string, std::pair<std::vector<int>, std::uint32_t>> by_key;
  std::vector<std::string> pattern_order;
  pattern_order.reserve(static_cast<std::size_t>(ratings.rows()));
  for (Eigen::Index i = 0; i < ratings.rows(); ++i) {
    std::string key;
    std::vector<int> row(static_cast<std::size_t>(in.R));
    for (int j = 0; j < in.R; ++j) {
      row[static_cast<std::size_t>(j)] = ratings(i, j);
      key += std::to_string(ratings(i, j));
      key += ',';
    }
    auto it = by_key.find(key);
    if (it == by_key.end()) {
      by_key.emplace(key, std::make_pair(std::move(row), std::uint32_t{1}));
      pattern_order.push_back(std::move(key));
    } else {
      it->second.second += 1;
    }
  }

  // For each unique observed pattern, enumerate compatible completions.
  std::uint32_t cum_offset = 0;
  for (const auto& key : pattern_order) {
    const auto& entry = by_key.at(key);
    const std::vector<int>& observed = entry.first;
    const std::uint32_t count = entry.second;

    std::vector<int> missing_positions;
    missing_positions.reserve(observed.size());
    for (std::size_t j = 0; j < observed.size(); ++j) {
      if (observed[j] == na_code) missing_positions.push_back(static_cast<int>(j));
    }

    std::vector<std::vector<int>> completions;
    if (missing_positions.empty()) {
      completions.push_back(observed);
    } else {
      std::vector<int> scratch = observed;
      enumerate_completions(scratch, missing_positions, 0, c, completions);
    }

    in.group_n_subjects.push_back(count);
    in.group_offsets.push_back(cum_offset);
    in.group_n_completions.push_back(static_cast<std::uint32_t>(completions.size()));
    for (const auto& comp : completions) {
      in.completion_indices.push_back(rank_tuple(comp, c));
    }
    cum_offset += static_cast<std::uint32_t>(completions.size());
  }

  // Total pattern space: c^R.
  std::uint64_t total = 1;
  for (int k = 0; k < in.R; ++k) total *= static_cast<std::uint64_t>(c);
  in.n_total_patterns = static_cast<std::size_t>(total);

  return in;
}

// --- EM core loop ---

Eigen::VectorXd initialise_theta(const EmInput& in, double start_alpha) {
  Eigen::VectorXd theta = Eigen::VectorXd::Constant(
      static_cast<Eigen::Index>(in.n_total_patterns), start_alpha);
  if (in.n_total_patterns == 0) return theta;
  for (std::size_t g = 0; g < in.group_n_subjects.size(); ++g) {
    const std::uint32_t n_subs = in.group_n_subjects[g];
    const std::uint32_t n_comps = in.group_n_completions[g];
    if (n_comps == 0) continue;
    const std::uint32_t offset = in.group_offsets[g];
    const double weight = static_cast<double>(n_subs) / static_cast<double>(n_comps);
    for (std::uint32_t k = 0; k < n_comps; ++k) {
      const std::uint64_t idx = in.completion_indices[offset + k];
      theta(static_cast<Eigen::Index>(idx)) += weight;
    }
  }
  const double sum = theta.sum();
  if (sum < zero_tol) return Eigen::VectorXd::Zero(theta.size());
  theta /= sum;
  return theta;
}

// In-place EM iteration; returns final iterations + convergence flag.
struct EmIterStatus {
  int iterations = 0;
  bool converged = false;
};

EmIterStatus run_em_iterations(
    Eigen::VectorXd& theta, const EmInput& in, EmOptions opts) {
  EmIterStatus status;
  std::size_t n_total_subjects = 0;
  for (std::uint32_t s : in.group_n_subjects) n_total_subjects += s;
  if (n_total_subjects == 0) {
    status.converged = true;
    return status;
  }
  Eigen::VectorXd expected = Eigen::VectorXd::Zero(theta.size());
  for (int it = 0; it < opts.max_iter; ++it) {
    status.iterations = it + 1;
    expected.setZero();
    for (std::size_t g = 0; g < in.group_n_subjects.size(); ++g) {
      const std::uint32_t n_comps = in.group_n_completions[g];
      if (n_comps == 0) continue;
      const std::uint32_t offset = in.group_offsets[g];
      // posterior over completions, proportional to theta on those completions.
      double sum_w = 0.0;
      for (std::uint32_t k = 0; k < n_comps; ++k) {
        sum_w += theta(static_cast<Eigen::Index>(in.completion_indices[offset + k]));
      }
      if (sum_w <= singular_tol) continue;
      const double scale = static_cast<double>(in.group_n_subjects[g]) / sum_w;
      for (std::uint32_t k = 0; k < n_comps; ++k) {
        const Eigen::Index idx = static_cast<Eigen::Index>(in.completion_indices[offset + k]);
        expected(idx) += scale * theta(idx);
      }
    }
    Eigen::VectorXd new_theta = expected / static_cast<double>(n_total_subjects);
    const double max_change = (new_theta - theta).cwiseAbs().maxCoeff();
    if (max_change < opts.tol) {
      theta = std::move(new_theta);
      status.converged = true;
      return status;
    }
    theta = std::move(new_theta);
  }
  return status;
}

// Louis observed-information variance for the pruned theta.
RealMat em_variance(
    const Eigen::VectorXd& theta_pruned,
    const std::vector<std::uint64_t>& pruned_ranks,
    const EmInput& in) {
  const Eigen::Index n_final = theta_pruned.size();
  if (n_final <= 1) return RealMat::Zero(n_final, n_final);

  // Identify the reference index (largest probability).
  Eigen::Index ref = 0;
  for (Eigen::Index i = 1; i < n_final; ++i) {
    if (theta_pruned(i) > theta_pruned(ref)) ref = i;
  }

  std::unordered_map<std::uint64_t, Eigen::Index> rank_to_pruned;
  rank_to_pruned.reserve(pruned_ranks.size());
  for (Eigen::Index i = 0; i < n_final; ++i) {
    rank_to_pruned.emplace(pruned_ranks[static_cast<std::size_t>(i)], i);
  }

  RealMat info_star = RealMat::Zero(n_final - 1, n_final - 1);

  for (std::size_t g = 0; g < in.group_n_subjects.size(); ++g) {
    const std::uint32_t n_comps = in.group_n_completions[g];
    const std::uint32_t offset = in.group_offsets[g];

    double sum_theta_subset = 0.0;
    bool has_pruned_completion = false;
    for (std::uint32_t k = 0; k < n_comps; ++k) {
      const std::uint64_t rank = in.completion_indices[offset + k];
      auto it = rank_to_pruned.find(rank);
      if (it != rank_to_pruned.end()) {
        has_pruned_completion = true;
        sum_theta_subset += theta_pruned(it->second);
      }
    }
    if (!has_pruned_completion || sum_theta_subset <= singular_tol) continue;

    // Score for theta_pruned: posterior / theta_pruned. In raw FIML this is
    // constant across all compatible retained completions.
    const double retained_score = 1.0 / sum_theta_subset;
    double ref_score = 0.0;
    Eigen::VectorXd s_reduced = Eigen::VectorXd::Zero(n_final - 1);
    for (std::uint32_t k = 0; k < n_comps; ++k) {
      const std::uint64_t rank = in.completion_indices[offset + k];
      auto it = rank_to_pruned.find(rank);
      if (it == rank_to_pruned.end()) continue;
      const Eigen::Index idx = it->second;
      if (idx == ref) {
        ref_score = retained_score;
      } else {
        const Eigen::Index reduced_idx = (idx < ref) ? idx : idx - 1;
        s_reduced(reduced_idx) = retained_score;
      }
    }
    if (ref_score != 0.0) s_reduced.array() -= ref_score;

    info_star.noalias() += static_cast<double>(in.group_n_subjects[g])
                          * s_reduced * s_reduced.transpose();
  }

  // Symmetrise (numerical safety) and invert.
  info_star = 0.5 * (info_star + info_star.transpose());
  Eigen::CompleteOrthogonalDecomposition<RealMat> cod(info_star);
  RealMat var_star = cod.pseudoInverse();

  // Expand back to the full n_final dimension via the constraint that theta
  // sums to 1: theta_ref = 1 - sum of others.
  RealMat J = RealMat::Zero(n_final, n_final - 1);
  Eigen::Index w = 0;
  for (Eigen::Index i = 0; i < n_final; ++i) {
    if (i == ref) continue;
    J(i, w++) = 1.0;
  }
  for (Eigen::Index j = 0; j < n_final - 1; ++j) J(ref, j) = -1.0;

  return J * var_star * J.transpose();
}

Result<EmRunResult> run_em(IntMatView ratings, int c, EmOptions opts) {
  if (ratings.rows() == 0 || ratings.cols() == 0) {
    return std::unexpected(Error::invalid_argument);
  }
  auto in_res = preprocess_raw(ratings, c);
  if (!in_res) return std::unexpected(in_res.error());
  const EmInput& in = *in_res;

  EmRunResult out;
  out.c = c;
  out.R = in.R;
  for (std::uint32_t s : in.group_n_subjects) out.n_subjects += s;

  if (in.n_total_patterns == 0) {
    out.converged = true;
    return out;
  }
  Eigen::VectorXd theta = initialise_theta(in, opts.start_alpha);
  const EmIterStatus status = run_em_iterations(theta, in, opts);
  out.iterations = status.iterations;
  out.converged = status.converged;
  if (!out.converged) return std::unexpected(Error::not_converged);

  // Prune patterns with negligible probability mass.
  std::vector<std::uint64_t> pruned_ranks;
  std::vector<double> pruned_theta;
  pruned_ranks.reserve(static_cast<std::size_t>(theta.size()));
  pruned_theta.reserve(static_cast<std::size_t>(theta.size()));
  for (Eigen::Index i = 0; i < theta.size(); ++i) {
    if (theta(i) > opts.prune_tol) {
      pruned_ranks.push_back(static_cast<std::uint64_t>(i));
      pruned_theta.push_back(theta(i));
    }
  }
  if (pruned_ranks.empty()) return std::unexpected(Error::numerical_error);

  Eigen::VectorXd theta_pruned = Eigen::Map<Eigen::VectorXd>(
      pruned_theta.data(), static_cast<Eigen::Index>(pruned_theta.size()));
  const double s = theta_pruned.sum();
  if (s < zero_tol) return std::unexpected(Error::numerical_error);
  theta_pruned /= s;

  out.theta_hat = std::move(theta_pruned);
  out.pattern_indices = std::move(pruned_ranks);
  out.vcov = em_variance(out.theta_hat, out.pattern_indices, in);
  return out;
}

// --- Mapping from theta to (Conger, Fleiss, Brennan-Prediger) ---

struct KappaMap {
  Eigen::VectorXd d_vec;       // length n_final
  RealMat Qed_conger;          // n_final x n_final
  RealMat Qed_fleiss;          // n_final x n_final
  double d_bp = 0.0;
};

KappaMap build_kappa_map(const EmRunResult& em, RealMatView weights) {
  // Use the AGREEMENT-matrix convention on input; convert internally.
  const int C = static_cast<int>(weights.rows());
  const RealMat L = RealMat::Constant(C, C, 1.0) - weights;

  const int R = em.R;
  const Eigen::Index n_final = em.theta_hat.size();
  const double n_pairs = static_cast<double>(R) * (R - 1.0) / 2.0;
  const double n_pairs_full = static_cast<double>(R) * (R - 1.0);

  // Per-pattern average within-subject disagreement.
  KappaMap m;
  m.d_vec = Eigen::VectorXd::Zero(n_final);
  std::vector<std::vector<int>> pattern_matrix;
  pattern_matrix.reserve(static_cast<std::size_t>(n_final));
  for (Eigen::Index j = 0; j < n_final; ++j) {
    auto row = unrank_tuple(em.pattern_indices[static_cast<std::size_t>(j)], R, em.c);
    double sum = 0.0;
    for (int c1 = 0; c1 < R - 1; ++c1) {
      for (int c2 = c1 + 1; c2 < R; ++c2) {
        sum += L(row[static_cast<std::size_t>(c1)], row[static_cast<std::size_t>(c2)]);
      }
    }
    m.d_vec(j) = sum / n_pairs;
    pattern_matrix.push_back(std::move(row));
  }

  // Conger: Qed_conger[u, v] = (1 / n_pairs) * sum_{c1 < c2} L(x_u[c1], x_v[c2]).
  m.Qed_conger = RealMat::Zero(n_final, n_final);
  for (Eigen::Index u = 0; u < n_final; ++u) {
    for (Eigen::Index v = 0; v < n_final; ++v) {
      double sum = 0.0;
      for (int c1 = 0; c1 < R - 1; ++c1) {
        for (int c2 = c1 + 1; c2 < R; ++c2) {
          sum += L(pattern_matrix[static_cast<std::size_t>(u)][static_cast<std::size_t>(c1)],
                   pattern_matrix[static_cast<std::size_t>(v)][static_cast<std::size_t>(c2)]);
        }
      }
      m.Qed_conger(u, v) = sum / n_pairs;
    }
  }

  // Fleiss: M_fleiss[k, u] = (1 / R) * (# of raters in pattern u giving category k).
  // Qed_fleiss = M^T L M.
  RealMat M_fleiss = RealMat::Zero(em.c, n_final);
  for (Eigen::Index u = 0; u < n_final; ++u) {
    for (int j = 0; j < R; ++j) {
      M_fleiss(pattern_matrix[static_cast<std::size_t>(u)][static_cast<std::size_t>(j)], u) += 1.0;
    }
  }
  M_fleiss /= static_cast<double>(R);
  m.Qed_fleiss = M_fleiss.transpose() * L * M_fleiss;

  // Brennan-Prediger.
  m.d_bp = L.sum() / (static_cast<double>(em.c) * static_cast<double>(em.c));

  (void)n_pairs_full;
  return m;
}

}  // namespace

Result<Estimation> estimate_fiml(
    IntMatView ratings, RealMatView weights, EmOptions opts) {
  const int n = static_cast<int>(ratings.rows());
  const int R = static_cast<int>(ratings.cols());
  if (n < 1) return std::unexpected(Error::invalid_argument);
  if (R < 2) return std::unexpected(Error::invalid_argument);
  const int C = static_cast<int>(weights.rows());
  if (weights.cols() != C) return std::unexpected(Error::dimension_mismatch);
  if (C < 1) return std::unexpected(Error::invalid_argument);

  auto em = run_em(ratings, C, opts);
  if (!em) return std::unexpected(em.error());
  if (em->theta_hat.size() == 0) return std::unexpected(Error::numerical_error);

  const KappaMap m = build_kappa_map(*em, weights);

  const Eigen::VectorXd& theta = em->theta_hat;
  const double pd  = m.d_vec.dot(theta);
  const double pec = (theta.transpose() * m.Qed_conger * theta).value();
  const double pef = (theta.transpose() * m.Qed_fleiss * theta).value();

  RealVec estimates(3);
  estimates(0) = (std::abs(pec) > singular_tol)
                     ? 1.0 - pd / pec
                     : std::numeric_limits<double>::quiet_NaN();
  estimates(1) = (std::abs(pef) > singular_tol)
                     ? 1.0 - pd / pef
                     : std::numeric_limits<double>::quiet_NaN();
  estimates(2) = (std::abs(m.d_bp) > singular_tol)
                     ? 1.0 - pd / m.d_bp
                     : std::numeric_limits<double>::quiet_NaN();

  // Jacobian of (kappa_C, kappa_F, kappa_BP) wrt theta.
  RealMat jacobian = RealMat::Zero(3, theta.size());
  const Eigen::VectorXd grad_pd = m.d_vec;
  if (std::abs(pec) > singular_tol) {
    const RealMat sym = 0.5 * (m.Qed_conger + m.Qed_conger.transpose());
    const Eigen::VectorXd grad_pec = 2.0 * sym * theta;
    jacobian.row(0) = (-(grad_pd * pec - pd * grad_pec).transpose()) / (pec * pec);
  }
  if (std::abs(pef) > singular_tol) {
    const RealMat sym = 0.5 * (m.Qed_fleiss + m.Qed_fleiss.transpose());
    const Eigen::VectorXd grad_pef = 2.0 * sym * theta;
    jacobian.row(1) = (-(grad_pd * pef - pd * grad_pef).transpose()) / (pef * pef);
  }
  if (std::abs(m.d_bp) > singular_tol) {
    jacobian.row(2) = -grad_pd.transpose() / m.d_bp;
  }

  RealMat vcov = jacobian * em->vcov * jacobian.transpose();
  return Estimation{std::move(estimates), std::move(vcov)};
}

}  // namespace misskappa
