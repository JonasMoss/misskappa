// FIML / EM estimator for categorical raw ratings under MAR.
//
// Algorithm preserved from the legacy emdiscrete::run_em and
// misskappa::kappaml::kappa pipeline in dev/legacy/misskappa/src/, ported
// onto Eigen + Result<T> with no exceptions. Only the "raw" categorical
// path is implemented; counts-format is out of Phase 1 scope.

#include "misskappa/estimate.hpp"
#include "misskappa/diagnostics.hpp"

#include "detail_pattern_checks.hpp"
#include "detail_psd_inverse.hpp"
#include "prof.hpp"

#include <Eigen/Eigenvalues>
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
  std::size_t n_active = 0;
  int c = 0;
  int R = 0;
  // Each "group" is one observed pattern (unique row of the input). Within a
  // group, completion_indices[group_offsets[g] .. +group_n_completions[g])
  // lists positions in the ACTIVE set (union of all groups' completions);
  // active_ranks maps an active position back to its rank in the C^R space.
  std::vector<std::uint32_t> group_n_subjects;
  std::vector<std::uint32_t> group_offsets;
  std::vector<std::uint32_t> group_n_completions;
  std::vector<std::uint32_t> completion_indices;
  std::vector<std::uint64_t> active_ranks;
  std::vector<std::uint32_t> subject_groups;
};

struct EmRunResult {
  Eigen::VectorXd theta_hat;    // pruned, normalised
  RealMat vcov;                 // covariance of theta_hat
  std::vector<std::uint64_t> pattern_indices;  // ranks of surviving patterns
  std::vector<std::uint32_t> pruned_active;    // active positions of survivors
  int c = 0;
  int R = 0;
  int iterations = 0;
  bool converged = false;
  std::size_t n_subjects = 0;
};

struct LouisReducedInfo {
  RealMat info_star;
  RealMat jacobian;
  RealMat group_scores;
  Eigen::Index ref = 0;
};

struct FimlVarianceCache {
  RealMat theta_vcov;
  RealMat var_star;
  RealMat info_star;
  RealMat theta_jacobian;
  RealMat group_scores;
  double info_rcond = 5e-5;  // truncation used for var_star, reused by the
                             // null-fraction diagnostic so both agree on what
                             // counts as a null direction.
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

template <typename Visitor>
void visit_category_tuples_rec(
    int C, int g, int depth, std::vector<int>& current, Visitor& visitor) {
  if (depth == g) {
    visitor(current);
    return;
  }
  for (int cat = 0; cat < C; ++cat) {
    current[static_cast<std::size_t>(depth)] = cat;
    visit_category_tuples_rec(C, g, depth + 1, current, visitor);
  }
}

template <typename Visitor>
void visit_category_tuples(int C, int g, Visitor& visitor) {
  std::vector<int> current(static_cast<std::size_t>(g), 0);
  visit_category_tuples_rec(C, g, 0, current, visitor);
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
      if (v < 0 || v >= c) return misskappa::unexpected(Error::invalid_argument);
    }
  }

  // Bucket subjects by observed pattern (as a string key for stable ordering).
  std::map<std::string, std::pair<std::vector<int>, std::uint32_t>> by_key;
  std::vector<std::string> pattern_order;
  std::vector<std::string> subject_keys;
  pattern_order.reserve(static_cast<std::size_t>(ratings.rows()));
  subject_keys.reserve(static_cast<std::size_t>(ratings.rows()));
  for (Eigen::Index i = 0; i < ratings.rows(); ++i) {
    std::string key;
    std::vector<int> row(static_cast<std::size_t>(in.R));
    for (int j = 0; j < in.R; ++j) {
      row[static_cast<std::size_t>(j)] = ratings(i, j);
      key += std::to_string(ratings(i, j));
      key += ',';
    }
    subject_keys.push_back(key);
    auto it = by_key.find(key);
    if (it == by_key.end()) {
      by_key.emplace(key, std::make_pair(std::move(row), std::uint32_t{1}));
      pattern_order.push_back(key);
    } else {
      it->second.second += 1;
    }
  }

  // For each unique observed pattern, enumerate compatible completions.
  std::uint32_t cum_offset = 0;
  std::vector<std::uint64_t> completion_ranks;
  std::unordered_map<std::string, std::uint32_t> key_to_group;
  key_to_group.reserve(pattern_order.size());
  for (std::size_t g = 0; g < pattern_order.size(); ++g) {
    const auto& key = pattern_order[g];
    key_to_group.emplace(key, static_cast<std::uint32_t>(g));
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
      completion_ranks.push_back(rank_tuple(comp, c));
    }
    cum_offset += static_cast<std::uint32_t>(completions.size());
  }
  in.subject_groups.reserve(subject_keys.size());
  for (const auto& key : subject_keys) {
    in.subject_groups.push_back(key_to_group.at(key));
  }

  // Compress the completion rank lists onto the active union: theta only ever
  // carries mass on cells reachable from at least one observed pattern (under
  // strict ML every other cell is exactly zero after the first M-step; under
  // flattening every other cell sits at the constant prior floor).
  in.active_ranks = completion_ranks;
  std::sort(in.active_ranks.begin(), in.active_ranks.end());
  in.active_ranks.erase(
      std::unique(in.active_ranks.begin(), in.active_ranks.end()),
      in.active_ranks.end());
  in.n_active = in.active_ranks.size();
  in.completion_indices.reserve(completion_ranks.size());
  for (std::uint64_t rk : completion_ranks) {
    const auto it = std::lower_bound(
        in.active_ranks.begin(), in.active_ranks.end(), rk);
    in.completion_indices.push_back(
        static_cast<std::uint32_t>(it - in.active_ranks.begin()));
  }

  // Total pattern space: c^R.
  std::uint64_t total = 1;
  for (int k = 0; k < in.R; ++k) total *= static_cast<std::uint64_t>(c);
  in.n_total_patterns = static_cast<std::size_t>(total);

  return in;
}

// --- EM core loop ---

// Returns theta over the ACTIVE set; `inactive_value` receives the shared
// value of every cell outside the active union (start_alpha / normaliser),
// so the result matches the full-table initialisation restricted to the
// active cells, with the same normalising constant.
Eigen::VectorXd initialise_theta(
    const EmInput& in, double start_alpha, double& inactive_value) {
  Eigen::VectorXd theta = Eigen::VectorXd::Constant(
      static_cast<Eigen::Index>(in.n_active), start_alpha);
  inactive_value = 0.0;
  if (in.n_active == 0) return theta;
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
  const double n_inactive =
      static_cast<double>(in.n_total_patterns) - static_cast<double>(in.n_active);
  const double sum = theta.sum() + start_alpha * n_inactive;
  if (sum < zero_tol) return Eigen::VectorXd::Zero(theta.size());
  theta /= sum;
  inactive_value = start_alpha / sum;
  return theta;
}

// In-place EM iteration; returns final iterations + convergence flag.
struct EmIterStatus {
  int iterations = 0;
  bool converged = false;
};

EmIterStatus run_em_iterations(
    Eigen::VectorXd& theta, double inactive_value,
    const EmInput& in, EmOptions opts) {
  EmIterStatus status;
  std::size_t n_total_subjects = 0;
  for (std::uint32_t s : in.group_n_subjects) n_total_subjects += s;
  if (n_total_subjects == 0) {
    status.converged = true;
    return status;
  }
  const double n_subj = static_cast<double>(n_total_subjects);
  // Dirichlet flattening: each cell gains delta = flatten / C^R pseudo-count
  // in the M-step, turning EM into posterior-mode iteration. Likelihood-flat
  // directions then contract only at rate ~ n / (n + flatten) per iteration.
  // Two adjustments keep that drift phase from tripping not_converged:
  // the iteration cap scales with (n + flatten) / flatten (bounded so
  // runaway fits stay finite), and the convergence tolerance is floored at
  // the per-iteration step a face-position error of face_tol produces —
  // the analytic-center position only needs resolving to face_tol, since
  // identified functionals are flat (to face-width order) along those
  // directions. flatten = 0 leaves strict ML and the user's settings
  // untouched.
  const double flatten = (opts.flatten > 0.0) ? opts.flatten : 0.0;
  const double delta =
      (flatten > 0.0 && theta.size() > 0)
          ? flatten / static_cast<double>(in.n_total_patterns)
          : 0.0;
  int iter_cap = opts.max_iter;
  double tol = opts.tol;
  if (flatten > 0.0) {
    const double scaled = 300.0 * (n_subj + flatten) / flatten;
    const double capped = std::min(scaled, 1.0e6);
    iter_cap = std::max(iter_cap, static_cast<int>(capped));
    constexpr double face_tol = 1e-4;
    tol = std::max(tol, face_tol * flatten / (n_subj + flatten));
  }
  // Cells outside the active union never receive expected mass, so they all
  // share one value: inactive_value at start, then delta / (n + flatten)
  // after every M-step (0 under strict ML). Track that as a scalar so the
  // dense state covers active cells only.
  const double inactive_next = delta / (n_subj + flatten);
  const bool any_inactive = in.n_total_patterns > in.n_active;
  Eigen::VectorXd next = Eigen::VectorXd::Zero(theta.size());
  std::uint32_t max_comps = 0;
  for (std::uint32_t nc : in.group_n_completions) max_comps = std::max(max_comps, nc);
  std::vector<double> wbuf(max_comps);
  for (int it = 0; it < iter_cap; ++it) {
    status.iterations = it + 1;
    next.setZero();
    const std::uint32_t* idx_base = in.completion_indices.data();
    for (std::size_t g = 0; g < in.group_n_subjects.size(); ++g) {
      const std::uint32_t n_comps = in.group_n_completions[g];
      if (n_comps == 0) continue;
      const std::uint32_t* idx = idx_base + in.group_offsets[g];
      // posterior over completions, proportional to theta on those completions.
      double sum_w = 0.0;
      for (std::uint32_t k = 0; k < n_comps; ++k) {
        wbuf[k] = theta(static_cast<Eigen::Index>(idx[k]));
        sum_w += wbuf[k];
      }
      if (sum_w <= singular_tol) continue;
      const double scale = static_cast<double>(in.group_n_subjects[g]) / sum_w;
      for (std::uint32_t k = 0; k < n_comps; ++k) {
        next(static_cast<Eigen::Index>(idx[k])) += scale * wbuf[k];
      }
    }
    next = (next.array() + delta).matrix() / (n_subj + flatten);
    double max_change = (next - theta).cwiseAbs().maxCoeff();
    if (any_inactive) {
      max_change = std::max(max_change, std::abs(inactive_next - inactive_value));
    }
    inactive_value = inactive_next;
    theta.swap(next);
    if (max_change < tol) {
      status.converged = true;
      return status;
    }
  }
  return status;
}

RealMat constraint_jacobian(Eigen::Index n_final, Eigen::Index ref) {
  if (n_final <= 1) return RealMat::Zero(n_final, 0);

  RealMat J = RealMat::Zero(n_final, n_final - 1);
  Eigen::Index w = 0;
  for (Eigen::Index i = 0; i < n_final; ++i) {
    if (i == ref) continue;
    J(i, w++) = 1.0;
  }
  for (Eigen::Index j = 0; j < n_final - 1; ++j) J(ref, j) = -1.0;
  return J;
}

LouisReducedInfo build_louis_reduced_info(
    const Eigen::VectorXd& theta_pruned,
    const std::vector<std::uint32_t>& pruned_active,
    const EmInput& in) {
  LouisReducedInfo out;
  const Eigen::Index n_final = theta_pruned.size();
  if (n_final <= 1) {
    out.info_star = RealMat::Zero(0, 0);
    out.jacobian = RealMat::Zero(n_final, 0);
    out.group_scores = RealMat::Zero(
        static_cast<Eigen::Index>(in.group_n_subjects.size()), 0);
    return out;
  }

  // Identify the reference index (largest probability).
  Eigen::Index ref = 0;
  for (Eigen::Index i = 1; i < n_final; ++i) {
    if (theta_pruned(i) > theta_pruned(ref)) ref = i;
  }
  out.ref = ref;

  std::vector<Eigen::Index> active_to_pruned(in.n_active, Eigen::Index{-1});
  for (Eigen::Index i = 0; i < n_final; ++i) {
    active_to_pruned[static_cast<std::size_t>(
        pruned_active[static_cast<std::size_t>(i)])] = i;
  }

  RealMat info_star = RealMat::Zero(n_final - 1, n_final - 1);
  RealMat group_scores = RealMat::Zero(
      static_cast<Eigen::Index>(in.group_n_subjects.size()), n_final - 1);

  for (std::size_t g = 0; g < in.group_n_subjects.size(); ++g) {
    const std::uint32_t n_comps = in.group_n_completions[g];
    const std::uint32_t offset = in.group_offsets[g];

    double sum_theta_subset = 0.0;
    bool has_pruned_completion = false;
    for (std::uint32_t k = 0; k < n_comps; ++k) {
      const Eigen::Index idx = active_to_pruned[in.completion_indices[offset + k]];
      if (idx >= 0) {
        has_pruned_completion = true;
        sum_theta_subset += theta_pruned(idx);
      }
    }
    if (!has_pruned_completion || sum_theta_subset <= singular_tol) continue;

    // Score for theta_pruned: posterior / theta_pruned. In raw FIML this is
    // constant across all compatible retained completions.
    const double retained_score = 1.0 / sum_theta_subset;
    double ref_score = 0.0;
    Eigen::VectorXd s_reduced = Eigen::VectorXd::Zero(n_final - 1);
    for (std::uint32_t k = 0; k < n_comps; ++k) {
      const Eigen::Index idx = active_to_pruned[in.completion_indices[offset + k]];
      if (idx < 0) continue;
      if (idx == ref) {
        ref_score = retained_score;
      } else {
        const Eigen::Index reduced_idx = (idx < ref) ? idx : idx - 1;
        s_reduced(reduced_idx) = retained_score;
      }
    }
    if (ref_score != 0.0) s_reduced.array() -= ref_score;
    group_scores.row(static_cast<Eigen::Index>(g)) = s_reduced.transpose();

    info_star.noalias() += static_cast<double>(in.group_n_subjects[g])
                          * s_reduced * s_reduced.transpose();
  }

  // Symmetrise for numerical safety before eigendecomposition.
  out.info_star = 0.5 * (info_star + info_star.transpose());
  out.jacobian = constraint_jacobian(n_final, ref);
  out.group_scores = std::move(group_scores);
  return out;
}

// Louis observed-information variance for the pruned theta.
RealMat em_variance(
    const Eigen::VectorXd& theta_pruned,
    const std::vector<std::uint32_t>& pruned_active,
    const EmInput& in,
    const EmOptions& opts) {
  const Eigen::Index n_final = theta_pruned.size();
  if (n_final <= 1) return RealMat::Zero(n_final, n_final);

  const LouisReducedInfo info = build_louis_reduced_info(theta_pruned, pruned_active, in);
  const RealMat var_star = detail::pseudo_inverse_psd(info.info_star, opts.info_rcond);
  return info.jacobian * var_star * info.jacobian.transpose();
}

FimlVarianceCache build_fiml_variance_cache(
    const EmRunResult& em, const EmInput& in, const EmOptions& opts) {
  FimlVarianceCache cache;
  const Eigen::Index n_final = em.theta_hat.size();
  cache.theta_vcov = RealMat::Zero(n_final, n_final);
  cache.var_star = RealMat::Zero(0, 0);
  cache.info_star = RealMat::Zero(0, 0);
  cache.theta_jacobian = RealMat::Zero(n_final, 0);
  cache.group_scores = RealMat::Zero(
      static_cast<Eigen::Index>(in.group_n_subjects.size()), 0);

  if (n_final <= 1) return cache;

  bench_prof::acc().n_final = static_cast<long>(n_final);
  LouisReducedInfo info;
  {
    bench_prof::Timer t(bench_prof::acc().louis);
    info = build_louis_reduced_info(em.theta_hat, em.pruned_active, in);
  }
  cache.info_rcond = opts.info_rcond;
  {
    bench_prof::Timer t(bench_prof::acc().pinv);
    cache.var_star = detail::pseudo_inverse_psd(info.info_star, opts.info_rcond);
  }
  cache.info_star = info.info_star;
  cache.theta_jacobian = info.jacobian;
  cache.group_scores = info.group_scores;
  {
    bench_prof::Timer t(bench_prof::acc().theta_vcov);
    cache.theta_vcov =
        cache.theta_jacobian * cache.var_star * cache.theta_jacobian.transpose();
  }
  return cache;
}

// Per-row fraction of the delta-method gradient lying in the truncated null
// space of the reduced Louis information (the directions pseudo_inverse_psd
// drops at the same rcond). This is a diagnostic, not a gate: when every
// rater pair is co-observed the coefficients are estimable functions of the
// identified pattern margins, so the saturated nuisance being flat does not
// break the point estimate; a fraction away from zero flags that the sample
// information is rank-deficient along directions the coefficient touches
// (sparse support), making the point estimate selection-dependent and the
// truncated SE optimistic about those directions.
Result<RealVec> reduced_gradient_null_fraction(
    const RealMat& info_star, const RealMat& jacobian_reduced, double rcond) {
  RealVec frac = RealVec::Zero(jacobian_reduced.rows());
  if (jacobian_reduced.cols() == 0 || jacobian_reduced.rows() == 0) return frac;
  if (info_star.rows() != jacobian_reduced.cols()
      || info_star.cols() != jacobian_reduced.cols()) {
    return misskappa::unexpected(Error::numerical_error);
  }

  Eigen::SelfAdjointEigenSolver<RealMat> es(0.5 * (info_star + info_star.transpose()));
  if (es.info() != Eigen::Success) {
    return misskappa::unexpected(Error::numerical_error);
  }

  const RealVec& evals = es.eigenvalues();
  const RealMat& evecs = es.eigenvectors();
  const double lambda_max = evals.maxCoeff();
  const double rc = (std::isfinite(rcond) && rcond > 0.0) ? rcond : 0.0;
  const double threshold =
      (std::isfinite(lambda_max) && lambda_max > 0.0) ? rc * lambda_max : 0.0;

  for (Eigen::Index row = 0; row < jacobian_reduced.rows(); ++row) {
    const Eigen::VectorXd grad = jacobian_reduced.row(row).transpose();
    const double norm_sq = grad.squaredNorm();
    if (!(norm_sq > 0.0) || !std::isfinite(norm_sq)) continue;
    double null_sq = 0.0;
    for (Eigen::Index k = 0; k < evals.size(); ++k) {
      if (evals(k) > threshold) continue;
      const double projection = evecs.col(k).dot(grad);
      if (!std::isfinite(projection)) {
        return misskappa::unexpected(Error::numerical_error);
      }
      null_sq += projection * projection;
    }
    frac(row) = std::sqrt(std::min(1.0, null_sq / norm_sq));
  }
  return frac;
}

Result<EmRunResult> run_em_preprocessed(
    const EmInput& in, int c, EmOptions opts, bool compute_vcov) {
  EmRunResult out;
  out.c = c;
  out.R = in.R;
  for (std::uint32_t s : in.group_n_subjects) out.n_subjects += s;

  if (in.n_total_patterns == 0) {
    out.converged = true;
    return out;
  }
  double inactive_value = 0.0;
  Eigen::VectorXd theta = initialise_theta(in, opts.start_alpha, inactive_value);
  bench_prof::acc().n_active = static_cast<long>(theta.size());
  bench_prof::acc().n_groups = static_cast<long>(in.group_n_subjects.size());
  EmIterStatus status;
  {
    bench_prof::Timer t(bench_prof::acc().em);
    status = run_em_iterations(theta, inactive_value, in, opts);
  }
  bench_prof::acc().em_iters += status.iterations;
  out.iterations = status.iterations;
  out.converged = status.converged;
  if (!out.converged) return misskappa::unexpected(Error::not_converged);

  // Prune patterns with negligible probability mass. With flattening every
  // cell sits at or above the prior floor delta / (n + flatten), so the
  // threshold moves just above that floor: cells whose mass is essentially
  // all prior (no data support) drop out and the variance machinery keeps
  // working on the data-supported cells only. The pruned prior mass totals
  // about flatten / (n + flatten), the same order as the flattening shift.
  double prune_threshold = opts.prune_tol;
  if (opts.flatten > 0.0 && in.n_total_patterns > 0) {
    const double floor =
        (opts.flatten / static_cast<double>(in.n_total_patterns))
        / (static_cast<double>(out.n_subjects) + opts.flatten);
    prune_threshold = std::max(prune_threshold, 2.0 * floor);
  }
  std::vector<std::uint64_t> pruned_ranks;
  std::vector<std::uint32_t> pruned_active;
  std::vector<double> pruned_theta;
  pruned_ranks.reserve(static_cast<std::size_t>(theta.size()));
  pruned_active.reserve(static_cast<std::size_t>(theta.size()));
  pruned_theta.reserve(static_cast<std::size_t>(theta.size()));
  for (Eigen::Index i = 0; i < theta.size(); ++i) {
    if (theta(i) > prune_threshold) {
      pruned_ranks.push_back(in.active_ranks[static_cast<std::size_t>(i)]);
      pruned_active.push_back(static_cast<std::uint32_t>(i));
      pruned_theta.push_back(theta(i));
    }
  }
  if (pruned_ranks.empty()) return misskappa::unexpected(Error::numerical_error);

  Eigen::VectorXd theta_pruned = Eigen::Map<Eigen::VectorXd>(
      pruned_theta.data(), static_cast<Eigen::Index>(pruned_theta.size()));
  const double s = theta_pruned.sum();
  if (s < zero_tol) return misskappa::unexpected(Error::numerical_error);
  theta_pruned /= s;

  out.theta_hat = std::move(theta_pruned);
  out.pattern_indices = std::move(pruned_ranks);
  out.pruned_active = std::move(pruned_active);
  if (compute_vcov) {
    out.vcov = em_variance(out.theta_hat, out.pruned_active, in, opts);
  }
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
  visit_category_tuples(distance.C, g, visitor);
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
  visit_category_tuples(distance.C, g - 1, visitor);
  return acc;
}

struct GwiseFimlMap {
  Eigen::VectorXd d_vec;
  Eigen::VectorXd grad_C;
  Eigen::VectorXd grad_F;
  double C_hat = 0.0;
  double F_hat = 0.0;
};

Result<GwiseFimlMap> build_gwise_fiml_map(
    const EmRunResult& em, loss::GwiseCategoricalDistance distance,
    GwiseOptions opts) {
  const int R = em.R;
  const int g = (opts.g <= 0) ? R : opts.g;
  if (g < 2 || g > R) return misskappa::unexpected(Error::invalid_argument);
  if (distance.C != em.c || distance.compute == nullptr) {
    return misskappa::unexpected(Error::invalid_argument);
  }

  std::int64_t category_tuples = 0;
  std::int64_t category_projection_tuples = 0;
  if (!checked_power(distance.C, g, opts.max_chance_tuples, category_tuples)
      || !checked_power(distance.C, g - 1, opts.max_chance_tuples, category_projection_tuples)) {
    return misskappa::unexpected(Error::not_supported);
  }

  const auto c_raters = combinations(R, g);
  if (c_raters.empty()) return misskappa::unexpected(Error::invalid_argument);

  const Eigen::Index n_final = em.theta_hat.size();
  const Eigen::VectorXd& theta = em.theta_hat;
  std::vector<std::vector<int>> pattern_matrix;
  pattern_matrix.reserve(static_cast<std::size_t>(n_final));

  RealMat probs = RealMat::Zero(R, distance.C);
  for (Eigen::Index u = 0; u < n_final; ++u) {
    auto row = unrank_tuple(em.pattern_indices[static_cast<std::size_t>(u)], R, em.c);
    for (int r = 0; r < R; ++r) {
      probs(r, row[static_cast<std::size_t>(r)]) += theta(u);
    }
    pattern_matrix.push_back(std::move(row));
  }
  const RealVec pooled = probs.colwise().mean().transpose();

  GwiseFimlMap m;
  m.d_vec = Eigen::VectorXd::Zero(n_final);
  m.grad_C = Eigen::VectorXd::Zero(n_final);
  m.grad_F = Eigen::VectorXd::Zero(n_final);

  std::vector<int> values(static_cast<std::size_t>(g), 0);
  for (Eigen::Index u = 0; u < n_final; ++u) {
    double acc = 0.0;
    for (const auto& raters : c_raters) {
      for (int pos = 0; pos < g; ++pos) {
        values[static_cast<std::size_t>(pos)] =
            pattern_matrix[static_cast<std::size_t>(u)][static_cast<std::size_t>(raters[pos])];
      }
      acc += distance.compute(values.data(), g, distance.C);
    }
    m.d_vec(u) = acc / static_cast<double>(c_raters.size());
  }

  std::vector<RealMat> fixed_by_subset;
  fixed_by_subset.reserve(c_raters.size());
  for (const auto& raters : c_raters) {
    m.C_hat += categorical_expectation(
        distance, g, [&](int pos, int cat) { return probs(raters[pos], cat); });

    RealMat fixed = RealMat::Zero(g, distance.C);
    for (int pos = 0; pos < g; ++pos) {
      for (int cat = 0; cat < distance.C; ++cat) {
        fixed(pos, cat) = categorical_expectation_fixed(
            distance, g, pos, cat,
            [&](int other_pos, int other_cat) { return probs(raters[other_pos], other_cat); });
      }
    }
    fixed_by_subset.push_back(std::move(fixed));
  }
  m.C_hat /= static_cast<double>(c_raters.size());

  for (Eigen::Index u = 0; u < n_final; ++u) {
    double grad = 0.0;
    for (std::size_t s = 0; s < c_raters.size(); ++s) {
      const auto& raters = c_raters[s];
      const RealMat& fixed = fixed_by_subset[s];
      for (int pos = 0; pos < g; ++pos) {
        const int r = raters[static_cast<std::size_t>(pos)];
        const int cat = pattern_matrix[static_cast<std::size_t>(u)][static_cast<std::size_t>(r)];
        grad += fixed(pos, cat);
      }
    }
    m.grad_C(u) = grad / static_cast<double>(c_raters.size());
  }

  m.F_hat = categorical_expectation(
      distance, g, [&](int /*pos*/, int cat) { return pooled(cat); });
  RealMat fleiss_fixed = RealMat::Zero(g, distance.C);
  for (int pos = 0; pos < g; ++pos) {
    for (int cat = 0; cat < distance.C; ++cat) {
      fleiss_fixed(pos, cat) = categorical_expectation_fixed(
          distance, g, pos, cat,
          [&](int /*other_pos*/, int other_cat) { return pooled(other_cat); });
    }
  }

  for (Eigen::Index u = 0; u < n_final; ++u) {
    double grad = 0.0;
    for (int pos = 0; pos < g; ++pos) {
      for (int r = 0; r < R; ++r) {
        const int cat = pattern_matrix[static_cast<std::size_t>(u)][static_cast<std::size_t>(r)];
        grad += fleiss_fixed(pos, cat) / static_cast<double>(R);
      }
    }
    m.grad_F(u) = grad;
  }

  return m;
}

Result<void> validate_alpha_values(const RealVec& values) {
  if (values.size() < 1) return misskappa::unexpected(Error::invalid_argument);
  for (Eigen::Index k = 0; k < values.size(); ++k) {
    if (!std::isfinite(values(k))) return misskappa::unexpected(Error::invalid_argument);
  }
  return {};
}

struct AlphaMap {
  double alpha = std::numeric_limits<double>::quiet_NaN();
  Eigen::VectorXd gradient;
};

AlphaMap build_alpha_map(const EmRunResult& em, const RealVec& values) {
  const int R = em.R;
  const Eigen::Index n_final = em.theta_hat.size();
  const Eigen::VectorXd& theta = em.theta_hat;
  const double factor = static_cast<double>(R) / static_cast<double>(R - 1);

  RealMat pattern_scores = RealMat::Zero(n_final, R);
  RealVec score_sum = RealVec::Zero(n_final);
  RealVec score_sq_sum = RealVec::Zero(n_final);
  RealVec mu = RealVec::Zero(R);
  double e_score_sum_sq = 0.0;
  double e_item_sq = 0.0;
  double mean_score_sum = 0.0;

  for (Eigen::Index u = 0; u < n_final; ++u) {
    const auto row = unrank_tuple(em.pattern_indices[static_cast<std::size_t>(u)], R, em.c);
    double s = 0.0;
    double ss = 0.0;
    for (int j = 0; j < R; ++j) {
      const double y = values(row[static_cast<std::size_t>(j)]);
      pattern_scores(u, j) = y;
      mu(j) += theta(u) * y;
      s += y;
      ss += y * y;
    }
    score_sum(u) = s;
    score_sq_sum(u) = ss;
    mean_score_sum += theta(u) * s;
    e_score_sum_sq += theta(u) * s * s;
    e_item_sq += theta(u) * ss;
  }

  const double t1 = e_score_sum_sq - mean_score_sum * mean_score_sum;
  const double t2 = e_item_sq - mu.squaredNorm();

  AlphaMap out;
  out.gradient = Eigen::VectorXd::Zero(n_final);
  if (std::abs(t1) <= singular_tol) return out;
  out.alpha = factor * (1.0 - t2 / t1);

  for (Eigen::Index u = 0; u < n_final; ++u) {
    const double grad_t1 = score_sum(u) * score_sum(u)
                           - 2.0 * mean_score_sum * score_sum(u);
    double mu_dot_y = 0.0;
    for (int j = 0; j < R; ++j) mu_dot_y += mu(j) * pattern_scores(u, j);
    const double grad_t2 = score_sq_sum(u) - 2.0 * mu_dot_y;
    out.gradient(u) = factor * (t2 * grad_t1 - t1 * grad_t2) / (t1 * t1);
  }
  return out;
}

Result<Estimation> map_fiml_kappa(
    const EmInput& in, const EmRunResult& em, RealMatView weights,
    const FimlVarianceCache& cache) {
  const int n = static_cast<int>(in.subject_groups.size());
  bench_prof::Timer t_rest(bench_prof::acc().rest);
  KappaMap m;
  {
    bench_prof::Timer t(bench_prof::acc().kappa_map);
    m = build_kappa_map(em, weights);
  }

  const Eigen::VectorXd& theta = em.theta_hat;
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

  const RealMat jacobian_reduced = jacobian * cache.theta_jacobian;
  Result<RealVec> null_frac = RealVec{};
  {
    bench_prof::Timer t(bench_prof::acc().null_frac);
    null_frac = reduced_gradient_null_fraction(
        cache.info_star, jacobian_reduced, cache.info_rcond);
  }
  if (!null_frac) return misskappa::unexpected(null_frac.error());

  RealMat vcov = jacobian * cache.theta_vcov * jacobian.transpose();
  RealMat psi = RealMat::Zero(n, 3);
  if (cache.var_star.rows() > 0) {
    const RealMat group_psi =
        static_cast<double>(n) * cache.group_scores * cache.var_star
        * jacobian_reduced.transpose();
    for (Eigen::Index i = 0; i < n; ++i) {
      psi.row(i) = group_psi.row(
          static_cast<Eigen::Index>(in.subject_groups[static_cast<std::size_t>(i)]));
    }
  }

  return Estimation{std::move(estimates), std::move(vcov), std::move(psi),
                    std::move(*null_frac)};
}

}  // namespace

Result<Estimation> estimate_fiml(
    IntMatView ratings, RealMatView weights, EmOptions opts) {
  const int n = static_cast<int>(ratings.rows());
  const int R = static_cast<int>(ratings.cols());
  if (n < 1) return misskappa::unexpected(Error::invalid_argument);
  if (R < 2) return misskappa::unexpected(Error::invalid_argument);
  const int C = static_cast<int>(weights.rows());
  if (weights.cols() != C) return misskappa::unexpected(Error::dimension_mismatch);
  if (C < 1) return misskappa::unexpected(Error::invalid_argument);

  Result<EmInput> in_res = EmInput{};
  {
    bench_prof::Timer t(bench_prof::acc().pre);
    in_res = preprocess_raw(ratings, C);
  }
  if (!in_res) return misskappa::unexpected(in_res.error());
  auto identified = detail::require_complete_pair_observation(detail::observed_mask(ratings));
  if (!identified) return misskappa::unexpected(identified.error());
  auto em = run_em_preprocessed(*in_res, C, opts, false);
  if (!em) return misskappa::unexpected(em.error());
  if (em->theta_hat.size() == 0) return misskappa::unexpected(Error::numerical_error);

  FimlVarianceCache cache = build_fiml_variance_cache(*em, *in_res, opts);
  return map_fiml_kappa(*in_res, *em, weights, cache);
}

Result<std::vector<Estimation>> estimate_fiml_many(
    IntMatView ratings, const std::vector<RealMat>& weights, EmOptions opts) {
  const int n = static_cast<int>(ratings.rows());
  const int R = static_cast<int>(ratings.cols());
  if (n < 1) return misskappa::unexpected(Error::invalid_argument);
  if (R < 2) return misskappa::unexpected(Error::invalid_argument);
  if (weights.empty()) return misskappa::unexpected(Error::invalid_argument);

  const int C = static_cast<int>(weights.front().rows());
  if (C < 1) return misskappa::unexpected(Error::invalid_argument);
  for (const auto& W : weights) {
    if (W.rows() != C || W.cols() != C) {
      return misskappa::unexpected(Error::dimension_mismatch);
    }
  }

  auto in_res = preprocess_raw(ratings, C);
  if (!in_res) return misskappa::unexpected(in_res.error());
  auto identified = detail::require_complete_pair_observation(detail::observed_mask(ratings));
  if (!identified) return misskappa::unexpected(identified.error());
  auto em = run_em_preprocessed(*in_res, C, opts, false);
  if (!em) return misskappa::unexpected(em.error());
  if (em->theta_hat.size() == 0) return misskappa::unexpected(Error::numerical_error);

  FimlVarianceCache cache = build_fiml_variance_cache(*em, *in_res, opts);
  std::vector<Estimation> out;
  out.reserve(weights.size());
  for (const auto& W : weights) {
    auto est = map_fiml_kappa(*in_res, *em, W, cache);
    if (!est) return misskappa::unexpected(est.error());
    out.push_back(std::move(*est));
  }
  return out;
}

Result<Estimation> estimate_fiml_gwise(
    IntMatView ratings, loss::GwiseCategoricalDistance distance,
    EmOptions em_opts, GwiseOptions gwise_opts) {
  const int n = static_cast<int>(ratings.rows());
  const int R = static_cast<int>(ratings.cols());
  const int g = (gwise_opts.g <= 0) ? R : gwise_opts.g;
  if (n < 1) return misskappa::unexpected(Error::invalid_argument);
  if (R < 2 || g < 2 || g > R) return misskappa::unexpected(Error::invalid_argument);
  if (distance.C < 1 || distance.compute == nullptr) {
    return misskappa::unexpected(Error::invalid_argument);
  }

  auto in_res = preprocess_raw(ratings, distance.C);
  if (!in_res) return misskappa::unexpected(in_res.error());
  auto tuple_identified =
      detail::require_complete_tuple_observation(detail::observed_mask(ratings), g);
  if (!tuple_identified) return misskappa::unexpected(tuple_identified.error());
  auto em = run_em_preprocessed(*in_res, distance.C, em_opts, true);
  if (!em) return misskappa::unexpected(em.error());
  if (em->theta_hat.size() == 0) return misskappa::unexpected(Error::numerical_error);

  auto map_res = build_gwise_fiml_map(*em, distance, gwise_opts);
  if (!map_res) return misskappa::unexpected(map_res.error());
  const GwiseFimlMap& m = *map_res;

  const Eigen::VectorXd& theta = em->theta_hat;
  const double pd = m.d_vec.dot(theta);
  const double pec = m.C_hat;
  const double pef = m.F_hat;

  RealVec estimates(2);
  estimates(0) = (std::abs(pec) > singular_tol)
                     ? 1.0 - pd / pec
                     : std::numeric_limits<double>::quiet_NaN();
  estimates(1) = (std::abs(pef) > singular_tol)
                     ? 1.0 - pd / pef
                     : std::numeric_limits<double>::quiet_NaN();

  RealMat jacobian = RealMat::Zero(2, theta.size());
  if (std::abs(pec) > singular_tol) {
    jacobian.row(0) =
        (-(m.d_vec * pec - pd * m.grad_C).transpose()) / (pec * pec);
  }
  if (std::abs(pef) > singular_tol) {
    jacobian.row(1) =
        (-(m.d_vec * pef - pd * m.grad_F).transpose()) / (pef * pef);
  }

  RealMat vcov = jacobian * em->vcov * jacobian.transpose();
  RealMat psi = RealMat::Zero(n, 2);
  RealVec null_frac;
  const LouisReducedInfo info =
      build_louis_reduced_info(em->theta_hat, em->pruned_active, *in_res);
  if (info.info_star.rows() > 0) {
    const RealMat var_star = detail::pseudo_inverse_psd(info.info_star, em_opts.info_rcond);
    const RealMat jacobian_reduced = jacobian * info.jacobian;
    auto frac = reduced_gradient_null_fraction(
        info.info_star, jacobian_reduced, em_opts.info_rcond);
    if (!frac) return misskappa::unexpected(frac.error());
    null_frac = std::move(*frac);
    const RealMat group_psi =
        static_cast<double>(n) * info.group_scores * var_star * jacobian_reduced.transpose();
    for (Eigen::Index i = 0; i < n; ++i) {
      psi.row(i) = group_psi.row(
          static_cast<Eigen::Index>(in_res->subject_groups[static_cast<std::size_t>(i)]));
    }
  }

  return Estimation{std::move(estimates), std::move(vcov), std::move(psi),
                    std::move(null_frac)};
}

Result<Estimation> estimate_alpha_fiml(
    IntMatView ratings, const RealVec& values, EmOptions opts) {
  const int n = static_cast<int>(ratings.rows());
  const int R = static_cast<int>(ratings.cols());
  if (n < 1) return misskappa::unexpected(Error::invalid_argument);
  if (R < 2) return misskappa::unexpected(Error::invalid_argument);
  auto valid = validate_alpha_values(values);
  if (!valid) return misskappa::unexpected(valid.error());
  const int C = static_cast<int>(values.size());

  auto in_res = preprocess_raw(ratings, C);
  if (!in_res) return misskappa::unexpected(in_res.error());
  auto pattern_identified =
      detail::require_complete_pair_observation(detail::observed_mask(ratings));
  if (!pattern_identified) return misskappa::unexpected(pattern_identified.error());
  auto em = run_em_preprocessed(*in_res, C, opts, true);
  if (!em) return misskappa::unexpected(em.error());
  if (em->theta_hat.size() == 0) return misskappa::unexpected(Error::numerical_error);

  const AlphaMap m = build_alpha_map(*em, values);
  RealVec estimates(1);
  estimates(0) = m.alpha;

  RealMat jacobian(1, em->theta_hat.size());
  jacobian.row(0) = m.gradient.transpose();
  RealMat vcov = jacobian * em->vcov * jacobian.transpose();

  RealMat psi = RealMat::Zero(n, 1);
  RealVec null_frac;
  const LouisReducedInfo info =
      build_louis_reduced_info(em->theta_hat, em->pruned_active, *in_res);
  if (info.info_star.rows() > 0) {
    const RealMat var_star = detail::pseudo_inverse_psd(info.info_star, opts.info_rcond);
    const RealMat jacobian_reduced = jacobian * info.jacobian;
    auto frac = reduced_gradient_null_fraction(
        info.info_star, jacobian_reduced, opts.info_rcond);
    if (!frac) return misskappa::unexpected(frac.error());
    null_frac = std::move(*frac);
    const RealMat group_psi =
        static_cast<double>(n) * info.group_scores * var_star * jacobian_reduced.transpose();
    for (Eigen::Index i = 0; i < n; ++i) {
      psi.row(i) = group_psi.row(
          static_cast<Eigen::Index>(in_res->subject_groups[static_cast<std::size_t>(i)]));
    }
  }

  return Estimation{std::move(estimates), std::move(vcov), std::move(psi),
                    std::move(null_frac)};
}

Result<FimlLouisDiagnostic> diagnose_fiml_louis(
    IntMatView ratings, RealMatView weights, EmOptions opts) {
  const int n = static_cast<int>(ratings.rows());
  const int R = static_cast<int>(ratings.cols());
  if (n < 1) return misskappa::unexpected(Error::invalid_argument);
  if (R < 2) return misskappa::unexpected(Error::invalid_argument);
  const int C = static_cast<int>(weights.rows());
  if (weights.cols() != C) return misskappa::unexpected(Error::dimension_mismatch);
  if (C < 1) return misskappa::unexpected(Error::invalid_argument);

  auto in_res = preprocess_raw(ratings, C);
  if (!in_res) return misskappa::unexpected(in_res.error());
  auto em = run_em_preprocessed(*in_res, C, opts, false);
  if (!em) return misskappa::unexpected(em.error());
  if (em->theta_hat.size() == 0) return misskappa::unexpected(Error::numerical_error);

  const KappaMap m = build_kappa_map(*em, weights);
  const Eigen::VectorXd& theta = em->theta_hat;
  const double pd = m.d_vec.dot(theta);
  const double pec = (theta.transpose() * m.Qed_conger * theta).value();
  if (std::abs(pec) <= singular_tol) return misskappa::unexpected(Error::numerical_error);

  const RealMat sym = 0.5 * (m.Qed_conger + m.Qed_conger.transpose());
  const RealVec grad_pec = 2.0 * sym * theta;
  const RealVec grad_full = -(m.d_vec * pec - pd * grad_pec) / (pec * pec);

  const LouisReducedInfo info =
      build_louis_reduced_info(em->theta_hat, em->pruned_active, *in_res);

  FimlLouisDiagnostic out;
  out.c = C;
  out.R = R;
  out.n_subjects = em->n_subjects;
  out.n_patterns = static_cast<std::size_t>(em->theta_hat.size());
  out.kappa_conger = 1.0 - pd / pec;

  const RealVec grad_reduced = info.jacobian.transpose() * grad_full;
  if (info.info_star.rows() == 0) {
    out.eigenvalues = RealVec::Zero(0);
    out.gradient_projection = RealVec::Zero(0);
    out.variance_contribution = RealVec::Zero(0);
    return out;
  }

  Eigen::SelfAdjointEigenSolver<RealMat> es(info.info_star);
  if (es.info() != Eigen::Success) return misskappa::unexpected(Error::numerical_error);
  const RealVec& evals = es.eigenvalues();
  const RealMat& evecs = es.eigenvectors();
  out.lambda_max = evals.maxCoeff();
  const double rc = (std::isfinite(opts.info_rcond) && opts.info_rcond > 0.0)
                        ? opts.info_rcond
                        : 0.0;
  out.threshold = out.lambda_max * rc;
  out.eigenvalues = RealVec::Zero(evals.size());
  out.gradient_projection = RealVec::Zero(evals.size());
  out.variance_contribution = RealVec::Zero(evals.size());

  for (Eigen::Index i = 0; i < evals.size(); ++i) {
    const Eigen::Index src = evals.size() - 1 - i;
    const double lambda = evals(src);
    const double projection = evecs.col(src).dot(grad_reduced);
    const bool retained = lambda > out.threshold;
    const double contribution = retained ? (projection * projection) / lambda : 0.0;
    out.eigenvalues(i) = lambda;
    out.gradient_projection(i) = projection;
    out.variance_contribution(i) = contribution;
    out.variance += contribution;
    if (retained) ++out.retained_rank;
  }
  return out;
}

}  // namespace misskappa
