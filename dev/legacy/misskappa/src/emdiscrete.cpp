#include "emdiscrete.h"
#include <map>
#include <unordered_map>
#include <numeric>

namespace emdiscrete {

// Anonymous namespace for implementation details internal to this file
namespace {
namespace preprocessing_internals {
// --- Helpers for "raw" (tuple) model ---
inline uint64_t rank_tuple(const uvec& rt, int c) {
  uint64_t rank = 0; uint64_t power = 1;
  for (int i = rt.n_elem - 1; i >= 0; --i) {
    rank += rt(i) * power; if (i > 0) power *= c;
  } return rank;
}
inline void generate_tuple_completions_recursive(
    uvec& p, const uvec& mi, arma::uword k, int c, std::vector<uvec>& res) {
  if (k == mi.n_elem) { res.push_back(p); return; }
  arma::uword rater_idx_to_fill = mi(k);
  for (int cat = 0; cat < c; ++cat) {
    p(rater_idx_to_fill) = cat;
    generate_tuple_completions_recursive(p, mi, k + 1, c, res);
  }
}
// --- Helpers for "counts" (composition) model ---
struct NChooseK_Memoized {
  std::vector<std::vector<uint64_t>> memo;
  NChooseK_Memoized(int max_n) { memo.resize(max_n + 1, std::vector<uint64_t>(max_n + 1, 0)); }
  uint64_t get(int n, int k) {
    if (k < 0 || k > n) return 0;
    if (k == 0 || k == n) return 1;
    if (k > n / 2) k = n - k;
    if (memo[n][k] != 0) return memo[n][k];
    uint64_t res = get(n - 1, k - 1) + get(n - 1, k);
    memo[n][k] = res; return res;
  }
};
inline uint64_t rank_composition(const uvec& v, int r, NChooseK_Memoized& C) {
  arma::uword c = v.n_elem; uint64_t rank = 0; int current_r = r;
  for (arma::uword i = 0; i < c - 1; ++i) {
    if (v(i) > 0) for (arma::uword j = 0; j < v(i); ++j) rank += C.get(current_r - j + (c - 1 - i), (c - 1 - i));
    current_r -= v(i);
  } return rank;
}
inline void generate_compositions_recursive(int n, int k, uvec& curr, std::vector<uvec>& res) {
  if (k == 1) { curr(curr.n_elem - k) = n; res.push_back(curr); return; }
  for (int i = 0; i <= n; ++i) {
    curr(curr.n_elem - k) = i;
    generate_compositions_recursive(n - i, k - 1, curr, res);
  }
}
inline std::vector<uvec> generate_compositions(int n, int k) {
  std::vector<uvec> res; if (k <= 0) return res;
  uvec current(k, arma::fill::zeros);
  generate_compositions_recursive(n, k, current, res); return res;
}
} // namespace preprocessing_internals

namespace em_internals {
inline Result<arma::vec> initialize_theta(const EM_Input& em_input, double start_alpha) {
  arma::vec theta(em_input.n_total_patterns, arma::fill::value(start_alpha));
  if (em_input.n_total_patterns == 0) return {Status::kOk, theta, ""};
  for (arma::uword i = 0; i < em_input.group_n_subjects.n_elem; ++i) {
    arma::uword n_subs = em_input.group_n_subjects(i);
    arma::uword n_comps = em_input.group_n_completions(i);
    if (n_comps == 0) continue;
    arma::uword offset = em_input.group_offsets(i);
    const uvec& comp_indices = em_input.completion_indices.subvec(offset, offset + n_comps - 1);
    double weight = static_cast<double>(n_subs) / n_comps;
    theta.elem(comp_indices) += weight;
  }
  double sum_theta = arma::sum(theta);
  if (sum_theta < 1e-9) return {Status::kError, std::nullopt, "Failed to initialize theta."};
  theta /= sum_theta;
  return {Status::kOk, theta, ""};
}

inline Result<bool> run_em_iterations(
    arma::vec& theta, int& iterations, bool& converged,
    const EM_Input& em_input, const EM_Options& options) {
  arma::uword n_total_subjects = arma::sum(em_input.group_n_subjects);
  if (n_total_subjects == 0) {
    iterations = 0; converged = true; return {Status::kOk, true, ""};
  }
  for (int it = 0; it < options.max_iter; ++it) {
    iterations = it + 1;
    arma::vec expected_counts(em_input.n_total_patterns, arma::fill::zeros);
    arma::vec new_theta;
    if (em_input.type == "counts") {
      for (arma::uword i = 0; i < em_input.group_n_subjects.n_elem; ++i) {
        arma::uword offset = em_input.group_offsets(i);
        arma::uword n_comps = em_input.group_n_completions(i);
        if (n_comps == 0) continue;
        const uvec& comp_indices = em_input.completion_indices.subvec(offset, offset + n_comps - 1);
        arma::vec theta_subset = theta.elem(comp_indices);
        const arma::vec& mvhg_coeffs_slice = em_input.multivariate_hypergeom_coeffs.subvec(offset, offset + n_comps - 1);
        arma::vec weights = mvhg_coeffs_slice % theta_subset;
        double sum_w = arma::sum(weights);
        if (sum_w > 1e-12) {
          expected_counts.elem(comp_indices) += em_input.group_n_subjects(i) * (weights / sum_w);
        }
      }
      arma::vec new_theta_unnormalized = expected_counts;
      double total_unnorm = arma::sum(new_theta_unnormalized);
      if (total_unnorm < 1e-12) {
        new_theta.set_size(theta.n_elem);
        new_theta.fill(1.0 / theta.n_elem);
      } else {
        new_theta = new_theta_unnormalized / total_unnorm;
      }
    } else { // "raw" model
      for (arma::uword i = 0; i < em_input.group_n_subjects.n_elem; ++i) {
        arma::uword offset = em_input.group_offsets(i);
        arma::uword n_comps = em_input.group_n_completions(i);
        if (n_comps == 0) continue;
        const uvec& comp_indices = em_input.completion_indices.subvec(offset, offset + n_comps - 1);
        arma::vec theta_subset = theta.elem(comp_indices);
        double sum_w = arma::sum(theta_subset);
        if (sum_w > 1e-12) {
          expected_counts.elem(comp_indices) += em_input.group_n_subjects(i) * (theta_subset / sum_w);
        }
      }
      new_theta = expected_counts / static_cast<double>(n_total_subjects);
    }
    if (arma::abs(new_theta - theta).max() < options.tol) {
      converged = true; theta = new_theta; return {Status::kOk, true, ""};
    }
    theta = new_theta;
  }
  converged = false; return {Status::kOk, true, ""};
}

inline Result<arma::mat> calculate_em_variance(
    const arma::vec& theta_hat, const uvec& pattern_indices, const EM_Input& em_input) {
  const arma::uword n_final = theta_hat.n_elem;
  if (n_final <= 1) return {Status::kOk, arma::mat(n_final, n_final, arma::fill::zeros), ""};
  arma::uword ref_idx_final = theta_hat.index_max();
  arma::mat info_star(n_final - 1, n_final - 1, arma::fill::zeros);
  std::unordered_map<arma::uword, arma::uword> original_to_pruned_map;
  for (arma::uword i = 0; i < n_final; ++i) original_to_pruned_map[pattern_indices(i)] = i;
  for (arma::uword i = 0; i < em_input.group_n_subjects.n_elem; ++i) {
    arma::uword offset = em_input.group_offsets(i);
    arma::uword n_comps = em_input.group_n_completions(i);
    std::vector<arma::uword> final_comps_indices_in_group;
    std::vector<arma::uword> final_comps_indices_in_theta;
    for(arma::uword j=0; j<n_comps; ++j){
      arma::uword original_idx = em_input.completion_indices(offset + j);
      if(original_to_pruned_map.count(original_idx)){
        final_comps_indices_in_group.push_back(j);
        final_comps_indices_in_theta.push_back(original_to_pruned_map.at(original_idx));
      }
    }
    if (final_comps_indices_in_group.empty()) continue;
    uvec group_final_indices = arma::conv_to<uvec>::from(final_comps_indices_in_theta);
    arma::vec theta_subset = theta_hat.elem(group_final_indices);
    arma::vec posterior_probs;
    if (em_input.type == "counts") {
      uvec absolute_indices_to_select = offset + uvec(final_comps_indices_in_group);
      arma::vec mvhg_coeffs_slice = em_input.multivariate_hypergeom_coeffs.elem(absolute_indices_to_select);
      arma::vec weights = mvhg_coeffs_slice % theta_subset;
      double sum_weights = arma::sum(weights);
      if (sum_weights > 1e-12) posterior_probs = weights / sum_weights;
      else continue;
    } else { // "raw" model
      double sum_theta_subset = arma::sum(theta_subset);
      if(sum_theta_subset > 1e-12) posterior_probs = theta_subset / sum_theta_subset;
      else continue;
    }
    if (posterior_probs.has_nan()) continue;
    arma::vec s_star_final(n_final, arma::fill::zeros);
    s_star_final.elem(group_final_indices) = posterior_probs / theta_subset;
    arma::vec s_star_reduced(n_final - 1);
    arma::uword counter = 0;
    for (arma::uword k = 0; k < n_final; ++k) {
      if (k != ref_idx_final) {
        s_star_reduced(counter++) = s_star_final(k) - s_star_final(ref_idx_final);
      }
    }
    info_star += em_input.group_n_subjects(i) * (s_star_reduced * s_star_reduced.t());
  }
  info_star = 0.5 * (info_star + info_star.t());
  arma::mat var_star = arma::pinv(info_star);
  arma::mat J(n_final, n_final - 1, arma::fill::zeros);
  arma::uword counter = 0;
  for (arma::uword i = 0; i < n_final; ++i) {
    if (i != ref_idx_final) J(i, counter++) = 1.0;
  }
  J.row(ref_idx_final).fill(-1.0);
  return {Status::kOk, J * var_star * J.t(), ""};
}
} // namespace em_internals
} // anonymous namespace


// --- Function Implementations ---
Result<EM_Input> preprocess_raw(const arma::imat& x, int c) {
  if (x.n_rows == 0 || x.n_cols == 0) return {Status::kError, std::nullopt, "Input data empty."};
  if (c <= 0) return {Status::kError, std::nullopt, "'c' must be positive."};
  EM_Input res; res.r = x.n_cols; res.c = c; res.type = "raw";
  arma::imat x_zb(x.n_rows, x.n_cols);
  for (arma::uword i = 0; i < x.n_elem; ++i) {
    const int val = x(i);
    if (val == kNaInteger) x_zb(i) = -1;
    else {
      int zb_cat = val;
      if (zb_cat < 0 || zb_cat >= c) return {Status::kError, std::nullopt, "Category " + std::to_string(val) + " out of range."};
      x_zb(i) = zb_cat;
    }
  }
  std::map<std::string, int> pc; std::map<std::string, arma::ivec> up; std::vector<std::string> uo;
  for (arma::uword i = 0; i < x_zb.n_rows; ++i) {
    std::string key = ""; arma::ivec row = x_zb.row(i).t();
    for (int val : row) key += std::to_string(val) + ",";
    if (pc.find(key) == pc.end()) { up[key] = row; uo.push_back(key); }
    pc[key]++;
  }
  std::vector<arma::uword> acv, gov, gncv, gnsv; arma::uword co = 0;
  for (const auto& key : uo) {
    arma::ivec ob_p = up.at(key); uvec mi = arma::find(ob_p == -1); std::vector<uvec> comps;
    uvec p_u = arma::conv_to<uvec>::from(ob_p);
    if (mi.n_elem > 0) preprocessing_internals::generate_tuple_completions_recursive(p_u, mi, 0, res.c, comps);
    else comps.push_back(p_u);
    gnsv.push_back(pc.at(key)); gov.push_back(co); gncv.push_back(comps.size());
    for (const auto& ct : comps) acv.push_back(preprocessing_internals::rank_tuple(ct, res.c));
    co += comps.size();
  }
  res.n_total_patterns = static_cast<size_t>(round(pow(res.c, res.r)));
  res.group_n_subjects = uvec(gnsv); res.group_offsets = uvec(gov);
  res.group_n_completions = uvec(gncv); res.completion_indices = uvec(acv);
  return {Status::kOk, res, ""};
}

Result<EM_Input> preprocess_counts(const arma::umat& x, int r) {
  if (x.n_rows == 0 || x.n_cols == 0) return {Status::kError, std::nullopt, "Input data empty."};
  const int c = x.n_cols;

  preprocessing_internals::NChooseK_Memoized C(r + c);

  std::map<std::string, arma::uword> pattern_counts; std::vector<uvec> unique_patterns;
  for (arma::uword i = 0; i < x.n_rows; ++i) {
    std::string key = ""; uvec row = x.row(i).t();
    for (arma::uword val : row) key += std::to_string(val) + ",";
    if (pattern_counts.find(key) == pattern_counts.end()) unique_patterns.push_back(row);
    pattern_counts[key]++;
  }

  std::unordered_map<int, std::vector<uvec>> comp_cache;
  std::map<uint64_t, bool> active_ranks_set;
  std::vector<std::vector<uint64_t>> completions_by_group;
  std::vector<std::vector<double>> mvhg_coeffs_by_group;

  for (const auto& p : unique_patterns) {
    int r_obs = arma::sum(p); int s = r - r_obs;
    if (s < 0) return {Status::kError, std::nullopt, "Subject has more ratings than 'r'."};
    if (comp_cache.find(s) == comp_cache.end()) comp_cache[s] = preprocessing_internals::generate_compositions(s, c);

    const auto& exts = comp_cache.at(s);
    std::vector<uint64_t> group_comps;
    std::vector<double> group_coeffs;
    for (const auto& ext : exts) {
      uvec z_comp = p + ext;
      uint64_t rank = preprocessing_internals::rank_composition(z_comp, r, C);
      active_ranks_set[rank] = true;
      group_comps.push_back(rank);

      double product_coeff = 1.0;
      for (arma::uword k=0; k < (unsigned int)c; ++k) {
        product_coeff *= C.get(z_comp(k), p(k));
      }
      group_coeffs.push_back(product_coeff);
    }
    completions_by_group.push_back(group_comps);
    mvhg_coeffs_by_group.push_back(group_coeffs);
  }

  std::vector<arma::uword> sorted_ranks; for (const auto& pair : active_ranks_set) sorted_ranks.push_back(pair.first);
  std::unordered_map<uint64_t, arma::uword> rank_to_idx;
  for (arma::uword i = 0; i < sorted_ranks.size(); ++i) rank_to_idx[sorted_ranks[i]] = i;

  EM_Input res_active; res_active.r = r; res_active.c = c; res_active.type = "counts";
  res_active.n_total_patterns = sorted_ranks.size();
  res_active.pattern_indices_map = uvec(sorted_ranks);

  std::vector<arma::uword> acv, gov, gncv, gnsv;
  std::vector<double> accv_mvhg;
  arma::uword co = 0;
  for (size_t i = 0; i < unique_patterns.size(); ++i) {
    std::string key = ""; for(arma::uword val : unique_patterns[i]) key += std::to_string(val) + ",";
    gnsv.push_back(pattern_counts.at(key));
    gov.push_back(co);

    const auto& comps = completions_by_group[i];
    const auto& coeffs = mvhg_coeffs_by_group[i];

    gncv.push_back(comps.size());
    for (size_t k = 0; k < comps.size(); ++k) {
      acv.push_back(rank_to_idx.at(comps[k]));
      accv_mvhg.push_back(coeffs[k]);
    }
    co += comps.size();
  }
  res_active.group_n_subjects = uvec(gnsv); res_active.group_offsets = uvec(gov);
  res_active.group_n_completions = uvec(gncv);
  res_active.completion_indices = uvec(acv);
  res_active.multivariate_hypergeom_coeffs = arma::vec(accv_mvhg);

  return {Status::kOk, res_active, ""};
}

Result<EM_Result> run_em(const EM_Input& em_input, const EM_Options& options) {
  EM_Result result;
  result.type = em_input.type; result.c = em_input.c; result.r = em_input.r;
  result.n_subjects = arma::sum(em_input.group_n_subjects);

  if (em_input.n_total_patterns == 0) {
    result.converged = true; result.iterations = 0;
    return {Status::kOk, result, ""};
  }
  auto theta_res = em_internals::initialize_theta(em_input, options.start_alpha);
  if (!theta_res.IsOk()) return {Status::kError, std::nullopt, theta_res.error_message};
  arma::vec theta = theta_res.value.value();
  result.converged = false;
  auto iter_res = em_internals::run_em_iterations(theta, result.iterations, result.converged, em_input, options);
  if (!iter_res.IsOk()) return {Status::kError, std::nullopt, iter_res.error_message};

  uvec pruned_indices = arma::find(theta > options.prune_tol);
  if (pruned_indices.is_empty()) {
    result.theta_hat.clear();
    result.pattern_indices.clear();
    result.var.clear();
    return {Status::kOk, result, ""};
  }

  if(em_input.type == "counts") {
    result.pattern_indices = em_input.pattern_indices_map.elem(pruned_indices);
  } else {
    result.pattern_indices = pruned_indices;
  }

  result.theta_hat = theta.elem(pruned_indices);
  if (arma::sum(result.theta_hat) > 1e-9) result.theta_hat /= arma::sum(result.theta_hat);

  uvec variance_indices = arma::find(theta > options.prune_tol);
  auto var_res = em_internals::calculate_em_variance(result.theta_hat, variance_indices, em_input);
  if (!var_res.IsOk()) return {Status::kError, std::nullopt, var_res.error_message};
  result.var = var_res.value.value();

  return {Status::kOk, result, ""};
}

} // namespace emdiscrete
