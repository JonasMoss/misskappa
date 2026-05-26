#include "emdiscrete.h"
#include "misskappa.h"
#include <map>
#include <unordered_map>

// Helper functions for unranking patterns. They are not part of the public API.
namespace {
emdiscrete::uvec unrank_tuple(uint64_t rank, int r, int c) {
  emdiscrete::uvec rt(r, arma::fill::zeros); uint64_t temp_rank = rank;
  for (int i = r - 1; i >= 0; --i) {
    rt(i) = temp_rank % c; temp_rank /= c;
  } return rt;
}
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
  emdiscrete::uvec unrank_composition(uint64_t rank, int r, int c, NChooseK_Memoized& C) {
    emdiscrete::uvec v(c, arma::fill::zeros); int current_r = r; uint64_t temp_rank = rank;
    for (int i = 0; i < c - 1; ++i) {
      if (current_r == 0) break;
      uint64_t combinations = 0; int count = 0;
      while (true) {
        combinations = C.get(current_r - count + c - 1 - i, c - 1 - i);
        if (temp_rank < combinations) break;
        temp_rank -= combinations; count++;
      }
      v(i) = count; current_r -= count;
    }
    if (current_r > 0) v(c - 1) = current_r;
    return v;
  }
}

namespace misskappa {
namespace kappaml {

// Using a C++14-style Estimator type definition
using Estimator = std::function<Result<std::pair<arma::vec, arma::mat>>(const arma::vec&)>;

Result<Estimator> create_estimator_raw(const emdiscrete::EM_Result& em_result, const arma::mat& loss_matrix) {
  const arma::uword n_final = em_result.theta_hat.n_elem;
  const int r = em_result.r; const int c = em_result.c;
  const double n_pairs = static_cast<double>(r) * (r - 1.0) / 2.0;

  arma::umat pattern_matrix(r, n_final);
  for (arma::uword j = 0; j < n_final; ++j) {
    pattern_matrix.col(j) = unrank_tuple(em_result.pattern_indices(j), r, c);
  }

  arma::vec d_vec(n_final, arma::fill::zeros);
  for (arma::uword j = 0; j < n_final; ++j) {
    double disagreement_sum = 0;
    for (int c1 = 0; c1 < r - 1; ++c1) for (int c2 = c1 + 1; c2 < r; ++c2) {
      disagreement_sum += loss_matrix(pattern_matrix(c1, j), pattern_matrix(c2, j));
    }
    d_vec(j) = disagreement_sum / n_pairs;
  }

  arma::mat M_conger(n_final, r * c, arma::fill::zeros);
  for (int rater = 0; rater < r; ++rater) for (int cat = 0; cat < c; ++cat) {
    uvec matching_patterns = arma::find(pattern_matrix.row(rater) == cat);
    if (!matching_patterns.is_empty()) M_conger.submat(matching_patterns, uvec{static_cast<arma::uword>(rater*c+cat)}).ones();
  }
  arma::mat Qed_conger(n_final, n_final, arma::fill::zeros);
  for (int c1 = 0; c1 < r - 1; ++c1) for (int c2 = c1 + 1; c2 < r; ++c2) {
    arma::mat M1 = M_conger.cols(c1 * c, (c1 + 1) * c - 1);
    arma::mat M2 = M_conger.cols(c2 * c, (c2 + 1) * c - 1);
    Qed_conger += M1 * loss_matrix * M2.t();
  }
  Qed_conger /= n_pairs;

  arma::mat M_fleiss(c, n_final, arma::fill::zeros);
  for(arma::uword j=0; j < n_final; ++j) for(int k=0; k < c; ++k) {
    M_fleiss(k, j) = arma::accu(pattern_matrix.col(j) == k);
  }
  M_fleiss /= static_cast<double>(r);
  arma::mat Qed_fleiss = M_fleiss.t() * loss_matrix * M_fleiss;

  const double Ped_bp = arma::accu(loss_matrix) / (c * c);

  Estimator estimator = [=](const arma::vec& theta) -> Result<std::pair<arma::vec, arma::mat>> {
    double Pd = arma::dot(d_vec, theta);
    double Ped_c = arma::as_scalar(theta.t() * Qed_conger * theta);
    double Ped_f = arma::as_scalar(theta.t() * Qed_fleiss * theta);
    arma::vec estimates(3);
    estimates(0) = (std::abs(Ped_c) < 1e-12) ? arma::datum::nan : 1.0 - Pd / Ped_c;
    estimates(1) = (std::abs(Ped_f) < 1e-12) ? arma::datum::nan : 1.0 - Pd / Ped_f;
    estimates(2) = (std::abs(Ped_bp) < 1e-12) ? arma::datum::nan : 1.0 - Pd / Ped_bp;
    arma::mat jacobian(3, theta.n_elem, arma::fill::zeros);
    arma::vec grad_Pd = d_vec;

    if (std::abs(Ped_c) > 1e-12) {
      arma::mat Qed_conger_sym = 0.5 * (Qed_conger + Qed_conger.t());
      arma::vec grad_Ped_c = 2.0 * Qed_conger_sym * theta;
      jacobian.row(0) = - (grad_Pd * Ped_c - Pd * grad_Ped_c).t() / (Ped_c * Ped_c);
    }
    if (std::abs(Ped_f) > 1e-12) {
      arma::mat Qed_fleiss_sym = 0.5 * (Qed_fleiss + Qed_fleiss.t());
      arma::vec grad_Ped_f = 2.0 * Qed_fleiss_sym * theta;
      jacobian.row(1) = - (grad_Pd * Ped_f - Pd * grad_Ped_f).t() / (Ped_f * Ped_f);
    }
    if (std::abs(Ped_bp) > 1e-12) {
      jacobian.row(2) = - grad_Pd.t() / Ped_bp;
    }
    return {emdiscrete::Status::kOk, std::make_pair(estimates, jacobian), ""};
  };
  return {emdiscrete::Status::kOk, estimator, ""};
}

Result<Estimator> create_estimator_counts(const emdiscrete::EM_Result& em_result, const arma::mat& loss_matrix) {
  const arma::uword n_final = em_result.theta_hat.n_elem;
  const int r = em_result.r; const int c = em_result.c;
  double r_pair = static_cast<double>(r) * (r - 1.0);

  NChooseK_Memoized C(r + c);
  arma::mat pattern_matrix(c, n_final);
  for (arma::uword j = 0; j < n_final; ++j) {
    pattern_matrix.col(j) = arma::conv_to<arma::vec>::from(
      unrank_composition(em_result.pattern_indices(j), r, c, C));
  }

  arma::mat agreement_matrix = arma::ones(c, c) - loss_matrix;
  arma::vec pa_vec(n_final);
  for (arma::uword i = 0; i < n_final; ++i) {
    arma::vec k = pattern_matrix.col(i);
    pa_vec(i) = (arma::as_scalar(k.t() * agreement_matrix * k) - arma::dot(agreement_matrix.diag(), k)) / r_pair;
  }
  arma::mat A_map = (1.0 / r) * pattern_matrix.t();
  const double Pe_bp = arma::accu(agreement_matrix) / (c * c);

  Estimator estimator = [=](const arma::vec& theta) -> Result<std::pair<arma::vec, arma::mat>> {
    double Pa = arma::dot(pa_vec, theta);
    arma::vec p_hat = A_map.t() * theta;
    double Pe_f = arma::as_scalar(p_hat.t() * agreement_matrix * p_hat);

    arma::vec estimates(2);
    estimates(0) = (std::abs(1.0 - Pe_f) < 1e-12) ? arma::datum::nan : (Pa - Pe_f) / (1.0 - Pe_f);
    estimates(1) = (std::abs(1.0 - Pe_bp) < 1e-12) ? arma::datum::nan : (Pa - Pe_bp) / (1.0 - Pe_bp);
    arma::mat jacobian(2, theta.n_elem, arma::fill::zeros);
    arma::vec grad_Pa = pa_vec;

    if (std::abs(1.0 - Pe_f) > 1e-12) {
      arma::vec grad_Pe_f = 2.0 * A_map * agreement_matrix * p_hat;
      jacobian.row(0) = ((grad_Pa - grad_Pe_f)*(1.0-Pe_f) - (Pa-Pe_f)*(-grad_Pe_f)).t() / std::pow(1.0-Pe_f, 2);
    }
    if (std::abs(1.0 - Pe_bp) > 1e-12) {
      jacobian.row(1) = grad_Pa.t() / (1.0 - Pe_bp);
    }
    return {emdiscrete::Status::kOk, std::make_pair(estimates, jacobian), ""};
  };
  return {emdiscrete::Status::kOk, estimator, ""};
}

Result<Estimation> kappa(const emdiscrete::EM_Result& em_res, const arma::mat& loss_matrix) {
  if (em_res.theta_hat.is_empty()) return {emdiscrete::Status::kError, std::nullopt, "EM result has empty theta."};
  auto factory = create_estimator_raw(em_res, loss_matrix);
  if (!factory.IsOk()) return {emdiscrete::Status::kError, std::nullopt, factory.error_message};

  auto est_res = factory.value.value()(em_res.theta_hat);
  if (!est_res.IsOk()) return {emdiscrete::Status::kError, std::nullopt, est_res.error_message};

  arma::mat vcov = est_res.value.value().second * em_res.var * est_res.value.value().second.t();
  return {emdiscrete::Status::kOk, Estimation{est_res.value.value().first, vcov}, ""};
}

Result<Estimation> kappa_counts(const emdiscrete::EM_Result& em_res, const arma::mat& loss_matrix) {
  if (em_res.theta_hat.is_empty()) return {emdiscrete::Status::kError, std::nullopt, "EM result has empty theta."};
  auto factory = create_estimator_counts(em_res, loss_matrix);
  if (!factory.IsOk()) return {emdiscrete::Status::kError, std::nullopt, factory.error_message};

  auto est_res = factory.value.value()(em_res.theta_hat);
  if (!est_res.IsOk()) return {emdiscrete::Status::kError, std::nullopt, est_res.error_message};

  arma::mat vcov = est_res.value.value().second * em_res.var * est_res.value.value().second.t();
  return {emdiscrete::Status::kOk, Estimation{est_res.value.value().first, vcov}, ""};
}

} // namespace kappaml
} // namespace misskappa
