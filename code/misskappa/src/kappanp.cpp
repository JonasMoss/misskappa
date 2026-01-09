#include "misskappa.h"
#include <set>
#include <map>
#include <algorithm> // For std::lexicographical_compare
#include <tuple>     // For std::tuple

namespace misskappa {
namespace kappanp {

struct UrowvecComparator {
  bool operator()(const arma::urowvec& a, const arma::urowvec& b) const {
    return std::lexicographical_compare(a.begin(), a.end(), b.begin(), b.end());
  }
};

struct IvecComparator {
  bool operator()(const arma::ivec& a, const arma::ivec& b) const {
    return std::lexicographical_compare(a.begin(), a.end(), b.begin(), b.end());
  }
};

Result<std::tuple<arma::vec, arma::mat>> calculate_inverse_weights(
    const arma::umat& M,
    int n,
    int R,
    bool use_ipw,
    bool use_gwet
) {
  arma::vec pi_j_inv(R);
  arma::mat pi_jk_inv(R, R);

  if (use_gwet) {
    // Gwet-style weighting (overrides IPW).
    arma::vec pi_j = arma::conv_to<arma::vec>::from(arma::sum(M, 0)) / n;
    if (arma::any(pi_j < 1e-9)) {
      return {emdiscrete::Status::kError, std::nullopt, "One or more raters has no ratings. Cannot use Gwet weighting."};
    }
    pi_j_inv = 1.0 / pi_j;
    pi_jk_inv.ones(); // Pairwise weights are not used for observed disagreement.
  } else if (use_ipw) {
    // Full IPW weighting.
    arma::vec pi_j = arma::conv_to<arma::vec>::from(arma::sum(M, 0)) / n;
    if (arma::any(pi_j < 1e-9)) {
      return {emdiscrete::Status::kError, std::nullopt, "One or more raters has no ratings. Cannot use IPW."};
    }
    pi_j_inv = 1.0 / pi_j;

    pi_jk_inv = arma::conv_to<arma::mat>::from(M.t() * M) / n;
    pi_jk_inv.for_each([](arma::mat::elem_type& val) {
      val = (val > 1e-9) ? 1.0 / val : 0.0;
    });
  } else {
    // No weighting: all weights are 1.0.
    pi_j_inv.ones();
    pi_jk_inv.ones();
  }

  return {emdiscrete::Status::kOk, std::make_tuple(std::move(pi_j_inv), std::move(pi_jk_inv)), ""};
}

Result<Estimation> kappa(
    const arma::imat& ratings,
    const arma::mat& loss_matrix,
    bool use_ipw,
    bool use_gwet // overrides use_ipw.
) {
  const int n = ratings.n_rows;
  const int R = ratings.n_cols;
  if (R < 2) return {emdiscrete::Status::kError, std::nullopt, "Input 'ratings' must have at least 2 raters (columns)."};

  std::map<arma::ivec, int, IvecComparator> row_map;
  std::vector<arma::ivec> uniq_rows_vec;
  arma::uvec counts_vec;
  arma::uvec inv_idx_vec(n);

  for (int i = 0; i < n; ++i) {
    arma::ivec key = ratings.row(i).t();
    auto it = row_map.find(key);
    if (it == row_map.end()) {
      int id = uniq_rows_vec.size();
      row_map[key] = id;
      uniq_rows_vec.push_back(key);
      counts_vec.resize(id + 1);
      counts_vec(id) = 1;
      inv_idx_vec(i) = id;
    } else {
      inv_idx_vec(i) = it->second;
      counts_vec(it->second)++;
    }
  }
  const int U = uniq_rows_vec.size();
  if (U == 0) return {emdiscrete::Status::kError, std::nullopt, "Input 'ratings' is empty or contains no valid patterns."};

  arma::imat xU(U, R);
  for(int u = 0; u < U; ++u) xU.row(u) = uniq_rows_vec[u].t();

  std::set<int> cat_set;
  for (arma::uword i = 0; i < ratings.n_elem; ++i) {
    if (ratings(i) != emdiscrete::kNaInteger) cat_set.insert(ratings(i));
  }
  const int C = cat_set.size();
  if (C == 0) return {emdiscrete::Status::kError, std::nullopt, "All ratings are missing."};
  if (static_cast<int>(loss_matrix.n_rows) != C) return {emdiscrete::Status::kError, std::nullopt, "Loss matrix dimensions do not match number of categories."};

  std::map<int, int> cat_to_idx;
  int current_idx = 0;
  for (int cat : cat_set) cat_to_idx[cat] = current_idx++;

  const arma::umat M = (ratings != emdiscrete::kNaInteger);
  const arma::umat MU = (xU != emdiscrete::kNaInteger);

  auto weights_result = calculate_inverse_weights(M, n, R, use_ipw, use_gwet);
  if (weights_result.status != emdiscrete::Status::kOk) {
    return {weights_result.status, std::nullopt, weights_result.error_message};
  }
  const auto& [pi_j_inv, pi_jk_inv] = *weights_result.value;

  arma::vec h_dN_u(U, arma::fill::zeros), h_dD_u(U, arma::fill::zeros);
  for (int u = 0; u < U; ++u) {
    for (int j = 0; j < R - 1; ++j) {
      for (int k = j + 1; k < R; ++k) {
        if (MU(u, j) && MU(u, k)) {
          double w = pi_jk_inv(j, k);
          int cat1_idx = cat_to_idx.at(xU(u, j));
          int cat2_idx = cat_to_idx.at(xU(u, k));
          h_dN_u(u) += loss_matrix(cat1_idx, cat2_idx) * w;
          h_dD_u(u) += w;
        }
      }
    }
  }
  const double psi_dN_sum = arma::dot(counts_vec, h_dN_u);
  const double psi_dD_sum = arma::dot(counts_vec, h_dD_u);

  arma::mat kernel_CN_uv(U, U), kernel_CD_uv(U, U), kernel_FN_uv(U, U), kernel_FD_uv(U, U);
  for (int u=0; u<U; ++u) for (int v=0; v<U; ++v) {
    double h_cn=0, h_cd=0, h_fn=0, h_fd=0;
    for (int j=0; j<R-1; ++j) for (int k=j+1; k<R; ++k) if(MU(u,j) && MU(v,k)) {
      double w = pi_j_inv(j) * pi_j_inv(k);
      h_cn += loss_matrix(cat_to_idx.at(xU(u,j)), cat_to_idx.at(xU(v,k))) * w;
      h_cd += w;
    }
    for (int j=0; j<R; ++j) for (int k=0; k<R; ++k) if(MU(u,j) && MU(v,k)) {
      double w = pi_j_inv(j) * pi_j_inv(k);
      h_fn += loss_matrix(cat_to_idx.at(xU(u,j)), cat_to_idx.at(xU(v,k))) * w;
      h_fd += w;
    }
    kernel_CN_uv(u,v)=h_cn; kernel_CD_uv(u,v)=h_cd;
    kernel_FN_uv(u,v)=h_fn; kernel_FD_uv(u,v)=h_fd;
  }
  arma::vec counts_d = arma::conv_to<arma::vec>::from(counts_vec);
  const double psi_CN_sum = arma::as_scalar(counts_d.t() * kernel_CN_uv * counts_d);
  const double psi_CD_sum = arma::as_scalar(counts_d.t() * kernel_CD_uv * counts_d);
  const double psi_FN_sum = arma::as_scalar(counts_d.t() * kernel_FN_uv * counts_d);
  const double psi_FD_sum = arma::as_scalar(counts_d.t() * kernel_FD_uv * counts_d);

  const double d_hat = (psi_dD_sum > 1e-9) ? psi_dN_sum / psi_dD_sum : 0.0;
  const double d_C_hat = (psi_CD_sum > 1e-9) ? psi_CN_sum / psi_CD_sum : 0.0;
  const double d_F_hat = (psi_FD_sum > 1e-9) ? psi_FN_sum / psi_FD_sum : 0.0;
  const double d_BP = arma::accu(loss_matrix) / (C * C);

  const double psi_dN_hat = psi_dN_sum / n;
  const double psi_dD_hat = psi_dD_sum / n;
  const double psi_CN_hat = psi_CN_sum / (n * n);
  const double psi_CD_hat = psi_CD_sum / (n * n);
  const double psi_FN_hat = psi_FN_sum / (n * n);
  const double psi_FD_hat = psi_FD_sum / (n * n);

  arma::vec phi_dN = h_dN_u.elem(inv_idx_vec) - psi_dN_hat;
  arma::vec phi_dD = h_dD_u.elem(inv_idx_vec) - psi_dD_hat;

  arma::vec h_CN1_u = (kernel_CN_uv * counts_d) / n;
  arma::vec h_CN2_u = (kernel_CN_uv.t() * counts_d) / n;
  arma::vec phi_CN = (h_CN1_u.elem(inv_idx_vec) - psi_CN_hat) + (h_CN2_u.elem(inv_idx_vec) - psi_CN_hat);

  arma::vec h_CD1_u = (kernel_CD_uv * counts_d) / n;
  arma::vec h_CD2_u = (kernel_CD_uv.t() * counts_d) / n;
  arma::vec phi_CD = (h_CD1_u.elem(inv_idx_vec) - psi_CD_hat) + (h_CD2_u.elem(inv_idx_vec) - psi_CD_hat);

  arma::vec h_FN1_u = (kernel_FN_uv * counts_d) / n;
  arma::vec h_FN2_u = (kernel_FN_uv.t() * counts_d) / n;
  arma::vec phi_FN = (h_FN1_u.elem(inv_idx_vec) - psi_FN_hat) + (h_FN2_u.elem(inv_idx_vec) - psi_FN_hat);

  arma::vec h_FD1_u = (kernel_FD_uv * counts_d) / n;
  arma::vec h_FD2_u = (kernel_FD_uv.t() * counts_d) / n;
  arma::vec phi_FD = (h_FD1_u.elem(inv_idx_vec) - psi_FD_hat) + (h_FD2_u.elem(inv_idx_vec) - psi_FD_hat);

  arma::mat phi_matrix(n, 6);
  phi_matrix.col(0) = phi_dN; phi_matrix.col(1) = phi_dD;
  phi_matrix.col(2) = phi_CN; phi_matrix.col(3) = phi_CD;
  phi_matrix.col(4) = phi_FN; phi_matrix.col(5) = phi_FD;

  arma::mat Gamma_hat = (1.0 / n) * phi_matrix.t() * phi_matrix;

  arma::mat J_d(3, 6, arma::fill::zeros);
  if (psi_dD_hat > 1e-9) { J_d(0, 0)=1/psi_dD_hat; J_d(0, 1)=-psi_dN_hat/std::pow(psi_dD_hat,2); }
  if (psi_CD_hat > 1e-9) { J_d(1, 2)=1/psi_CD_hat; J_d(1, 3)=-psi_CN_hat/std::pow(psi_CD_hat,2); }
  if (psi_FD_hat > 1e-9) { J_d(2, 4)=1/psi_FD_hat; J_d(2, 5)=-psi_FN_hat/std::pow(psi_FD_hat,2); }
  arma::mat disagreement_cov_matrix = (J_d * Gamma_hat * J_d.t()) / n;

  arma::mat J_kappa(3, 3, arma::fill::zeros);
  if (d_C_hat > 1e-9) { J_kappa(0, 0)=-1/d_C_hat; J_kappa(0, 1)=d_hat/std::pow(d_C_hat,2); }
  if (d_F_hat > 1e-9) { J_kappa(1, 0)=-1/d_F_hat; J_kappa(1, 2)=d_hat/std::pow(d_F_hat,2); }
  if (d_BP > 1e-9)    { J_kappa(2, 0)=-1/d_BP; }

  arma::mat kappa_cov_matrix = J_kappa * disagreement_cov_matrix * J_kappa.t();

  arma::vec estimates(3);
  estimates(0) = (d_C_hat > 1e-9) ? 1.0 - d_hat / d_C_hat : arma::datum::nan;
  estimates(1) = (d_F_hat > 1e-9) ? 1.0 - d_hat / d_F_hat : arma::datum::nan;
  estimates(2) = (d_BP > 1e-9) ? 1.0 - d_hat / d_BP : arma::datum::nan;

  return {emdiscrete::Status::kOk, Estimation{estimates, kappa_cov_matrix}, ""};
}

Result<Estimation> kappa_counts(
    const arma::umat& counts,
    const arma::mat& loss_matrix
) {
  const int n = counts.n_rows;
  const int C = counts.n_cols;
  if (n == 0) return {emdiscrete::Status::kError, std::nullopt, "Input counts data is empty."};
  if (static_cast<int>(loss_matrix.n_rows) != C) return {emdiscrete::Status::kError, std::nullopt, "Loss matrix dimensions do not match number of categories."};

  // A. Setup & Binning
  std::map<arma::urowvec, int, UrowvecComparator> row_map;
  arma::umat xU_counts;
  arma::uvec counts_of_counts_vec;
  arma::uvec inv_idx_vec(n);

  for (int i = 0; i < n; ++i) {
    arma::urowvec key = counts.row(i);
    auto it = row_map.find(key);
    if (it == row_map.end()) {
      int id = xU_counts.n_rows;
      xU_counts.insert_rows(id, key);
      row_map[key] = id;
      counts_of_counts_vec.resize(id + 1);
      counts_of_counts_vec(id) = 1;
      inv_idx_vec(i) = id;
    } else {
      inv_idx_vec(i) = it->second;
      counts_of_counts_vec(it->second)++;
    }
  }
  const int U = xU_counts.n_rows;

  // B. Calculate psi Sums
  double psi_dN_sum = 0, psi_dD_sum = 0, psi_FN_sum = 0, psi_FD_sum = 0;

  // U-statistic part
  arma::vec h_dN_u(U, arma::fill::zeros), h_dD_u(U, arma::fill::zeros);
  for(int u = 0; u < U; ++u) {
    arma::urowvec N_u = xU_counts.row(u);
    double r_u = arma::sum(N_u);
    h_dD_u(u) = r_u * (r_u - 1.0) / 2.0;
    // Sum of losses for subject u
    h_dN_u(u) = 0.5 * (arma::as_scalar(N_u * loss_matrix * N_u.t()) - arma::dot(N_u, loss_matrix.diag()));
  }
  psi_dN_sum = arma::dot(counts_of_counts_vec, h_dN_u);
  psi_dD_sum = arma::dot(counts_of_counts_vec, h_dD_u);

  // V-statistic part
  arma::mat kernel_FN_uv = xU_counts * loss_matrix * xU_counts.t();
  arma::vec r_u_vec = arma::conv_to<arma::vec>::from(arma::sum(xU_counts, 1));
  arma::mat kernel_FD_uv = r_u_vec * r_u_vec.t();

  arma::vec counts_d = arma::conv_to<arma::vec>::from(counts_of_counts_vec);
  psi_FN_sum = arma::as_scalar(counts_d.t() * kernel_FN_uv * counts_d);
  psi_FD_sum = arma::as_scalar(counts_d.t() * kernel_FD_uv * counts_d);

  // C. Calculate Point Estimates
  const double psi_dN_hat = psi_dN_sum / n;
  const double psi_dD_hat = psi_dD_sum / n;
  const double psi_FN_hat = psi_FN_sum / (n * n);
  const double psi_FD_hat = psi_FD_sum / (n * n);

  const double d_hat = (psi_dD_hat > 1e-9) ? psi_dN_hat / psi_dD_hat : 0.0;
  const double d_F_hat = (psi_FD_hat > 1e-9) ? psi_FN_hat / psi_FD_hat : 0.0;
  const double d_BP = arma::accu(loss_matrix) / (C * C);

  arma::vec estimates(2);
  estimates(0) = (d_F_hat > 1e-9) ? 1.0 - d_hat / d_F_hat : arma::datum::nan; // Fleiss
  estimates(1) = (d_BP > 1e-9) ? 1.0 - d_hat / d_BP : arma::datum::nan;       // BP

  // D. Calculate Influence Functions
  arma::vec phi_dN = h_dN_u.elem(inv_idx_vec) - psi_dN_hat;
  arma::vec phi_dD = h_dD_u.elem(inv_idx_vec) - psi_dD_hat;

  arma::vec g1_FN_u = (kernel_FN_uv * counts_d) / n;
  arma::vec g2_FN_u = (kernel_FN_uv.t() * counts_d) / n;
  arma::vec phi_FN = (g1_FN_u.elem(inv_idx_vec) - psi_FN_hat) + (g2_FN_u.elem(inv_idx_vec) - psi_FN_hat);

  arma::vec g1_FD_u = (kernel_FD_uv * counts_d) / n;
  arma::vec g2_FD_u = (kernel_FD_uv.t() * counts_d) / n;
  arma::vec phi_FD = (g1_FD_u.elem(inv_idx_vec) - psi_FD_hat) + (g2_FD_u.elem(inv_idx_vec) - psi_FD_hat);

  // E. Delta Method
  arma::mat phi_matrix(n, 4);
  phi_matrix.col(0) = phi_dN;
  phi_matrix.col(1) = phi_dD;
  phi_matrix.col(2) = phi_FN;
  phi_matrix.col(3) = phi_FD;

  arma::mat Gamma_hat = (1.0 / n) * phi_matrix.t() * phi_matrix;

  arma::mat J_d(2, 4, arma::fill::zeros);
  if (psi_dD_hat > 1e-9) { J_d(0, 0) = 1.0 / psi_dD_hat; J_d(0, 1) = -psi_dN_hat / std::pow(psi_dD_hat, 2); }
  if (psi_FD_hat > 1e-9) { J_d(1, 2) = 1.0 / psi_FD_hat; J_d(1, 3) = -psi_FN_hat / std::pow(psi_FD_hat, 2); }
  arma::mat disagreement_cov_matrix = (J_d * Gamma_hat * J_d.t()) / n;

  arma::mat J_kappa(2, 2, arma::fill::zeros);
  if (d_F_hat > 1e-9) { J_kappa(0, 0) = -1.0 / d_F_hat; J_kappa(0, 1) = d_hat / std::pow(d_F_hat, 2); }
  if (d_BP > 1e-9)    { J_kappa(1, 0) = -1.0 / d_BP; }

  arma::mat kappa_cov_matrix = J_kappa * disagreement_cov_matrix * J_kappa.t();

  // F. Return Result
  return {emdiscrete::Status::kOk, Estimation{estimates, kappa_cov_matrix}, ""};
}

Result<Estimation> kappa_continuous(
    const arma::mat& ratings,
    const loss::LossFunction& loss_func,
    bool use_ipw,
    bool use_gwet // overrides use_ipw.
) {
  // 1. Setup
  const int n = ratings.n_rows;
  const int R = ratings.n_cols;
  if (R < 2) return {emdiscrete::Status::kError, std::nullopt, "Input 'ratings' must have at least 2 raters (columns)."};

  arma::umat M(ratings.n_rows, ratings.n_cols, arma::fill::zeros);
  M.elem(arma::find_finite(ratings)).ones();

  if (M.n_elem == 0 || arma::accu(M) == 0) return {emdiscrete::Status::kError, std::nullopt, "All ratings are missing."};

  auto weights_result = calculate_inverse_weights(M, n, R, use_ipw, use_gwet);
  if (weights_result.status != emdiscrete::Status::kOk) {
    return {weights_result.status, std::nullopt, weights_result.error_message};
  }
  const auto& [pi_j_inv, pi_jk_inv] = *weights_result.value;

  // 3. Calculate Point Estimates of Psi Components (O(n^2) loop)
  double psi_dN_sum = 0, psi_dD_sum = 0, psi_CN_sum = 0,
    psi_CD_sum = 0, psi_FN_sum = 0, psi_FD_sum = 0;

  for (int i = 0; i < n; ++i) {
    // d_hat components (U-statistic kernel part)
    for (int j = 0; j < R - 1; ++j) {
      for (int k = j + 1; k < R; ++k) {
        if (M(i, j) && M(i, k)) {
          psi_dN_sum += loss_func(ratings(i, j), ratings(i, k)) * pi_jk_inv(j, k);
          psi_dD_sum += pi_jk_inv(j, k);
        }
      }
    }
    // d_C and d_F components (V-statistic kernel part)
    for (int ip = 0; ip < n; ++ip) {
      // Conger
      for (int j = 0; j < R - 1; ++j) {
        for (int k = j + 1; k < R; ++k) {
          if (M(i, j) && M(ip, k)) {
            double w = pi_j_inv(j) * pi_j_inv(k);
            psi_CN_sum += loss_func(ratings(i, j), ratings(ip, k)) * w;
            psi_CD_sum += w;
          }
        }
      }
      // Fleiss
      for (int j = 0; j < R; ++j) {
        for (int k = 0; k < R; ++k) {
          if (M(i, j) && M(ip, k)) {
            double w = pi_j_inv(j) * pi_j_inv(k);
            psi_FN_sum += loss_func(ratings(i, j), ratings(ip, k)) * w;
            psi_FD_sum += w;
          }
        }
      }
    }
  }

  const double psi_dN_hat = psi_dN_sum / n;
  const double psi_dD_hat = psi_dD_sum / n;
  const double psi_CN_hat = psi_CN_sum / (n * n);
  const double psi_CD_hat = psi_CD_sum / (n * n);
  const double psi_FN_hat = psi_FN_sum / (n * n);
  const double psi_FD_hat = psi_FD_sum / (n * n);

  // Calculate Disagreement and Kappa Estimates
  const double d_hat = (psi_dD_hat > 1e-9) ? psi_dN_hat / psi_dD_hat : 0.0;
  const double d_C_hat = (psi_CD_hat > 1e-9) ? psi_CN_hat / psi_CD_hat : 0.0;
  const double d_F_hat = (psi_FD_hat > 1e-9) ? psi_FN_hat / psi_FD_hat : 0.0;
  // No d_BP for continuous data

  arma::vec kappa_estimates(2);
  kappa_estimates(0) = (d_C_hat > 1e-9) ? 1.0 - d_hat / d_C_hat : arma::datum::nan; // Conger
  kappa_estimates(1) = (d_F_hat > 1e-9) ? 1.0 - d_hat / d_F_hat : arma::datum::nan; // Fleiss

  // 4. Calculate Influence Functions
  arma::vec h_dN(n, arma::fill::zeros), h_dD(n, arma::fill::zeros);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < R - 1; ++j) {
      for (int k = j + 1; k < R; ++k) {
        if (M(i, j) && M(i, k)) {
          h_dN(i) += loss_func(ratings(i, j), ratings(i, k)) * pi_jk_inv(j, k);
          h_dD(i) += pi_jk_inv(j, k);
        }
      }
    }
  }
  arma::vec phi_dN = h_dN - psi_dN_hat;
  arma::vec phi_dD = h_dD - psi_dD_hat;

  arma::mat kernel_CN(n, n), kernel_CD(n, n), kernel_FN(n, n), kernel_FD(n, n);
  for (int i = 0; i < n; ++i) {
    for (int ip = 0; ip < n; ++ip) {
      double h_cn = 0, h_cd = 0, h_fn = 0, h_fd = 0;
      // Conger kernel
      for (int j = 0; j < R - 1; ++j) {
        for (int k = j + 1; k < R; ++k) {
          if (M(i, j) && M(ip, k)) {
            double w = pi_j_inv(j) * pi_j_inv(k);
            h_cn += loss_func(ratings(i, j), ratings(ip, k)) * w;
            h_cd += w;
          }
        }
      }
      // Fleiss kernel
      for (int j = 0; j < R; ++j) {
        for (int k = 0; k < R; ++k) {
          if (M(i, j) && M(ip, k)) {
            double w = pi_j_inv(j) * pi_j_inv(k);
            h_fn += loss_func(ratings(i, j), ratings(ip, k)) * w;
            h_fd += w;
          }
        }
      }
      kernel_CN(i, ip) = h_cn; kernel_CD(i, ip) = h_cd;
      kernel_FN(i, ip) = h_fn; kernel_FD(i, ip) = h_fd;
    }
  }

  arma::vec h_CN1 = arma::mean(kernel_CN, 1);
  arma::vec h_CN2 = arma::mean(kernel_CN, 0).t();
  arma::vec h_CD1 = arma::mean(kernel_CD, 1);
  arma::vec h_CD2 = arma::mean(kernel_CD, 0).t();
  arma::vec h_FN1 = arma::mean(kernel_FN, 1);
  arma::vec h_FN2 = arma::mean(kernel_FN, 0).t();
  arma::vec h_FD1 = arma::mean(kernel_FD, 1);
  arma::vec h_FD2 = arma::mean(kernel_FD, 0).t();

  arma::vec phi_CN = (h_CN1 - psi_CN_hat) + (h_CN2 - psi_CN_hat);
  arma::vec phi_CD = (h_CD1 - psi_CD_hat) + (h_CD2 - psi_CD_hat);
  arma::vec phi_FN = (h_FN1 - psi_FN_hat) + (h_FN2 - psi_FN_hat);
  arma::vec phi_FD = (h_FD1 - psi_FD_hat) + (h_FD2 - psi_FD_hat);

  arma::mat phi_matrix(n, 6);
  phi_matrix.col(0) = phi_dN; phi_matrix.col(1) = phi_dD;
  phi_matrix.col(2) = phi_CN; phi_matrix.col(3) = phi_CD;
  phi_matrix.col(4) = phi_FN; phi_matrix.col(5) = phi_FD;

  arma::mat Gamma_hat = (1.0 / n) * phi_matrix.t() * phi_matrix;

  // 5. Delta Method 1: Asymptotic Covariance of Disagreements (Sigma_hat)
  arma::mat J_d(3, 6, arma::fill::zeros);
  if (psi_dD_hat > 1e-9) { J_d(0, 0) = 1.0 / psi_dD_hat; J_d(0, 1) = -psi_dN_hat / (psi_dD_hat * psi_dD_hat); }
  if (psi_CD_hat > 1e-9) { J_d(1, 2) = 1.0 / psi_CD_hat; J_d(1, 3) = -psi_CN_hat / (psi_CD_hat * psi_CD_hat); }
  if (psi_FD_hat > 1e-9) { J_d(2, 4) = 1.0 / psi_FD_hat; J_d(2, 5) = -psi_FN_hat / (psi_FD_hat * psi_FD_hat); }
  arma::mat Sigma_hat = J_d * Gamma_hat * J_d.t();

  // 6. Delta Method 2: Asymptotic Covariance of Kappas
  arma::mat J_kappa(2, 3, arma::fill::zeros);
  if (d_C_hat > 1e-9) { J_kappa(0, 0) = -1.0 / d_C_hat; J_kappa(0, 1) = d_hat / (d_C_hat * d_C_hat); }
  if (d_F_hat > 1e-9) { J_kappa(1, 0) = -1.0 / d_F_hat; J_kappa(1, 2) = d_hat / (d_F_hat * d_F_hat); }

  arma::mat kappa_cov_asymptotic = J_kappa * Sigma_hat * J_kappa.t();
  arma::mat kappa_cov_matrix = kappa_cov_asymptotic / n;

  return {emdiscrete::Status::kOk, Estimation{kappa_estimates, kappa_cov_matrix}, ""};
}

} // namespace kappanp
} // namespace misskappa
