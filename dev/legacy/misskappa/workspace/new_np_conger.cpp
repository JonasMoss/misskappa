// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//' @title Calculate Conger's Kappa for Discrete Data with Missingness
//' @description
//' A high-performance C++ implementation to calculate Conger's kappa, its
//' components (observed and chance disagreement), and their joint asymptotic
//' covariance matrix using influence functions. This function is designed for
//' discrete (binned) data and handles missing values under the MCAR assumption.
//'
//' @param ratings_mat An integer matrix of ratings (subjects x raters).
//'   Ratings should be 1-indexed (e.g., 1, 2, 3, ...). `NA`s are treated as missing.
//' @param loss_mat A numeric matrix representing the loss function. `loss_mat(c1, c2)`
//'   is the disagreement between category `c1` and `c2`. For 0-1 loss, this is
//'   a matrix of 1s with a 0 diagonal.
//'
//' @return A list containing:
  //'   \item{estimates}{A named vector with estimates for d (observed disagreement),
//'     dC (Conger's chance disagreement), and kappaC (Conger's kappa).}
//'   \item{vcov}{The 3x3 asymptotic covariance matrix for (d, dC, kappaC).}
//'   \item{influence_functions}{A matrix of the estimated influence functions for (d, dC) for each subject.}
//' @export
// [[Rcpp::export]]
Rcpp::List conger_kappa_discrete(const arma::imat& ratings_mat, const arma::mat& loss_mat) {
    // ---- SETUP ----
    int n_subjects = ratings_mat.n_rows;
    int n_raters = ratings_mat.n_cols;
    int n_categories = loss_mat.n_rows;

    // Create a boolean matrix for missingness (M)
    // and a zero-indexed rating matrix for C++
    arma::umat M(n_subjects, n_raters, arma::fill::zeros);
    arma::imat Y = ratings_mat; // Copy to modify
    for (int i = 0; i < n_subjects; ++i) {
        for (int j = 0; j < n_raters; ++j) {
            if (arma::is_finite(Y(i, j))) {
                M(i, j) = 1;
                Y(i, j) -= 1; // Convert to 0-indexed
            }
        }
    }

    // ---- PART A: PRE-COMPUTATION (NUISANCE PARAMETERS) ----

    // 1. Rater-specific observation counts and marginal probabilities
    arma::vec n_j(n_raters);
    arma::mat p_jc(n_raters, n_categories);
    for (int j = 0; j < n_raters; ++j) {
        n_j(j) = arma::accu(M.col(j));
        if (n_j(j) == 0) {
            Rcpp::stop("Rater %d has no observations.", j + 1);
        }
        for (int c = 0; c < n_categories; ++c) {
            double count_jc = 0;
            for (int i = 0; i < n_subjects; ++i) {
                if (M(i, j) && Y(i, j) == c) {
                    count_jc++;
                }
            }
            p_jc(j, c) = count_jc / n_j(j);
        }
    }

    // 2. Pairwise quantities
    int n_pairs = n_raters * (n_raters - 1) / 2;
    arma::vec d_jk_hat(n_pairs);
    arma::vec pi_jk_hat(n_pairs);

    int pair_idx = 0;
    for (int j = 0; j < n_raters; ++j) {
        for (int k = j + 1; k < n_raters; ++k) {
            double num_d_jk = 0;
            double den_d_jk = 0;
            for (int i = 0; i < n_subjects; ++i) {
                if (M(i, j) && M(i, k)) {
                    num_d_jk += loss_mat(Y(i, j), Y(i, k));
                    den_d_jk++;
                }
            }
            d_jk_hat(pair_idx) = (den_d_jk > 0) ? num_d_jk / den_d_jk : 0;
            pi_jk_hat(pair_idx) = den_d_jk / n_subjects;
            pair_idx++;
        }
    }

    // 3. Rater-specific expected losses (lambda)
    arma::mat lambda_jc(n_raters, n_categories);
    for (int j = 0; j < n_raters; ++j) {
        for (int c = 0; c < n_categories; ++c) {
            double expected_loss = 0;
            for (int b = 0; b < n_categories; ++b) {
                expected_loss += loss_mat(c, b) * p_jc(j, b);
            }
            lambda_jc(j, c) = expected_loss;
        }
    }

    // ---- PART B: INFLUENCE FUNCTION LOOP ----

    arma::mat psi_matrix(n_subjects, 2, arma::fill::zeros); // Columns for d, dC

    for (int i = 0; i < n_subjects; ++i) {
        // --- Calculate Fundamental Influence Values for subject i ---
        arma::vec psi_p_ic(n_raters * n_categories, arma::fill::zeros);
        for (int j = 0; j < n_raters; ++j) {
            if (M(i, j)) {
                for (int c = 0; c < n_categories; ++c) {
                    double indicator = (Y(i, j) == c) ? 1.0 : 0.0;
                    psi_p_ic(j * n_categories + c) = (indicator - p_jc(j, c)) / (n_j(j) / n_subjects);
                }
            }
        }

        // --- Assemble Final Influence Vector for subject i ---
        double psi_d_i = 0;
        double psi_dC_i = 0;

        pair_idx = 0;
        for (int j = 0; j < n_raters; ++j) {
            for (int k = j + 1; k < n_raters; ++k) {
                // For observed disagreement d
                if (M(i, j) && M(i, k)) {
                    double loss_ijk = loss_mat(Y(i, j), Y(i, k));
                    psi_d_i += (loss_ijk - d_jk_hat(pair_idx)) / pi_jk_hat(pair_idx);
                }

                // For chance disagreement dC
                double psi_dF_jk_i = 0;
                for (int c = 0; c < n_categories; ++c) {
                    psi_dF_jk_i += lambda_jc(k, c) * psi_p_ic[j * n_categories + c];
                    psi_dF_jk_i += lambda_jc(j, c) * psi_p_ic[k * n_categories + c];
                }
                psi_dC_i += psi_dF_jk_i;

                pair_idx++;
            }
        }

        psi_matrix(i, 0) = psi_d_i / n_pairs;
        psi_matrix(i, 1) = psi_dC_i / n_pairs;
    }

    // ---- PART C: FINAL ASSEMBLY ----

    // 1. Point Estimates
    double d_hat = arma::mean(d_jk_hat);
    double dC_hat = 0;
    pair_idx = 0;
    for(int j=0; j<n_raters; ++j){
      for(int k=j+1; k<n_raters; ++k){
        double dF_jk_hat = 0;
        for(int c1=0; c1<n_categories; ++c1){
          for(int c2=0; c2<n_categories; ++c2){
            dF_jk_hat += loss_mat(c1,c2) * p_jc(j,c1) * p_jc(k,c2);
          }
        }
        dC_hat += dF_jk_hat;
        pair_idx++;
      }
    }
    dC_hat /= n_pairs;

    double kappaC_hat = 1.0 - d_hat / dC_hat;

    Rcpp::NumericVector estimates = Rcpp::NumericVector::create(
        Rcpp::Named("d") = d_hat,
        Rcpp::Named("dC") = dC_hat,
        Rcpp::Named("kappaC") = kappaC_hat
    );

    // 2. Covariance Matrix of (d, dC)
    arma::mat vcov_d_dC = arma::cov(psi_matrix) / n_subjects;

    // 3. Delta Method for Kappa
    arma::vec grad(2);
    grad(0) = -1.0 / dC_hat;        // d(kappa)/d(d)
    grad(1) = d_hat / (dC_hat * dC_hat); // d(kappa)/d(dC)

    double var_kappaC = arma::as_scalar(grad.t() * vcov_d_dC * grad);

    // 4. Assemble final 3x3 vcov matrix
    arma::mat vcov_final(3, 3, arma::fill::zeros);
    vcov_final.submat(0, 0, 1, 1) = vcov_d_dC;

    // Cov(d, kappaC) and Cov(dC, kappaC) via delta method
    arma::vec cov_d_kappa = vcov_d_dC * grad;
    vcov_final(0, 2) = cov_d_kappa(0);
    vcov_final(2, 0) = cov_d_kappa(0);
    vcov_final(1, 2) = cov_d_kappa(1);
    vcov_final(2, 1) = cov_d_kappa(1);
    vcov_final(2, 2) = var_kappaC;

    return Rcpp::List::create(
        Rcpp::Named("estimates") = estimates,
        Rcpp::Named("vcov") = vcov_final,
        Rcpp::Named("influence_functions") = psi_matrix
    );
}
