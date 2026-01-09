// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "misskappa.h" // The single public header for our kappa library
#include <set>
#include <map>

// ... rcpp_generate_loss_matrix helper remains unchanged ...
arma::mat rcpp_generate_loss_matrix(
    std::string weight_type, int c, const arma::vec& values
) {
  std::optional<arma::vec> values_opt = values;

  misskappa::Result<arma::mat> w_res;
  if (weight_type == "identity" || weight_type == "unweighted") w_res = misskappa::loss::identity_weights(c);
  else if (weight_type == "linear")    w_res = misskappa::loss::linear_weights(c, values_opt.value());
  else if (weight_type == "quadratic") w_res = misskappa::loss::quadratic_weights(c, values_opt.value());
  else if (weight_type == "ordinal")   w_res = misskappa::loss::ordinal_weights(c);
  else if (weight_type == "radical")   w_res = misskappa::loss::radical_weights(c, values_opt.value());
  else if (weight_type == "ratio")     w_res = misskappa::loss::ratio_weights(c, values_opt.value());
  else if (weight_type == "circular")  w_res = misskappa::loss::circular_weights(c, values_opt.value());
  else if (weight_type == "bipolar")   w_res = misskappa::loss::bipolar_weights(c, values_opt.value());
  else Rcpp::stop("Unknown weight type: " + weight_type);

  if (!w_res.IsOk()) Rcpp::stop(w_res.error_message);
  return arma::ones<arma::mat>(c, c) - w_res.value.value();
}


//' Unified C++ Backend for Continuous Data Kappa
//' @keywords internal
// [[Rcpp::export]]
 Rcpp::List unified_kappa_continuous_rcpp(
     const Rcpp::NumericMatrix& x_r,
     std::string method,
     std::string weight_type
 ) {
   arma::mat x_scores = Rcpp::as<arma::mat>(x_r);
   misskappa::Result<misskappa::Estimation> kappa_res;

   if (method == "quadratic") {
     arma::vec dummy_values; // Not needed for kappaqp's continuous method
     kappa_res = misskappa::kappaqp::kappa(x_scores, dummy_values);
   } else if (method == "available" || method == "ipw") {
     // This block does NOT perform binning.
     // It calculates the data range solely to parameterize the loss function.
     arma::mat x_finite = x_scores;
     x_finite.elem(arma::find_nonfinite(x_finite)).zeros();
     double min_val = x_finite.min(), max_val = x_finite.max();

     // Create the appropriate continuous loss FUNCTION from the factory.
     misskappa::Result<misskappa::loss::LossFunction> loss_factory_res;
     if(weight_type == "identity") loss_factory_res = misskappa::loss::create_identity_loss();
     else if(weight_type == "linear") loss_factory_res = misskappa::loss::create_linear_loss(min_val, max_val);
     else if(weight_type == "quadratic") loss_factory_res = misskappa::loss::create_quadratic_loss(min_val, max_val);
     else if(weight_type == "radical") loss_factory_res = misskappa::loss::create_radical_loss(min_val, max_val);
     else if(weight_type == "ratio") loss_factory_res = misskappa::loss::create_ratio_loss(min_val, max_val);
     else Rcpp::stop("Weight type '" + weight_type + "' is not supported for continuous 'np'/'ipw' method.");

     if(!loss_factory_res.IsOk()) Rcpp::stop(loss_factory_res.error_message);

     // Call the continuous non-parametric backend with the UNMODIFIED data
     // and the generated loss function.
     kappa_res = misskappa::kappanp::kappa_continuous(
       x_scores,
       loss_factory_res.value.value(),
       (method == "ipw"),
       (method == "gwet")
     );

   } else {
     Rcpp::stop("Unknown method for continuous data: " + method);
   }

   if (!kappa_res.IsOk()) Rcpp::stop(kappa_res.error_message);

   arma::vec all_estimates = kappa_res.value.value().estimates;
   arma::mat all_vcov = kappa_res.value.value().vcov;

   arma::vec final_estimates;
   arma::mat final_vcov;

   if (method == "quadratic") {
     // kappaqp returns [Conger, Fleiss, BP], we want [Conger, Fleiss]
     final_estimates = all_estimates.head(2);
     final_vcov = all_vcov.submat(0, 0, 1, 1);
   } else {
     // kappanp_continuous returns [Conger, Fleiss]
     final_estimates = all_estimates;
     final_vcov = all_vcov;
   }

   return Rcpp::List::create(
     Rcpp::_["estimates"] = final_estimates,
     Rcpp::_["vcov"] = final_vcov
   );
 }


//' Unified C++ Backend for Raw Categorical Data Kappa
//' @keywords internal
// [[Rcpp::export]]
 Rcpp::List unified_kappa_raw_rcpp(
     const Rcpp::IntegerMatrix& x_r,
     std::string method,
     std::string weight_type,
     Rcpp::Nullable<Rcpp::NumericVector> values,
     Rcpp::List options
 ) {
   std::set<int> unique_cats_set;
   for(int i=0; i<x_r.nrow(); ++i) for(int j=0; j<x_r.ncol(); ++j) {
     if(x_r(i,j) != NA_INTEGER) unique_cats_set.insert(x_r(i,j));
   }
   std::vector<int> unique_cats(unique_cats_set.begin(), unique_cats_set.end());
   int c = unique_cats.size();
   if (c == 0) Rcpp::stop("No valid ratings provided.");

   arma::vec cat_values_arma;
   if (values.isNotNull()) {
     cat_values_arma = Rcpp::as<arma::vec>(values);
     if (static_cast<int>(cat_values_arma.n_elem) != c) {
       Rcpp::stop("Length of 'values' does not match the number of unique categories in 'x'.");
     }
   } else {
     cat_values_arma.set_size(c);
     for(int i=0; i<c; ++i) cat_values_arma(i) = unique_cats[i];
   }

   arma::mat loss = rcpp_generate_loss_matrix(weight_type, c, cat_values_arma);

   misskappa::Result<misskappa::Estimation> kappa_res;

   arma::imat x_arma_original(x_r.nrow(), x_r.ncol());
   for(int i=0; i<x_r.nrow(); ++i) for(int j=0; j<x_r.ncol(); ++j) {
     x_arma_original(i,j) = (x_r(i,j) == NA_INTEGER) ? emdiscrete::kNaInteger : x_r(i,j);
   }

   if (method == "ml") {
     std::map<int, int> cat_to_idx;
     for(int i=0; i<c; ++i) cat_to_idx[unique_cats[i]] = i;
     arma::imat x_arma_indexed(x_r.nrow(), x_r.ncol());
     for(int i=0; i<x_r.nrow(); ++i) for(int j=0; j<x_r.ncol(); ++j) {
       if(x_r(i,j) == NA_INTEGER) {
         x_arma_indexed(i,j) = emdiscrete::kNaInteger;
       } else {
         x_arma_indexed(i,j) = cat_to_idx.at(x_r(i,j));
       }
     }

     emdiscrete::EM_Options em_opts;
     Rcpp::List em_opts_r = Rcpp::as<Rcpp::List>(options["em_options"]);
     em_opts.tol = Rcpp::as<double>(em_opts_r["tol"]);
     em_opts.max_iter = Rcpp::as<int>(em_opts_r["max_iter"]);

     auto prep_res = emdiscrete::preprocess_raw(x_arma_indexed, c);
     if (!prep_res.IsOk()) Rcpp::stop(prep_res.error_message);

     auto em_res = emdiscrete::run_em(prep_res.value.value(), em_opts);
     if (!em_res.IsOk()) Rcpp::stop(em_res.error_message);

     kappa_res = misskappa::kappaml::kappa(em_res.value.value(), loss);

   } else if (method == "ipw" || method == "available" || method == "gwet") {
     kappa_res = misskappa::kappanp::kappa(x_arma_original, loss, (method == "ipw"), (method == "gwet"));

   } else if (method == "quadratic") {
     std::map<int, int> cat_to_idx;
     for(int i=0; i<c; ++i) cat_to_idx[unique_cats[i]] = i;
     arma::mat x_scores(x_r.nrow(), x_r.ncol());
     for(int i=0; i<x_r.nrow(); ++i) for(int j=0; j<x_r.ncol(); ++j) {
       if (x_r(i,j) == NA_INTEGER) {
         x_scores(i,j) = arma::datum::nan;
       } else {
         x_scores(i,j) = cat_values_arma(cat_to_idx.at(x_r(i,j)));
       }
     }
     kappa_res = misskappa::kappaqp::kappa(x_scores, cat_values_arma);

   } else {
     Rcpp::stop("Unknown method: " + method);
   }

   if (!kappa_res.IsOk()) Rcpp::stop(kappa_res.error_message);
   return Rcpp::List::create(
     Rcpp::_["estimates"] = kappa_res.value.value().estimates,
     Rcpp::_["vcov"] = kappa_res.value.value().vcov
   );
 }


//' Unified C++ Backend for Counts Data Kappa
//' @keywords internal
// [[Rcpp::export]]
 Rcpp::List unified_kappa_counts_rcpp(
     const Rcpp::IntegerMatrix& x_r, int r, std::string method,
     std::string weight_type, Rcpp::Nullable<Rcpp::NumericVector> values, Rcpp::List options
 ) {
   arma::umat x_arma = Rcpp::as<arma::umat>(x_r);
   int c = x_arma.n_cols;

   arma::vec cat_values_arma;
   if (values.isNotNull()) {
     cat_values_arma = Rcpp::as<arma::vec>(values);
   } else {
     cat_values_arma = arma::regspace(1, c);
   }
   if(static_cast<int>(cat_values_arma.n_elem) != c) Rcpp::stop("Length of 'values' does not match number of categories.");

   arma::mat loss = rcpp_generate_loss_matrix(weight_type, c, cat_values_arma);
   misskappa::Result<misskappa::Estimation> kappa_res;

   if (method == "ml") {
     emdiscrete::EM_Options em_opts;
     Rcpp::List em_opts_r = Rcpp::as<Rcpp::List>(options["em_options"]);
     em_opts.tol = Rcpp::as<double>(em_opts_r["tol"]);
     em_opts.max_iter = Rcpp::as<int>(em_opts_r["max_iter"]);

     auto prep_res = emdiscrete::preprocess_counts(x_arma, r);
     if (!prep_res.IsOk()) Rcpp::stop(prep_res.error_message);

     auto em_res = emdiscrete::run_em(prep_res.value.value(), em_opts);
     if (!em_res.IsOk()) Rcpp::stop(em_res.error_message);

     kappa_res = misskappa::kappaml::kappa_counts(em_res.value.value(), loss);
   } else if (method == "quadratic") {
     kappa_res = misskappa::kappaqp::kappa_counts(Rcpp::as<arma::mat>(x_r), cat_values_arma, r);
   } else if (method == "available") {
     // --- NEW CASE ---
     kappa_res = misskappa::kappanp::kappa_counts(x_arma, loss);
   } else {
     Rcpp::stop("Unknown or unsupported method for counts data: " + method);
   }

   if (!kappa_res.IsOk()) Rcpp::stop(kappa_res.error_message);
   return Rcpp::List::create(
     Rcpp::_["estimates"] = kappa_res.value.value().estimates,
     Rcpp::_["vcov"] = kappa_res.value.value().vcov
   );
 }

