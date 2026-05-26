// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <quadagree.hpp>

// --- Helper functions to convert R strings to C++ enums ---
quadagree::Transform StringToTransform(Rcpp::String s) {
  if (s == "fisher") return quadagree::Transform::kFisher;
  if (s == "log") return quadagree::Transform::kLog;
  if (s == "arcsin") return quadagree::Transform::kArcsin;
  return quadagree::Transform::kNone;
}

quadagree::Alternative StringToAlternative(Rcpp::String s) {
  if (s == "greater") return quadagree::Alternative::kGreater;
  if (s == "less") return quadagree::Alternative::kLess;
  return quadagree::Alternative::kTwoSided;
}

//' @title Raw Kappa Calculation (C++ Core)
//' @description Computes estimates, CIs, and covariance for Fleiss', Conger's, and BP kappa.
//' @param c1 The constant for Brennan-Prediger chance agreement.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List KappasRawCpp(const arma::mat& x, double c1,
                       Rcpp::String transform_str = "none",
                       double conf_level = 0.95,
                       Rcpp::String alternative_str = "two.sided",
                       bool bootstrap = false,
                       int n_reps = 1000,
                       Rcpp::Nullable<int> seed = R_NilValue) {
 std::mt19937 rng;
 if (seed.isNotNull()) { rng.seed(Rcpp::as<int>(seed)); }
 else { std::random_device rd; rng.seed(rd()); }

 quadagree::Result<quadagree::AllCIResults> res = quadagree::QuadagreeRaw(
   x, c1, StringToTransform(transform_str), conf_level,
   StringToAlternative(alternative_str), bootstrap, n_reps, rng);

 if (res.status != quadagree::Status::kOk) {
   Rcpp::stop("C++ core calculation failed.");
 }

 const auto& val = *res.value;
 Rcpp::NumericVector interval = Rcpp::NumericVector::create(val.fleiss.conf_low, val.fleiss.conf_high);
 Rcpp::List fleiss_list = Rcpp::List::create(Rcpp::_["estimate"] = val.fleiss.estimate,
                                             Rcpp::_["std_err"] = val.fleiss.std_err,
                                             Rcpp::_["interval"] = interval,
                                             Rcpp::_["n_eff"] = val.fleiss.n_eff);
 interval = Rcpp::NumericVector::create(val.conger.conf_low, val.conger.conf_high);
 Rcpp::List conger_list = Rcpp::List::create(Rcpp::_["estimate"] = val.conger.estimate,
                                             Rcpp::_["std_err"] = val.conger.std_err,
                                             Rcpp::_["interval"] = interval,
                                             Rcpp::_["n_eff"] = val.conger.n_eff);
 interval = Rcpp::NumericVector::create(val.bp.conf_low, val.bp.conf_high);
 Rcpp::List bp_list = Rcpp::List::create(Rcpp::_["estimate"] = val.bp.estimate,
                                         Rcpp::_["std_err"] = val.bp.std_err,
                                         Rcpp::_["interval"] = interval,
                                         Rcpp::_["n_eff"] = val.bp.n_eff);

 return Rcpp::List::create(Rcpp::_["fleiss"] = fleiss_list, Rcpp::_["conger"] = conger_list, Rcpp::_["bp"] = bp_list, Rcpp::_["scaled_acov"] = val.scaled_acov);
}

//' @title Aggregated Kappa Calculation (C++ Core)
//' @description Computes estimates and CIs for Fleiss' and BP kappa from aggregated data.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List KappasAggrCpp(const arma::mat& x, const arma::vec& values, int R, double c1,
                        Rcpp::String transform_str = "none",
                        double conf_level = 0.95,
                        Rcpp::String alternative_str = "two.sided",
                        bool bootstrap = false,
                        int n_reps = 1000,
                        Rcpp::Nullable<int> seed = R_NilValue) {
 std::mt19937 rng;
 if (seed.isNotNull()) { rng.seed(Rcpp::as<int>(seed)); }
 else { std::random_device rd; rng.seed(rd()); }

 quadagree::Result<quadagree::AggrCIResults> res = quadagree::QuadagreeAggr(
   x, values, R, c1, StringToTransform(transform_str), conf_level,
   StringToAlternative(alternative_str), bootstrap, n_reps, rng);

 if (res.status != quadagree::Status::kOk) {
   Rcpp::stop("C++ core calculation failed.");
 }

 const auto& val = *res.value;
 Rcpp::NumericVector interval = Rcpp::NumericVector::create(val.fleiss.conf_low, val.fleiss.conf_high);
 Rcpp::List fleiss_list = Rcpp::List::create(Rcpp::_["estimate"] = val.fleiss.estimate,
                                             Rcpp::_["std_err"] = val.fleiss.std_err,
                                             Rcpp::_["interval"] = interval,
                                             Rcpp::_["n_eff"] = val.fleiss.n_eff);

 interval = Rcpp::NumericVector::create(val.bp.conf_low, val.bp.conf_high);
 Rcpp::List bp_list = Rcpp::List::create(Rcpp::_["estimate"] = val.bp.estimate,
                                         Rcpp::_["std_err"] = val.bp.std_err,
                                         Rcpp::_["interval"] = interval,
                                         Rcpp::_["n_eff"] = val.bp.n_eff);

 return Rcpp::List::create(Rcpp::_["fleiss"] = fleiss_list, Rcpp::_["bp"] = bp_list);
}
