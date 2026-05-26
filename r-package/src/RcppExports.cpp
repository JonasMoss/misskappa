// Hand-written counterpart to Rcpp::compileAttributes output, so the
// package builds in one shot without invoking compileAttributes.

#include <Rcpp.h>

using namespace Rcpp;

// Forward declaration of the actual implementation in rcpp_glue.cpp.
Rcpp::List rcpp_kappa_raw(
    const Rcpp::IntegerMatrix& x,
    std::string method,
    std::string weight_type,
    Rcpp::Nullable<Rcpp::NumericVector> values,
    Rcpp::List em_options);

RcppExport SEXP _misskappa_rcpp_kappa_raw(
    SEXP xSEXP, SEXP methodSEXP, SEXP weight_typeSEXP,
    SEXP valuesSEXP, SEXP em_optionsSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter<const Rcpp::IntegerMatrix&>::type x(xSEXP);
  Rcpp::traits::input_parameter<std::string>::type method(methodSEXP);
  Rcpp::traits::input_parameter<std::string>::type weight_type(weight_typeSEXP);
  Rcpp::traits::input_parameter<Rcpp::Nullable<Rcpp::NumericVector>>::type values(valuesSEXP);
  Rcpp::traits::input_parameter<Rcpp::List>::type em_options(em_optionsSEXP);
  rcpp_result_gen = Rcpp::wrap(rcpp_kappa_raw(x, method, weight_type, values, em_options));
  return rcpp_result_gen;
  END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_misskappa_rcpp_kappa_raw", (DL_FUNC)&_misskappa_rcpp_kappa_raw, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_misskappa(DllInfo* dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
