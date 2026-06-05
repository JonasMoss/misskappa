// Hand-written counterpart to Rcpp::compileAttributes output, so the
// package builds in one shot without invoking compileAttributes.

#include <Rcpp.h>

using namespace Rcpp;

// Forward declarations of the actual implementations in rcpp_glue.cpp.
Rcpp::List rcpp_kappa_raw(
    const Rcpp::IntegerMatrix& x,
    std::string method,
    std::string weight_type,
    Rcpp::Nullable<Rcpp::NumericVector> values,
    Rcpp::List em_options);

Rcpp::List rcpp_kappa_continuous(
    const Rcpp::NumericMatrix& x,
    std::string method,
    std::string weight_type);

Rcpp::List rcpp_alpha_raw(
    const Rcpp::IntegerMatrix& x,
    std::string method,
    Rcpp::Nullable<Rcpp::NumericVector> values,
    Rcpp::List em_options);

Rcpp::List rcpp_kappa_counts(
    const Rcpp::IntegerMatrix& x,
    std::string weight_type,
    Rcpp::Nullable<Rcpp::NumericVector> values);

Rcpp::List rcpp_kappa_quadratic(
    const Rcpp::NumericMatrix& x,
    Rcpp::NumericVector values,
    std::string vcov_type,
    double relative_kurtosis);

Rcpp::List rcpp_kappa_quadratic_counts(
    const Rcpp::IntegerMatrix& x,
    Rcpp::NumericVector values,
    int r_total);

Rcpp::List rcpp_kappa_fiml_counts(
    const Rcpp::IntegerMatrix& x,
    std::string weight_type,
    Rcpp::Nullable<Rcpp::NumericVector> values,
    int r_total,
    Rcpp::List em_options);

Rcpp::List rcpp_fiml_louis_spectrum(
    const Rcpp::IntegerMatrix& x,
    std::string weight_type,
    Rcpp::Nullable<Rcpp::NumericVector> values,
    Rcpp::List em_options);

Rcpp::List rcpp_kappa_gwise_categorical(
    const Rcpp::IntegerMatrix& x,
    std::string distance_type,
    int g,
    int max_chance_tuples);

Rcpp::List rcpp_kappa_gwise_continuous(
    const Rcpp::NumericMatrix& x,
    std::string distance_type,
    int g,
    int max_chance_tuples);

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

RcppExport SEXP _misskappa_rcpp_kappa_continuous(
    SEXP xSEXP, SEXP methodSEXP, SEXP weight_typeSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter<const Rcpp::NumericMatrix&>::type x(xSEXP);
  Rcpp::traits::input_parameter<std::string>::type method(methodSEXP);
  Rcpp::traits::input_parameter<std::string>::type weight_type(weight_typeSEXP);
  rcpp_result_gen = Rcpp::wrap(rcpp_kappa_continuous(x, method, weight_type));
  return rcpp_result_gen;
  END_RCPP
}

RcppExport SEXP _misskappa_rcpp_alpha_raw(
    SEXP xSEXP, SEXP methodSEXP, SEXP valuesSEXP, SEXP em_optionsSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter<const Rcpp::IntegerMatrix&>::type x(xSEXP);
  Rcpp::traits::input_parameter<std::string>::type method(methodSEXP);
  Rcpp::traits::input_parameter<Rcpp::Nullable<Rcpp::NumericVector>>::type values(valuesSEXP);
  Rcpp::traits::input_parameter<Rcpp::List>::type em_options(em_optionsSEXP);
  rcpp_result_gen = Rcpp::wrap(rcpp_alpha_raw(x, method, values, em_options));
  return rcpp_result_gen;
  END_RCPP
}

RcppExport SEXP _misskappa_rcpp_kappa_counts(
    SEXP xSEXP, SEXP weight_typeSEXP, SEXP valuesSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter<const Rcpp::IntegerMatrix&>::type x(xSEXP);
  Rcpp::traits::input_parameter<std::string>::type weight_type(weight_typeSEXP);
  Rcpp::traits::input_parameter<Rcpp::Nullable<Rcpp::NumericVector>>::type values(valuesSEXP);
  rcpp_result_gen = Rcpp::wrap(rcpp_kappa_counts(x, weight_type, values));
  return rcpp_result_gen;
  END_RCPP
}

RcppExport SEXP _misskappa_rcpp_kappa_quadratic(
    SEXP xSEXP, SEXP valuesSEXP, SEXP vcov_typeSEXP, SEXP relative_kurtosisSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter<const Rcpp::NumericMatrix&>::type x(xSEXP);
  Rcpp::traits::input_parameter<Rcpp::NumericVector>::type values(valuesSEXP);
  Rcpp::traits::input_parameter<std::string>::type vcov_type(vcov_typeSEXP);
  Rcpp::traits::input_parameter<double>::type relative_kurtosis(relative_kurtosisSEXP);
  rcpp_result_gen = Rcpp::wrap(rcpp_kappa_quadratic(x, values, vcov_type, relative_kurtosis));
  return rcpp_result_gen;
  END_RCPP
}

RcppExport SEXP _misskappa_rcpp_kappa_quadratic_counts(
    SEXP xSEXP, SEXP valuesSEXP, SEXP r_totalSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter<const Rcpp::IntegerMatrix&>::type x(xSEXP);
  Rcpp::traits::input_parameter<Rcpp::NumericVector>::type values(valuesSEXP);
  Rcpp::traits::input_parameter<int>::type r_total(r_totalSEXP);
  rcpp_result_gen = Rcpp::wrap(rcpp_kappa_quadratic_counts(x, values, r_total));
  return rcpp_result_gen;
  END_RCPP
}

RcppExport SEXP _misskappa_rcpp_kappa_fiml_counts(
    SEXP xSEXP, SEXP weight_typeSEXP, SEXP valuesSEXP,
    SEXP r_totalSEXP, SEXP em_optionsSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter<const Rcpp::IntegerMatrix&>::type x(xSEXP);
  Rcpp::traits::input_parameter<std::string>::type weight_type(weight_typeSEXP);
  Rcpp::traits::input_parameter<Rcpp::Nullable<Rcpp::NumericVector>>::type values(valuesSEXP);
  Rcpp::traits::input_parameter<int>::type r_total(r_totalSEXP);
  Rcpp::traits::input_parameter<Rcpp::List>::type em_options(em_optionsSEXP);
  rcpp_result_gen = Rcpp::wrap(rcpp_kappa_fiml_counts(x, weight_type, values, r_total, em_options));
  return rcpp_result_gen;
  END_RCPP
}

RcppExport SEXP _misskappa_rcpp_fiml_louis_spectrum(
    SEXP xSEXP, SEXP weight_typeSEXP, SEXP valuesSEXP, SEXP em_optionsSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter<const Rcpp::IntegerMatrix&>::type x(xSEXP);
  Rcpp::traits::input_parameter<std::string>::type weight_type(weight_typeSEXP);
  Rcpp::traits::input_parameter<Rcpp::Nullable<Rcpp::NumericVector>>::type values(valuesSEXP);
  Rcpp::traits::input_parameter<Rcpp::List>::type em_options(em_optionsSEXP);
  rcpp_result_gen = Rcpp::wrap(rcpp_fiml_louis_spectrum(x, weight_type, values, em_options));
  return rcpp_result_gen;
  END_RCPP
}

RcppExport SEXP _misskappa_rcpp_kappa_gwise_categorical(
    SEXP xSEXP, SEXP distance_typeSEXP, SEXP gSEXP, SEXP max_chance_tuplesSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter<const Rcpp::IntegerMatrix&>::type x(xSEXP);
  Rcpp::traits::input_parameter<std::string>::type distance_type(distance_typeSEXP);
  Rcpp::traits::input_parameter<int>::type g(gSEXP);
  Rcpp::traits::input_parameter<int>::type max_chance_tuples(max_chance_tuplesSEXP);
  rcpp_result_gen = Rcpp::wrap(rcpp_kappa_gwise_categorical(x, distance_type, g, max_chance_tuples));
  return rcpp_result_gen;
  END_RCPP
}

RcppExport SEXP _misskappa_rcpp_kappa_gwise_continuous(
    SEXP xSEXP, SEXP distance_typeSEXP, SEXP gSEXP, SEXP max_chance_tuplesSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter<const Rcpp::NumericMatrix&>::type x(xSEXP);
  Rcpp::traits::input_parameter<std::string>::type distance_type(distance_typeSEXP);
  Rcpp::traits::input_parameter<int>::type g(gSEXP);
  Rcpp::traits::input_parameter<int>::type max_chance_tuples(max_chance_tuplesSEXP);
  rcpp_result_gen = Rcpp::wrap(rcpp_kappa_gwise_continuous(x, distance_type, g, max_chance_tuples));
  return rcpp_result_gen;
  END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_misskappa_rcpp_kappa_raw", (DL_FUNC)&_misskappa_rcpp_kappa_raw, 5},
    {"_misskappa_rcpp_kappa_continuous", (DL_FUNC)&_misskappa_rcpp_kappa_continuous, 3},
    {"_misskappa_rcpp_alpha_raw", (DL_FUNC)&_misskappa_rcpp_alpha_raw, 4},
    {"_misskappa_rcpp_kappa_counts", (DL_FUNC)&_misskappa_rcpp_kappa_counts, 3},
    {"_misskappa_rcpp_kappa_quadratic", (DL_FUNC)&_misskappa_rcpp_kappa_quadratic, 4},
    {"_misskappa_rcpp_kappa_quadratic_counts", (DL_FUNC)&_misskappa_rcpp_kappa_quadratic_counts, 3},
    {"_misskappa_rcpp_kappa_fiml_counts", (DL_FUNC)&_misskappa_rcpp_kappa_fiml_counts, 5},
    {"_misskappa_rcpp_fiml_louis_spectrum", (DL_FUNC)&_misskappa_rcpp_fiml_louis_spectrum, 4},
    {"_misskappa_rcpp_kappa_gwise_categorical", (DL_FUNC)&_misskappa_rcpp_kappa_gwise_categorical, 4},
    {"_misskappa_rcpp_kappa_gwise_continuous", (DL_FUNC)&_misskappa_rcpp_kappa_gwise_continuous, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_misskappa(DllInfo* dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
