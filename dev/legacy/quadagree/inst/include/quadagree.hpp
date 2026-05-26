#ifndef QUADAGREE_QUADAGREE_H_
#define QUADAGREE_QUADAGREE_H_

#include <cmath>
#include <cstddef>

#include <array>
#include <limits>
#include <optional>
#include <random>
#include <string>
#include <vector>

#include <armadillo>

#if __cplusplus >= 201703L
// standard optional is available
#else
#include <experimental/optional>
namespace std {
using experimental::nullopt;
using experimental::optional;
}  // namespace std
#endif

namespace quadagree {

// =============================================================================
// Public API Structs & Enums
// These are the types that the user-facing API functions will return and use.
// =============================================================================

/// @brief Status code for fallible operations.
enum class Status { kOk, kError };

/// @brief A wrapper for operations that can fail.
template <typename T>
struct Result {
  Status status;
  std::optional<T> value;
};

// --- Public Enums ---
enum class Alternative { kTwoSided, kGreater, kLess };
enum class Transform { kNone, kFisher, kLog, kArcsin };

// --- Public Result Structs ---
struct CIResults {
  double estimate;
  double std_err;
  double conf_low;
  double conf_high;
  int n_eff;
};

struct AllCIResults {
  CIResults fleiss;
  CIResults conger;
  CIResults bp;
  arma::mat scaled_acov;
};

struct AggrCIResults {
  CIResults fleiss;
  CIResults bp;
};


// =============================================================================
// Internal Implementation Details
// Functions and types in this namespace are not part of the public API
// and are subject to change without notice.
// =============================================================================
namespace detail {

// --- Internal Result Structs ---
struct AcovResults { arma::mat psi; arma::vec mu_hat; arma::mat sigma_hat; int r; };
struct RawKappaEstimates { double fleiss_est; double conger_est; double bp_est; arma::mat scaled_acov; int n_eff; };
struct AggrKappaEstimates { double fleiss_est; double bp_est; double fleiss_var; double bp_var; int n_eff; };

// --- Transformer Functors ---
struct NoneTransform { static double Est(double v){return v;} static double Sd(double e, double s){return s;} static double Inv(double v){return v;} };
struct FisherTransform { static double Est(double v){return std::atanh(v);} static double Sd(double e, double s){return s/(1.0-e*e);} static double Inv(double v){return std::tanh(v);} };
struct LogTransform { static double Est(double v){return std::log(1.0-v);} static double Sd(double e, double s){return s/std::abs(1.0-e);} static double Inv(double v){return 1.0-std::exp(v);} };
struct ArcsinTransform { static double Est(double v){return std::asin(v);} static double Sd(double e, double s){return s/std::sqrt(1.0-e*e);} static double Inv(double v){return std::sin(v);} };

/// @brief Compute the inverse CDF (quantile) of the standard normal distribution.
inline double QNorm(double p) {
  using std::log;
  using std::min;
  using std::numeric_limits;
  using std::sqrt;

  if (p <= 0.0) {
    if (p == 0.0) return -numeric_limits<double>::infinity();
    return numeric_limits<double>::quiet_NaN();
  }
  if (p >= 1.0) {
    if (p == 1.0) return numeric_limits<double>::infinity();
    return numeric_limits<double>::quiet_NaN();
  }

  // Coefficients
  static constexpr double a1 = -3.969683028665376e+01;
  static constexpr double a2 = 2.209460984245205e+02;
  static constexpr double a3 = -2.759285104469687e+02;
  static constexpr double a4 = 1.383577518672690e+02;
  static constexpr double a5 = -3.066479806614716e+01;
  static constexpr double a6 = 2.506628277459239e+00;

  static constexpr double b1 = -5.447609879822406e+01;
  static constexpr double b2 = 1.615858368580409e+02;
  static constexpr double b3 = -1.556989798598866e+02;
  static constexpr double b4 = 6.680131188771972e+01;
  static constexpr double b5 = -1.328068155288572e+01;

  static constexpr double c1 = -7.784894002430293e-03;
  static constexpr double c2 = -3.223964580411365e-01;
  static constexpr double c3 = -2.400758277161838e+00;
  static constexpr double c4 = -2.549732539343734e+00;
  static constexpr double c5 = 4.374664141464968e+00;
  static constexpr double c6 = 2.938163982698783e+00;

  static constexpr double d1 = 7.784695709041462e-03;
  static constexpr double d2 = 3.224671290700398e-01;
  static constexpr double d3 = 2.445134137142996e+00;
  static constexpr double d4 = 3.754408661907416e+00;

  const double plow = 0.02425;

  // Use symmetry: work with the smaller tail
  double q = min(p, 1.0 - p);
  double u;

  if (q > plow) {
    // Central region
    double d = q - 0.5;
    double t = d * d;
    u = d * (((((a1 * t + a2) * t + a3) * t + a4) * t + a5) * t + a6) /
      (((((b1 * t + b2) * t + b3) * t + b4) * t + b5) * t + 1.0);
  } else {
    // Tail region
    double t = sqrt(-2.0 * log(q));
    u = (((((c1 * t + c2) * t + c3) * t + c4) * t + c5) * t + c6) /
      (((((d1 * t + d2) * t + d3) * t + d4) * t + 1.0));
  }

  // Restore the sign for the upper tail
  return (p > 0.5) ? -u : u;
}

inline arma::vec GetCiQuantiles(Alternative alt, double conf_level) {
  double half_alpha = (1.0 - conf_level) / 2.0;
  if (alt == Alternative::kGreater) return {1.0 - conf_level, 1.0};
  if (alt == Alternative::kLess) return {0.0, conf_level};
  return {half_alpha, 1.0 - half_alpha};
}

// --- Asymptotic and Bootstrap CI Calculations ---
template <typename T_Transformer>
inline arma::vec CiAsymptotic(double est, double sd, const arma::vec& quants) {
  double est_t = T_Transformer::Est(est);
  double sd_t = T_Transformer::Sd(est, sd);
  double z_lower = QNorm(quants(0));
  double z_upper = QNorm(quants(1));
  arma::vec ci = {T_Transformer::Inv(est_t + z_lower * sd_t), T_Transformer::Inv(est_t + z_upper * sd_t)};
  return arma::sort(ci);
}

// --- Core Calculation Engines ---
Result<AcovResults> CalculatePsiInternal(const arma::mat& x, const arma::mat& M, int n_eff);
Result<RawKappaEstimates> CalculateCompleteEstimatesInternal(const arma::mat& x, double c1);
Result<RawKappaEstimates> CalculateRawEstimatesInternal(const arma::mat& x, double c1);
Result<AggrKappaEstimates> CalculateAggrEstimatesInternal(const arma::mat& x, const arma::vec& values, int R, double c1);

inline Result<RawKappaEstimates> CalculateCompleteEstimatesInternal(const arma::mat& x, double c1) {
  const int n=x.n_rows, r=x.n_cols;
  const arma::vec mu_hat=arma::mean(x,0).t();
  const arma::mat sigma_hat=arma::cov(x,1);
  const arma::mat Y_hat=x.each_row()-mu_hat.t();
  arma::mat gamma_ss(r*(r+1)/2,r*(r+1)/2,arma::fill::zeros);
  int row_idx=0; for(int j=0; j<r; ++j) for(int i=j; i<r; ++i) {
    int col_idx=0; for(int l=0; l<r; ++l) for(int k=l; k<r; ++k) {
      if(col_idx>=row_idx) {
        double w_ij = arma::as_scalar(arma::mean(Y_hat.col(i)%Y_hat.col(j)));
        double w_kl = arma::as_scalar(arma::mean(Y_hat.col(k)%Y_hat.col(l)));
        double w_ijkl = arma::as_scalar(arma::mean(Y_hat.col(i)%Y_hat.col(j)%Y_hat.col(k)%Y_hat.col(l)));
        gamma_ss(row_idx, col_idx) = w_ijkl - w_ij * w_kl;
      } col_idx++;
    } row_idx++;
  }
  gamma_ss = arma::symmatu(gamma_ss);
  arma::mat gamma_s_mu(r*(r+1)/2,r,arma::fill::zeros);
  row_idx=0; for(int j=0; j<r; ++j) for(int i=j; i<r; ++i) {
    for(int k=0; k<r; ++k) { gamma_s_mu(row_idx,k)=arma::as_scalar(arma::mean(Y_hat.col(i)%Y_hat.col(j)%Y_hat.col(k))); }
    row_idx++;
  }
  arma::mat top=arma::join_rows(gamma_ss,gamma_s_mu), bot=arma::join_rows(gamma_s_mu.t(),sigma_hat);
  arma::mat gamma_full=arma::join_cols(top,bot);
  arma::mat J(3,gamma_full.n_cols,arma::fill::zeros);
  int p=r*(r+1)/2; J(0,arma::span(0,p-1)).ones();
  row_idx=0; int current_diag=0; for(int j=0; j<r; ++j) {J(1,row_idx+current_diag)=1.0; row_idx+=(r-j); current_diag++;}
  J(2,arma::span(p,p+r-1)) = 2.0*(mu_hat-arma::mean(mu_hat)).t();
  arma::mat Psi = J*gamma_full*J.t();
  double t1=arma::accu(sigma_hat), t2=arma::trace(sigma_hat), t3=arma::accu(arma::pow(mu_hat-arma::mean(mu_hat),2));
  double NF=t1-t2-t3, DF=(r-1.0)*(t2+t3), NC=t1-t2, DC=(r-1.0)*t2+(double)r*t3;
  double fleiss_est=(DF!=0)?NF/DF:arma::datum::nan, conger_est=(DC!=0)?NC/DC:arma::datum::nan;
  double d_obs=(2.0/(r-1.0))*(t2+t3-(1.0/r)*t1); double bp_est=1.0-d_obs/c1;
  arma::vec gF(3,arma::fill::zeros), gC(3,arma::fill::zeros);
  if(DF!=0){gF(0)=1.0/DF; gF(1)=(-DF-NF*(r-1.0))/(DF*DF); gF(2)=(-DF-NF*(r-1.0))/(DF*DF);}
  if(DC!=0){gC(0)=1.0/DC; gC(1)=(-DC-NC*(r-1.0))/(DC*DC); gC(2)=(-NC*(double)r)/(DC*DC);}
  arma::vec gBP(3,arma::fill::zeros); if(c1!=0){double m=-(2.0/(c1*(r-1.0)));gBP={m*(-1.0/r),m,m};}
  arma::mat G(3,3); G.col(0)=gF; G.col(1)=gC; G.col(2)=gBP;
  return {Status::kOk, RawKappaEstimates{fleiss_est, conger_est, bp_est, G.t()*Psi*G, n}};
}

inline Result<RawKappaEstimates> CalculateRawEstimatesInternal(const arma::mat& x, double c1) {
  //if (!x.has_nan()) return CalculateCompleteEstimatesInternal(x, c1);
  arma::mat M(x.n_rows, x.n_cols, arma::fill::zeros); M.elem(arma::find_finite(x)).ones();
  arma::uvec eff_rows = arma::find(arma::sum(M, 1) > 0);
  if (eff_rows.n_elem < 2) return {Status::kError, std::nullopt};
  arma::mat x_eff = x.rows(eff_rows), M_eff = M.rows(eff_rows); int n_eff = x_eff.n_rows;
  Result<AcovResults> out_res = CalculatePsiInternal(x_eff, M_eff, n_eff);
  if (out_res.status != Status::kOk) return {Status::kError, std::nullopt};
  const AcovResults& out = *out_res.value;
  int r = out.r;
  double t1=arma::accu(out.sigma_hat), t2=arma::trace(out.sigma_hat), t3=arma::accu(arma::pow(out.mu_hat-arma::mean(out.mu_hat),2));
  double NF=t1-t2-t3, DF=(r-1.0)*(t2+t3), NC=t1-t2, DC=(r-1.0)*t2+(double)r*t3;
  double fleiss_est=(DF!=0)?NF/DF:arma::datum::nan, conger_est=(DC!=0)?NC/DC:arma::datum::nan;
  double d_obs=(2.0/(r-1.0))*(t2+t3-(1.0/r)*t1); double bp_est=1.0-d_obs/c1;
  arma::vec gF(3,arma::fill::zeros), gC(3,arma::fill::zeros);
  if(DF!=0) {gF(0)=1.0/DF; gF(1)=(-DF-NF*(r-1.0))/(DF*DF); gF(2)=(-DF-NF*(r-1.0))/(DF*DF);}
  if(DC!=0) {gC(0)=1.0/DC; gC(1)=(-DC-NC*(r-1.0))/(DC*DC); gC(2)=(-NC*(double)r)/(DC*DC);}
  arma::vec gBP(3,arma::fill::zeros); if(c1!=0){double m=-(2.0/(c1*(r-1.0)));gBP={m*(-1.0/r),m,m};}
  arma::mat G(3,3); G.col(0)=gF; G.col(1)=gC; G.col(2)=gBP;
  return {Status::kOk, RawKappaEstimates{fleiss_est, conger_est, bp_est, G.t()*out.psi*G, n_eff}};
}

inline Result<AggrKappaEstimates> CalculateAggrEstimatesInternal(const arma::mat& x, const arma::vec& values, int R, double c1) {
  arma::uvec eff_rows=arma::find(arma::sum(x,1)>0);
  if(eff_rows.n_elem<2) return {Status::kError, std::nullopt};
  arma::mat x_eff=x.rows(eff_rows); int n_eff=x_eff.n_rows;
  arma::vec r_i_eff=arma::sum(x_eff,1);
  arma::vec s_i=((double)R/r_i_eff)%(x_eff*values), q_i_sq=((double)R/r_i_eff)%(x_eff*arma::pow(values,2));
  arma::mat Z=arma::join_rows(s_i, arma::pow(s_i,2), q_i_sq);
  arma::vec phi=arma::mean(Z,0).t(); arma::mat Psi_dist=arma::cov(Z,1);
  double Df=phi(2)-(phi(0)*phi(0)/R), est_f=arma::datum::nan, var_f=0.0;
  if(Df>1e-9){
    double Nf=phi(1)-phi(0)*phi(0), r_inv=1.0/(R-1.0);
    arma::vec grad_f={r_inv*((-2.0*phi(0)/Df)+(Nf*(2.0*phi(0)/R))/(Df*Df)),r_inv/Df,r_inv*(-Nf/(Df*Df))};
    var_f=arma::as_scalar(grad_f.t()*Psi_dist*grad_f);
    est_f=r_inv*((phi(1)-phi(0)*phi(0))/Df - 1.0);
  }
  double d_bp=(2.0/(R-1.0))*(phi(2)-(1.0/R)*phi(1)), est_bp=1.0-d_bp/c1, var_bp=0.0;
  if(c1!=0){
    double m=-(2.0/(c1*(R-1.0)));
    arma::vec grad_bp={0.0, m*(-1.0/R), m};
    var_bp=arma::as_scalar(grad_bp.t()*Psi_dist*grad_bp);
  }
  return {Status::kOk, AggrKappaEstimates{est_f, est_bp, var_f, var_bp, n_eff}};
}

inline Result<AcovResults> CalculatePsiInternal(const arma::mat& x, const arma::mat& M, int n_eff) {
  int r=x.n_cols; arma::mat Xf=x; Xf.elem(arma::find_nonfinite(Xf)).zeros();
  arma::vec count=arma::sum(M,0).t(), mu_hat(r);
  for(int j=0; j<r; ++j){mu_hat(j)=(count(j)>0)?arma::accu(Xf.col(j))/double(count(j)):arma::datum::nan;}
  arma::mat Y_hat=x; Y_hat.each_row()-=mu_hat.t(); Y_hat.elem(arma::find_nonfinite(x)).zeros();
  arma::vec p1=count/double(n_eff); arma::mat p2=M.t()*M/double(n_eff), sigma_hat(r,r,arma::fill::zeros);
  for(int j=0; j<r; ++j) {for(int k=j; k<r; ++k) {
    arma::uvec idx=arma::find(M.col(j)==1&&M.col(k)==1);
    if(idx.n_elem>0) {
      double sum_j=0.0,sum_k=0.0; for(arma::uword i=0;i<idx.n_elem;++i){arma::uword r_idx=idx(i); sum_j+=x(r_idx,j); sum_k+=x(r_idx,k);}
      double mu_j_p=sum_j/idx.n_elem, mu_k_p=sum_k/idx.n_elem, sum_prod=0.0;
      for(arma::uword i=0;i<idx.n_elem;++i){arma::uword r_idx=idx(i); sum_prod+=(x(r_idx,j)-mu_j_p)*(x(r_idx,k)-mu_k_p);}
      sigma_hat(j,k)=sum_prod/idx.n_elem;
    }
    sigma_hat(k,j)=sigma_hat(j,k);
  }}
  arma::cube mu3(r,r,r,arma::fill::zeros), p3(r,r,r,arma::fill::zeros);
  arma::field<arma::mat> mu4(r,r), p4(r,r);
  for(int i=0;i<r;++i)for(int j=0;j<r;++j){mu4(i,j).zeros(r,r); p4(i,j).zeros(r,r);}
  for(int i=0;i<r;++i)for(int j=0;j<r;++j)for(int k=0;k<r;++k) {
    arma::uvec idx3=arma::find(M.col(i)%M.col(j)%M.col(k));
    p3(i,j,k)=double(idx3.n_elem)/double(n_eff);
    if(idx3.n_elem>0){arma::vec y_i=Y_hat.col(i),y_j=Y_hat.col(j),y_k=Y_hat.col(k); mu3(i,j,k)=arma::mean(y_i.elem(idx3)%y_j.elem(idx3)%y_k.elem(idx3));}
    for(int l=0;l<r;++l){
      arma::uvec idx4=arma::find(M.col(i)%M.col(j)%M.col(k)%M.col(l));
      p4(i,j)(k,l)=double(idx4.n_elem)/double(n_eff);
      if(idx4.n_elem>0){arma::vec y_i=Y_hat.col(i),y_j=Y_hat.col(j),y_k=Y_hat.col(k),y_l=Y_hat.col(l); mu4(i,j)(k,l)=arma::mean(y_i.elem(idx4)%y_j.elem(idx4)%y_k.elem(idx4)%y_l.elem(idx4));}
    }
  }
  arma::mat Psi(3,3,arma::fill::zeros); arma::vec v_dmu=2*(mu_hat-arma::mean(mu_hat));
  for(int i=0;i<r;++i)for(int j=0;j<r;++j)for(int k=0;k<r;++k)for(int l=0;l<r;++l){
    if(p2(i,j)>1e-9&&p2(k,l)>1e-9){
      double gamma=mu4(i,j)(k,l)-sigma_hat(i,j)*sigma_hat(k,l);
      double term=(p4(i,j)(k,l)/(p2(i,j)*p2(k,l)))*gamma;
      Psi(0,0)+=term; if(k==l)Psi(0,1)+=term; if(i==j&&k==l)Psi(1,1)+=term;
    }
  }
  Psi(1,0)=Psi(0,1);
  for(int i=0;i<r;++i)for(int j=0;j<r;++j){
    if(p1(i)>1e-9&&p1(j)>1e-9){Psi(2,2)+=v_dmu(i)*v_dmu(j)*(p2(i,j)/(p1(i)*p1(j)))*sigma_hat(i,j);}
  }
  for(int i=0;i<r;++i)for(int j=0;j<r;++j)for(int k=0;k<r;++k){
    if(p1(i)>1e-9&&p2(j,k)>1e-9){
      double omega=(p3(i,j,k)/(p1(i)*p2(j,k)))*mu3(i,j,k);
      Psi(0,2)+=v_dmu(i)*omega; if(j==k)Psi(1,2)+=v_dmu(i)*omega;
    }
  }
  Psi(2,0)=Psi(0,2); Psi(2,1)=Psi(1,2);
  return {Status::kOk, AcovResults{Psi, mu_hat, sigma_hat, r}};
}

} // namespace detail

// =============================================================================
// Public C++ API
// =============================================================================

/// @brief Orchestrates raw kappa calculations for Fleiss, Conger, and BP.
inline Result<AllCIResults> QuadagreeRaw(
    const arma::mat& x, double c1, Transform transform_type, double conf_level,
    Alternative alternative, bool bootstrap, int n_reps, std::mt19937& rng) {

  Result<detail::RawKappaEstimates> estimates_res = detail::CalculateRawEstimatesInternal(x, c1);
  if (estimates_res.status != Status::kOk) return {Status::kError, std::nullopt};

  const auto& estimates = *estimates_res.value;
  const double fleiss_sd = (estimates.scaled_acov(0,0)>0) ? std::sqrt(estimates.scaled_acov(0,0)/estimates.n_eff) : 0.0;
  const double conger_sd = (estimates.scaled_acov(1,1)>0) ? std::sqrt(estimates.scaled_acov(1,1)/estimates.n_eff) : 0.0;
  const double bp_sd = (estimates.scaled_acov(2,2)>0) ? std::sqrt(estimates.scaled_acov(2,2)/estimates.n_eff) : 0.0;

  const arma::vec quants = detail::GetCiQuantiles(alternative, conf_level);
  arma::vec fleiss_ci, conger_ci, bp_ci;

  auto calculate_cis = [&](auto transformer) -> Status {
    if (bootstrap) return Status::kError; // Bootstrap not yet implemented for 3-way CI
    fleiss_ci = detail::CiAsymptotic<decltype(transformer)>(estimates.fleiss_est, fleiss_sd, quants);
    conger_ci = detail::CiAsymptotic<decltype(transformer)>(estimates.conger_est, conger_sd, quants);
    bp_ci = detail::CiAsymptotic<decltype(transformer)>(estimates.bp_est, bp_sd, quants);
    return Status::kOk;
  };

  Status op_status = Status::kError;
  switch (transform_type) {
  case Transform::kFisher: op_status = calculate_cis(detail::FisherTransform{}); break;
  case Transform::kLog: op_status = calculate_cis(detail::LogTransform{}); break;
  case Transform::kArcsin: op_status = calculate_cis(detail::ArcsinTransform{}); break;
  default: op_status = calculate_cis(detail::NoneTransform{}); break;
  }
  if (op_status != Status::kOk) return {Status::kError, std::nullopt};

  CIResults fleiss_res = {estimates.fleiss_est, fleiss_sd, fleiss_ci(0), std::min(fleiss_ci(1),1.0), estimates.n_eff};
  CIResults conger_res = {estimates.conger_est, conger_sd, conger_ci(0), std::min(conger_ci(1),1.0), estimates.n_eff};
  CIResults bp_res = {estimates.bp_est, bp_sd, bp_ci(0), std::min(bp_ci(1),1.0), estimates.n_eff};

  return {Status::kOk, AllCIResults{fleiss_res, conger_res, bp_res, estimates.scaled_acov}};
}

/// @brief Orchestrates aggregated kappa calculations for Fleiss and BP.
inline Result<AggrCIResults> QuadagreeAggr(
    const arma::mat& x, const arma::vec& values, int R, double c1,
    Transform transform_type, double conf_level, Alternative alternative,
    bool bootstrap, int n_reps, std::mt19937& rng) {

  Result<detail::AggrKappaEstimates> estimates_res = detail::CalculateAggrEstimatesInternal(x, values, R, c1);
  if (estimates_res.status != Status::kOk) return {Status::kError, std::nullopt};

  const auto& estimates = *estimates_res.value;
  const double fleiss_sd = (estimates.fleiss_var > 0) ? std::sqrt(estimates.fleiss_var / estimates.n_eff) : 0.0;
  const double bp_sd = (estimates.bp_var > 0) ? std::sqrt(estimates.bp_var / estimates.n_eff) : 0.0;

  const arma::vec quants = detail::GetCiQuantiles(alternative, conf_level);
  arma::vec fleiss_ci, bp_ci;

  auto calculate_cis = [&](auto transformer) -> Status {
    if (bootstrap) return Status::kError; // Not yet implemented
    fleiss_ci = detail::CiAsymptotic<decltype(transformer)>(estimates.fleiss_est, fleiss_sd, quants);
    bp_ci = detail::CiAsymptotic<decltype(transformer)>(estimates.bp_est, bp_sd, quants);
    return Status::kOk;
  };

  Status op_status = Status::kError;
  switch (transform_type) {
  case Transform::kFisher: op_status = calculate_cis(detail::FisherTransform{}); break;
  case Transform::kLog: op_status = calculate_cis(detail::LogTransform{}); break;
  case Transform::kArcsin: op_status = calculate_cis(detail::ArcsinTransform{}); break;
  default: op_status = calculate_cis(detail::NoneTransform{}); break;
  }
  if (op_status != Status::kOk) return {Status::kError, std::nullopt};

  CIResults fleiss_res = {estimates.fleiss_est, fleiss_sd, fleiss_ci(0), std::min(fleiss_ci(1),1.0), estimates.n_eff};
  CIResults bp_res = {estimates.bp_est, bp_sd, bp_ci(0), std::min(bp_ci(1),1.0), estimates.n_eff};

  return {Status::kOk, AggrCIResults{fleiss_res, bp_res}};
}


} // namespace quadagree
#endif // QUADAGREE_QUADAGREE_H_
