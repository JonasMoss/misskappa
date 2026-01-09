#include "misskappa.h"
#include <vector>

namespace { // Anonymous namespace for internal implementation details

// --- Helper structs ported directly from quadagree.h ---
struct AcovResults {
  arma::mat psi;
  arma::vec mu_hat;
  arma::mat sigma_hat;
  int r;
};
  struct RawKappaEstimates {
    double fleiss_est;
    double conger_est;
    double bp_est;
    arma::mat scaled_acov;
    int n_eff;
  };


  // --- Core Calculation Engines ported directly from quadagree.h detail namespace ---

  // This is the "heavy-lifter" function that calculates the VCV of summary statistics
  // for data that may be missing.
  misskappa::Result<AcovResults> CalculatePsiInternal(const arma::mat& x, const arma::mat& M, int n_eff) {
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
    return {emdiscrete::Status::kOk, AcovResults{Psi, mu_hat, sigma_hat, r}};
  }

  // This is the orchestrator function that calls CalculatePsiInternal and computes the final estimates.
  misskappa::Result<RawKappaEstimates> CalculateRawEstimatesInternal(const arma::mat& x, double c1) {
    arma::mat M(x.n_rows, x.n_cols, arma::fill::zeros); M.elem(arma::find_finite(x)).ones();
    arma::uvec eff_rows = arma::find(arma::sum(M, 1) > 1); // Need at least 2 ratings for variance
    if (eff_rows.n_elem < 2) return {emdiscrete::Status::kError, std::nullopt, "Not enough subjects with at least two ratings."};
    arma::mat x_eff = x.rows(eff_rows), M_eff = M.rows(eff_rows); int n_eff = x_eff.n_rows;

    misskappa::Result<AcovResults> out_res = CalculatePsiInternal(x_eff, M_eff, n_eff);
    if (!out_res.IsOk()) return {emdiscrete::Status::kError, std::nullopt, "Internal Psi calculation failed."};

    const AcovResults& out = out_res.value.value();
    int r = out.r;
    double t1=arma::accu(out.sigma_hat), t2=arma::trace(out.sigma_hat), t3=arma::accu(arma::pow(out.mu_hat-arma::mean(out.mu_hat),2));
    double NF=t1-t2-t3, DF=(r-1.0)*(t2+t3), NC=t1-t2, DC=(r-1.0)*t2+(double)r*t3;
    double fleiss_est=(DF!=0)?NF/DF:arma::datum::nan, conger_est=(DC!=0)?NC/DC:arma::datum::nan;
    double d_obs=(2.0/(r-1.0))*(t2+t3-(1.0/r)*t1);
    double bp_est = (std::abs(c1) > 1e-9) ? 1.0 - d_obs/c1 : arma::datum::nan;

    arma::vec gF(3,arma::fill::zeros), gC(3,arma::fill::zeros);
    if(DF!=0) {gF(0)=1.0/DF; gF(1)=(-DF-NF*(r-1.0))/(DF*DF); gF(2)=(-DF-NF*(r-1.0))/(DF*DF);}
    if(DC!=0) {gC(0)=1.0/DC; gC(1)=(-DC-NC*(r-1.0))/(DC*DC); gC(2)=(-NC*(double)r)/(DC*DC);}

    arma::vec gBP(3,arma::fill::zeros);
    if(std::abs(c1) > 1e-9){double m=-(2.0/(c1*(r-1.0)));gBP={m*(-1.0/r),m,m};}

    arma::mat G(3,3); G.col(0)=gF; G.col(1)=gC; G.col(2)=gBP;

    // The scaled_acov is for [Fleiss, Conger, BP]
    return {emdiscrete::Status::kOk, RawKappaEstimates{fleiss_est, conger_est, bp_est, G.t()*out.psi*G, n_eff}};
  }

} // end anonymous namespace


namespace misskappa {
namespace kappaqp {

Result<Estimation> kappa(const arma::mat& x, const arma::vec& values) {
  int C = values.n_elem;
  double c1 = (std::abs(C) > 1e-9) ? (2.0 / (C * C)) * (C * arma::sum(arma::pow(values, 2)) - std::pow(arma::sum(values), 2)) : 0.0;

  // Call the new, general-purpose raw estimator
  auto estimates_res = CalculateRawEstimatesInternal(x, c1);
  if (!estimates_res.IsOk()) {
    return {emdiscrete::Status::kError, std::nullopt, estimates_res.error_message};
  }

  const auto& est_val = estimates_res.value.value();

  // The estimates in the returned struct are ordered [Fleiss, Conger, BP].
  // The misskappa raw output standard is [Conger, Fleiss, BP].
  arma::vec final_estimates = {est_val.conger_est, est_val.fleiss_est, est_val.bp_est};

  // The scaled_acov from the internal function is for [Fleiss, Conger, BP].
  // We must reorder it to match our final_estimates vector.
  arma::uvec reorder_indices = {1, 0, 2};
  arma::mat vcov_reordered = est_val.scaled_acov.submat(reorder_indices, reorder_indices);

  // Unscale the variance-covariance matrix by the number of effective subjects.
  arma::mat final_vcov = vcov_reordered / est_val.n_eff;

  return {emdiscrete::Status::kOk, Estimation{final_estimates, final_vcov}, ""};
}

Result<Estimation> kappa_counts(const arma::mat& x, const arma::vec& values, int R) {
  if (x.n_rows < 2) return {emdiscrete::Status::kError, std::nullopt, "Input matrix must have at least 2 rows."};

  // Compute c1 internally
  int C = values.n_elem;
  double c1 = (2.0 / (C * C)) * (C * arma::sum(arma::pow(values, 2)) - std::pow(arma::sum(values), 2));

  arma::uvec eff_rows=arma::find(arma::sum(x,1)>0);
  if(eff_rows.n_elem<2) return {emdiscrete::Status::kError, std::nullopt, "Not enough subjects with ratings."};
  arma::mat x_eff=x.rows(eff_rows); int n_eff=x_eff.n_rows;
  arma::vec r_i_eff=arma::sum(x_eff,1);
  arma::vec s_i=((double)R/r_i_eff)%(x_eff*values), q_i_sq=((double)R/r_i_eff)%(x_eff*arma::pow(values,2));
  arma::mat Z=arma::join_rows(s_i, arma::pow(s_i,2), q_i_sq);
  arma::vec phi=arma::mean(Z,0).t(); arma::mat Psi_dist=arma::cov(Z,1);

  // Fleiss' Kappa Calculations
  double est_f=arma::datum::nan;
  arma::vec grad_f(3, arma::fill::zeros);
  double Df=phi(2)-(phi(0)*phi(0)/R);
  if(Df>1e-9){
    double Nf=phi(1)-phi(0)*phi(0);
    double r_inv=1.0/(R-1.0);
    grad_f={r_inv*((-2.0*phi(0)/Df)+(Nf*(2.0*phi(0)/R))/(Df*Df)),r_inv/Df,r_inv*(-Nf/(Df*Df))};
    est_f=r_inv*(Nf/Df - 1.0);
  }

  // Brennan-Prediger Calculations
  double est_bp=arma::datum::nan;
  arma::vec grad_bp(3, arma::fill::zeros);
  if(std::abs(c1) > 1e-9){
    double d_bp=(2.0/(R-1.0))*(phi(2)-(1.0/R)*phi(1));
    est_bp=1.0-d_bp/c1;
    double m=-(2.0/(c1*(R-1.0)));
    grad_bp={0.0, m*(-1.0/R), m};
  }

  // Assemble gradient matrix and compute full V-Cov matrix
  arma::mat G(3,2);
  G.col(0) = grad_f;
  G.col(1) = grad_bp;

  arma::vec estimates = {est_f, est_bp};
  arma::mat vcov = (G.t() * Psi_dist * G) / n_eff;

  return {emdiscrete::Status::kOk, Estimation{estimates, vcov}, ""};
}

} // namespace kappaqp
} // namespace misskappa
