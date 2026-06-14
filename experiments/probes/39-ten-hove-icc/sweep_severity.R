# Probe 39 — benign -> severe gradient, reported as ABSOLUTE bias at a MEANINGFUL
# agreement level (shared component s2a=2 keeps truth in a normal range as the
# camp violation grows). delta = strength of the rater-camp structural violation.
#   struct = ten Hove structured (raters fixed); nt_fiml = saturated.
# Asymptotic n=30000 (MC SE ~0.004), MCAR 50% vs strong anchor-MAR.

suppressMessages({ library(misskappa); library(lme4) })
R <- 6; camp <- c(1,1,1,-1,-1,-1); mu <- 3; s2a <- 2; s2e <- 1; bj <- seq(-0.3,0.3,length.out=R)
gen <- function(n, delta, wtype) {
  a <- rnorm(n,0,sqrt(s2a)); w <- if (wtype=="gauss") rnorm(n) else sample(c(-1,1),n,replace=TRUE)
  mu + outer(a,rep(1,R)) + matrix(bj,n,R,byrow=TRUE) + outer(w, delta*camp) + matrix(rnorm(n*R,0,sqrt(s2e)),n,R)
}
conger_q <- function(m,S){R<-length(m);nu<-de<-0;for(i in 1:(R-1))for(j in (i+1):R){d2<-(m[i]-m[j])^2
  nu<-nu+S[i,i]+S[j,j]-2*S[i,j]+d2;de<-de+S[i,i]+S[j,j]+d2};1-nu/de}
truth_of <- function(delta){Sg<-s2a+delta^2*outer(camp,camp);diag(Sg)<-diag(Sg)+s2e;conger_q(mu+bj,Sg)}
ampute <- function(X, mech, rate=0.5){n<-nrow(X)
  if(mech=="mcar"){X[matrix(runif(n*R)<rate,n,R)]<-NA;return(X)}
  pr<-plogis(log(rate/(1-rate))+6*scale(X[,1])[,1]); for(j in 2:R) X[runif(n)<pr,j]<-NA; X}
est_nt <- function(X){f<-tryCatch(kappa(X,estimator="nt_fiml",weight="quadratic"),error=function(e)NULL);if(is.null(f))NA else unname(coef(f)["Conger"])}
est_st <- function(X){d<-data.frame(subj=factor(c(row(X))),rat=factor(c(col(X)),levels=seq_len(R)),y=c(X));d<-d[!is.na(d$y),]
  m<-tryCatch(suppressWarnings(lmer(y~rat+(1|subj),d,REML=FALSE,control=lmerControl(calc.derivs=FALSE))),error=function(e)NULL)
  if(is.null(m))return(NA);fe<-fixef(m);mr<-numeric(R);mr[1]<-fe[["(Intercept)"]];for(j in 2:R)mr[j]<-fe[["(Intercept)"]]+fe[[paste0("rat",j)]]
  s2s<-as.numeric(VarCorr(m)$subj[1]);s2sr<-sigma(m)^2;Sg<-matrix(s2s,R,R);diag(Sg)<-s2s+s2sr;conger_q(mr,Sg)}

rows <- list()
for (cfg in list(c("gauss","0"),c("gauss","1"),c("gauss","1.5"),c("gauss","2"),c("binary","1.5"))) {
  wt<-cfg[1]; delta<-as.numeric(cfg[2]); tru<-truth_of(delta)
  for (mc in c("mcar","mar_strong")) {
    set.seed(11); X<-ampute(gen(30000,delta,wt),mc); nt<-est_nt(X); st<-est_st(X)
    rows[[length(rows)+1]] <- data.frame(
      violation=sprintf("%s d=%.1f",wt,delta), truth=round(tru,3), missing=mc,
      nt_abs=round(nt-tru,3), struct_abs=round(st-tru,3),
      nt_kappa=round(nt,3), struct_kappa=round(st,3))
  }
}
cat("\nABSOLUTE bias for quadratic Conger kappa (n=30000, MC SE ~0.004)\n",
    "s2a=2 shared component keeps agreement meaningful; delta = camp violation\n\n",sep="")
print(do.call(rbind, rows), row.names=FALSE)
