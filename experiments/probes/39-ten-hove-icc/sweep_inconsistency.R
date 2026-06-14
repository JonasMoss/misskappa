# Probe 39 — MAP the inconsistency: when is the ten Hove structured estimator
# (and the saturated nt_fiml) biased for the quadratic Conger kappa, from benign
# to severe? Grid = structural violation x missingness mechanism. Asymptotic
# (n=20000) %bias, so any nonzero entry is genuine inconsistency, not finite-n.
#
# Rows (structural violation of the additive/compound-symmetry model):
#   compound  : delta=0, truly additive Gaussian          (struct CORRECTLY specified)
#   camps_G   : Gaussian latent attribute, rater camps +/- (block cov; nt_fiml still correct)
#   group_NG  : BINARY group attribute, rater camps +/-    (2-component mixture; nt_fiml ALSO wrong)
# Cols (missingness):
#   mcar      : 50% completely at random
#   mar_mild  : anchor-MAR on rater 1, slope 2
#   mar_strong: anchor-MAR on rater 1, slope 6

suppressMessages({ library(misskappa); library(lme4) })
R <- 6; camp <- c(1,1,1,-1,-1,-1); mu <- 3; s2a <- 1; s2e <- 1; bj <- seq(-0.3,0.3,length.out=R)

gen <- function(n, delta, wtype) {
  a <- rnorm(n,0,sqrt(s2a))
  w <- if (wtype=="gauss") rnorm(n) else sample(c(-1,1), n, replace=TRUE)   # binary => mixture
  mu + outer(a,rep(1,R)) + matrix(bj,n,R,byrow=TRUE) + outer(w, delta*camp) +
    matrix(rnorm(n*R,0,sqrt(s2e)),n,R)
}
conger_q <- function(m,S){R<-length(m);nu<-de<-0;for(i in 1:(R-1))for(j in (i+1):R){d2<-(m[i]-m[j])^2
  nu<-nu+S[i,i]+S[j,j]-2*S[i,j]+d2;de<-de+S[i,i]+S[j,j]+d2};1-nu/de}
truth_of <- function(delta){                         # Var(w)=1 either way (gauss or +/-1)
  Sg <- s2a + delta^2*outer(camp,camp); diag(Sg)<-diag(Sg)+s2e; conger_q(mu+bj, Sg) }

ampute <- function(X, mech, rate=0.5) {
  n<-nrow(X)
  if (mech=="mcar") { X[matrix(runif(n*R)<rate,n,R)] <- NA; return(X) }
  slope <- if (mech=="mar_mild") 2 else 6
  pr <- plogis(log(rate/(1-rate)) + slope*scale(X[,1])[,1])     # raters 2..R drop, MAR on rater 1
  for (j in 2:R) X[runif(n)<pr, j] <- NA
  X
}
est_nt <- function(X){f<-tryCatch(kappa(X,estimator="nt_fiml",weight="quadratic"),error=function(e)NULL)
  if(is.null(f))NA else unname(coef(f)["Conger"])}
est_st <- function(X){d<-data.frame(subj=factor(c(row(X))),rat=factor(c(col(X)),levels=seq_len(R)),y=c(X));d<-d[!is.na(d$y),]
  m<-tryCatch(suppressWarnings(lmer(y~rat+(1|subj),d,REML=FALSE,control=lmerControl(calc.derivs=FALSE))),error=function(e)NULL)
  if(is.null(m))return(NA);fe<-fixef(m);mr<-numeric(R);mr[1]<-fe[["(Intercept)"]]
  for(j in 2:R)mr[j]<-fe[["(Intercept)"]]+fe[[paste0("rat",j)]]
  s2s<-as.numeric(VarCorr(m)$subj[1]);s2sr<-sigma(m)^2;Sg<-matrix(s2s,R,R);diag(Sg)<-s2s+s2sr;conger_q(mr,Sg)}

viol <- list(compound=list(delta=0,wt="gauss"), camps_G=list(delta=2,wt="gauss"), group_NG=list(delta=2,wt="binary"))
mechs <- c("mcar","mar_mild","mar_strong")
rows <- list()
for (vn in names(viol)) {
  v <- viol[[vn]]; tru <- truth_of(v$delta)
  for (mc in mechs) {
    set.seed(7); X <- ampute(gen(20000, v$delta, v$wt), mc)
    nt <- est_nt(X); st <- est_st(X)
    rows[[length(rows)+1]] <- data.frame(violation=vn, truth=round(tru,3), missing=mc,
      nt_fiml_pct=round(100*(nt-tru)/tru,1), struct_pct=round(100*(st-tru)/tru,1))
  }
}
cat("\nAsymptotic %bias for quadratic Conger kappa (n=20000)\n",
    "rows: structural violation; cols within: missingness\n\n", sep="")
print(do.call(rbind, rows), row.names=FALSE)
