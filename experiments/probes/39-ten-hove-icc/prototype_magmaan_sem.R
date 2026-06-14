# Probe 39 — PROVE the SEM route: the ten Hove structured estimator as a
# constrained one-factor CFA in magmaan, quadratic Conger kappa as a := defined
# parameter. Check it reproduces the probe's lmer-based `struct` estimate, on the
# same data, and yields a robust SE (no crossed-clustering / merDeriv headache:
# subjects are i.i.d. rows in the wide form).

suppressMessages({ library(lme4) })
core <- magmaan::magmaan_core

## camps DGP (Gaussian block covariance) + truth -----------------------------
R <- 6; camp <- c(1,1,1,-1,-1,-1); mu <- 3; s2a <- 2; s2e <- 1; delta <- 1.2; bj <- seq(-0.3,0.3,length.out=R)
gen <- function(n){a<-rnorm(n,0,sqrt(s2a));w<-rnorm(n)
  mu+outer(a,rep(1,R))+matrix(bj,n,R,byrow=TRUE)+outer(w,delta*camp)+matrix(rnorm(n*R,0,sqrt(s2e)),n,R)}
conger_q<-function(m,S){R<-length(m);nu<-de<-0;for(i in 1:(R-1))for(j in (i+1):R){d2<-(m[i]-m[j])^2
  nu<-nu+S[i,i]+S[j,j]-2*S[i,j]+d2;de<-de+S[i,i]+S[j,j]+d2};1-nu/de}
Sig<-s2a+delta^2*outer(camp,camp);diag(Sig)<-diag(Sig)+s2e; truth<-conger_q(mu+bj,Sig)

set.seed(42); X <- gen(400); X[matrix(runif(length(X))<0.4,nrow(X),ncol(X))] <- NA   # MCAR 40%
X <- X[rowSums(!is.na(X)) > 0, ]                       # drop all-NA subjects (uninformative)
colnames(X) <- paste0("r",1:R); df <- as.data.frame(X)

## (a) probe's lmer-based struct estimate ------------------------------------
est_st_lmer <- function(X){R<-ncol(X);d<-data.frame(subj=factor(c(row(X))),rat=factor(c(col(X)),levels=1:R),y=c(X));d<-d[!is.na(d$y),]
  m<-suppressWarnings(lmer(y~rat+(1|subj),d,REML=FALSE,control=lmerControl(calc.derivs=FALSE)))
  fe<-fixef(m);mr<-numeric(R);mr[1]<-fe[["(Intercept)"]];for(j in 2:R)mr[j]<-fe[["(Intercept)"]]+fe[[paste0("rat",j)]]
  s2s<-as.numeric(VarCorr(m)$subj[1]);s2sr<-sigma(m)^2;Sg<-matrix(s2s,R,R);diag(Sg)<-s2s+s2sr;conger_q(mr,Sg)}
k_lmer <- est_st_lmer(X)

## (b) magmaan one-factor SEM with kappa := defined parameter -----------------
build_model <- function(R) {
  load_terms <- paste(sprintf("1*r%d", 1:R), collapse=" + ")
  resid <- paste(sprintf("r%d ~~ e*r%d", 1:R, 1:R), collapse="\n")
  inter <- paste(sprintf("r%d ~ m%d*1", 1:R, 1:R), collapse="\n")
  pairs <- combn(R, 2)
  sdd <- paste(sprintf("(m%d-m%d)^2", pairs[1,], pairs[2,]), collapse=" + ")
  RR <- R*(R-1)
  paste0("f =~ ", load_terms, "\n",
         "f ~~ psi*f\n", resid, "\n", inter, "\n",
         sprintf("kappa := (%d*psi) / (%d*(psi+e) + (%s))", RR, RR, sdd))
}
model <- build_model(R)
fit <- magmaan::magmaan(model, df, estimator = "FIML", meanstructure = TRUE, se = "none", test = "none")

## naive (expected-info) delta-method SE
info_naive <- core$inference_information_expected(fit)
vc_naive   <- core$inference_vcov_partable(info_naive, fit$partable)
defs_naive <- magmaan::compute_defined(model, fit, vc_naive)

## robust FIML MLR sandwich -> delta-method SE on the kappa defined parameter
mlr <- core$estimate_fiml_robust_mlr(fit)
defs_robust <- tryCatch(magmaan::compute_defined(model, fit, mlr$vcov), error = function(e) {message(conditionMessage(e)); NULL})

k_sem  <- defs_naive$est[defs_naive$lhs=="kappa"]
cat(sprintf("truth quadratic Conger kappa = %.4f\n\n", truth))
cat(sprintf("struct via lmer (probe)      = %.5f\n", k_lmer))
cat(sprintf("struct via magmaan SEM       = %.5f   (match diff %.2e)\n", k_sem, abs(k_sem-k_lmer)))
cat(sprintf("  naive  delta-method SE     = %.4f\n", defs_naive$se[defs_naive$lhs=="kappa"]))
if (!is.null(defs_robust)) cat(sprintf("  robust FIML-MLR sandwich SE= %.4f\n", defs_robust$se[defs_robust$lhs=="kappa"]))
