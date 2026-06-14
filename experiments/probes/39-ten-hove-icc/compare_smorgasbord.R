# Probe 39 — estimand-MATCHED comparison across a smorgasbord of canonical DGPs
# (reused from studies 23 & 29). BOTH estimators now target the SAME object:
# the population quadratic Conger kappa = conger_q(mu, Sigma).
#   nt_fiml : saturated normal (mu, Sigma) -> Conger      [misskappa]
#   struct  : ten Hove structured model: subjects random, raters FIXED
#             lmer(y ~ rat + (1|subj)) -> mu_r, Sigma = s2s*J + s2sr*I -> conger_q
# This removes the fixed-vs-random-rater R/(R-1) convention gap from the earlier
# run: any remaining difference is genuine estimation bias/efficiency.
#
# Reports, per DGP: truth; and per estimator the ASYMPTOTIC %bias (large-n,
# complete data -> consistency check) and the finite-sample %bias / SD / RMSE
# (S small, MCAR) -> efficiency.

suppressMessages({ library(misskappa); library(lme4) })

## ---- DGPs (verbatim from studies/23 and studies/29) ------------------------
make_Sigma <- function(R, rho = 0.6) { S <- matrix(rho,R,R); diag(S)<-1; outer(sqrt(seq(0.8,1.4,length.out=R)), sqrt(seq(0.8,1.4,length.out=R)))*S }
make_mu    <- function(R) seq(-0.45, 0.45, length.out = R)
gen_continuous <- function(n, R, dist) {
  Sigma <- make_Sigma(R); mu <- make_mu(R)
  Z <- if (dist=="normal") matrix(rnorm(n*R),n,R) else matrix((rchisq(n*R,df=3)-3)/sqrt(6),n,R)
  sweep(Z %*% chol(Sigma), 2L, mu, "+")
}
gen_mixture <- function(n, R, rho0=0.1, rho1=0.85, shift=1.0) {
  mu <- make_mu(R); d <- sqrt(seq(0.8,1.4,length.out=R)); z <- rbinom(n,1L,0.5); X <- matrix(0,n,R)
  for (cl in 0:1) { idx<-which(z==cl); if(!length(idx)) next
    rho<-if(cl==0)rho0 else rho1; sh<-if(cl==0)-shift else shift
    S<-matrix(rho,R,R); diag(S)<-1; S<-outer(d,d)*S
    X[idx,]<-sweep(matrix(rnorm(length(idx)*R),length(idx),R)%*%chol(S),2L,mu+sh,"+") }
  X
}
gen_categorical <- function(n, R) {           # study 23 skill-guess, C=5
  C<-5L; vals<-seq_len(C); s<-seq(0.65,0.78,length.out=R); centers<-seq(2,C-1,length.out=R)
  guess<-vapply(centers,function(ctr){p<-exp(-((vals-ctr)^2)/(2*1.2^2));p/sum(p)},numeric(C))
  Tcat<-sample(vals,n,replace=TRUE); X<-matrix(0,n,R)
  for (j in seq_len(R)){sk<-runif(n)<s[j]; g<-sample(vals,n,replace=TRUE,prob=guess[,j]); X[,j]<-ifelse(sk,Tcat,g)}
  X
}
simulate_copula <- function(n, R, corr, p) {
  C<-length(p); Sigma<-matrix(corr,R,R); diag(Sigma)<-1
  Z<-matrix(rnorm(n*R),n,R)%*%chol(Sigma); cuts<-c(-Inf,qnorm(cumsum(p)[-C]),Inf)
  matrix(apply(Z,2L,function(z) as.integer(cut(z,breaks=cuts,labels=FALSE))),n,R)
}
ds_row_prob <- function(C,truth,acc,sev){cats<-seq_len(C);ctr<-truth+sev;err<-exp(-0.5*((cats-ctr)/0.85)^2);err[truth]<-0;if(sum(err)<=0)err<-rep(1,C);err<-err/sum(err);p<-(1-acc)*err;p[truth]<-p[truth]+acc;p/sum(p)}
simulate_ds <- function(n, R) {               # study 29 ds5x6_severity
  C<-5L; p<-c(0.38,0.27,0.18,0.11,0.06); acc<-c(0.86,0.82,0.78,0.74,0.70,0.66); sev<-c(-0.65,-0.35,-0.10,0.15,0.40,0.70)
  truth<-sample.int(C,n,replace=TRUE,prob=p); X<-matrix(NA_integer_,n,R)
  for (j in seq_len(R)) for (i in seq_len(n)) X[i,j]<-sample.int(C,1L,prob=ds_row_prob(C,truth[i],acc[j],sev[j]))
  X
}
simulate_moss <- function(n, R, target) {     # study 29 knowledge model, C=5
  kp<-sqrt(target); truth<-sample.int(5L,n,replace=TRUE); X<-matrix(NA_integer_,n,R)
  for (j in seq_len(R)){k<-runif(n)<kp; g<-sample.int(5L,n,replace=TRUE); X[,j]<-ifelse(k,truth,g)}
  X
}

dgps <- list(
  gauss_normal  = list(gen=function(n) gen_continuous(n,6,"normal"),  R=6, truth_n=2e5, big=12000),
  gauss_skew    = list(gen=function(n) gen_continuous(n,6,"skewed"),  R=6, truth_n=2e5, big=12000),
  gauss_mixture = list(gen=function(n) gen_mixture(n,6),              R=6, truth_n=2e5, big=12000),
  cat_skill     = list(gen=function(n) gen_categorical(n,5),          R=5, truth_n=2e5, big=12000),
  cat_copula_hi = list(gen=function(n) simulate_copula(n,6,0.72,rep(0.2,5)), R=6, truth_n=2e5, big=12000),
  cat_copula_sk = list(gen=function(n) simulate_copula(n,6,0.55,c(.44,.25,.16,.10,.05)), R=6, truth_n=2e5, big=12000),
  cat_ds_sev    = list(gen=function(n) simulate_ds(n,6),              R=6, truth_n=6e4,  big=6000),
  cat_moss_80   = list(gen=function(n) simulate_moss(n,5,0.80),       R=5, truth_n=2e5, big=12000)
)

## ---- the SAME functional both estimators target ----------------------------
conger_q <- function(m, Sig) {
  R<-length(m); nu<-de<-0
  for (i in 1:(R-1)) for (j in (i+1):R){d2<-(m[i]-m[j])^2
    nu<-nu+Sig[i,i]+Sig[j,j]-2*Sig[i,j]+d2; de<-de+Sig[i,i]+Sig[j,j]+d2}
  1-nu/de
}
est_ntfiml <- function(X) {
  f <- tryCatch(kappa(X, estimator="nt_fiml", weight="quadratic"), error=function(e) NULL)
  if (is.null(f)) NA else unname(coef(f)["Conger"])
}
est_struct <- function(X) {                    # ten Hove structured, raters FIXED
  R <- ncol(X)
  d <- data.frame(subj=factor(c(row(X))), rat=factor(c(col(X)), levels=seq_len(R)), y=c(X))
  d <- d[!is.na(d$y), ]
  m <- tryCatch(suppressWarnings(lmer(y ~ rat + (1|subj), d, REML=FALSE,
                control=lmerControl(calc.derivs=FALSE))), error=function(e) NULL)
  if (is.null(m)) return(NA)
  fe <- fixef(m); mu_r <- numeric(R); mu_r[1] <- fe[["(Intercept)"]]
  for (j in 2:R) mu_r[j] <- fe[["(Intercept)"]] + fe[[paste0("rat", j)]]
  s2s <- as.numeric(VarCorr(m)$subj[1]); s2sr <- sigma(m)^2
  Sig <- matrix(s2s, R, R); diag(Sig) <- s2s + s2sr
  conger_q(mu_r, Sig)
}
ampute_mcar <- function(X, rate=0.4) { X[matrix(runif(length(X))<rate, nrow(X), ncol(X))] <- NA; X }

## ---- run -------------------------------------------------------------------
S <- 80; B <- if (exists("B_OVERRIDE")) B_OVERRIDE else 150
rows <- list()
for (nm in names(dgps)) {
  g <- dgps[[nm]]
  set.seed(9999); Yt <- g$gen(g$truth_n); tru <- conger_q(colMeans(Yt), cov(Yt))
  set.seed(2); Xb <- ampute_mcar(g$gen(g$big))               # asymptotic (large n, MCAR)
  asy_nt <- est_ntfiml(Xb); asy_st <- est_struct(Xb)
  M <- t(vapply(seq_len(B), function(k) {                    # finite n
    set.seed(5000+k); X <- ampute_mcar(g$gen(S))
    c(nt=est_ntfiml(X), st=est_struct(X))
  }, numeric(2)))
  M <- M[stats::complete.cases(M), , drop=FALSE]
  pe <- function(x) 100*(mean(x)-tru)/tru
  rmse <- function(x) sqrt(mean((x-tru)^2))
  message(sprintf("  %-14s truth=%.3f reps=%d", nm, tru, nrow(M)))
  rows[[nm]] <- data.frame(
    DGP=nm, truth=round(tru,3), R=g$R, reps=nrow(M),
    estimator=c("nt_fiml(sat)","struct(tenHove)"),
    asymp_pct = round(c(100*(asy_nt-tru)/tru, 100*(asy_st-tru)/tru),1),
    bias_pct  = round(c(pe(M[,"nt"]), pe(M[,"st"])),1),
    SD        = round(c(sd(M[,"nt"]), sd(M[,"st"])),4),
    RMSE      = round(c(rmse(M[,"nt"]), rmse(M[,"st"])),4))
}
out <- do.call(rbind, rows)
cat(sprintf("\nEstimand-matched: both -> quadratic Conger kappa via conger_q(mu,Sigma)\nS=%d, MCAR 40%%, B=%d; asymp = large-n complete-ish check\n\n", S, B))
print(out, row.names=FALSE)
