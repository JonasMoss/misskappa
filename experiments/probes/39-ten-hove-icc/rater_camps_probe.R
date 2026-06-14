# Probe 39 — the realistic misspecification: rater "camps" disagree on a subject
# attribute (e.g. some raters score a subgroup higher, others lower). This makes
# the rater covariance BLOCK-structured, not compound-symmetric, so the ten Hove
# structured model is misspecified -- while the saturated nt_fiml is not.
#
# Jointly Gaussian by construction (so nt_fiml is correctly specified):
#   y_sj = mu + b_j + a_s + beta_j * w_s + e_sj,
#   a_s ~ N(0,s2a)  (shared subject level),
#   w_s ~ N(0,s2w)  (latent subject attribute, e.g. group),
#   beta_j = +/- delta  (rater camp: half score w higher, half lower),
#   e_sj ~ N(0,s2e).
# Cov(y_i,y_j) = s2a + delta^2 * camp_i*camp_j * s2w  -> block, NOT compound sym.
#
# Prediction: under MCAR both consistent (kappa needs only pooled moments);
# under anchor-MAR the structured model breaks, the saturated nt_fiml holds.

suppressMessages({ library(misskappa); library(lme4) })

R <- 6; camp <- c(1, 1, 1, -1, -1, -1)            # three raters each way
mu <- 3; s2a <- 1; s2w <- 1; s2e <- 1; delta <- 1
bj <- seq(-0.3, 0.3, length.out = R)              # mild rater main effects

gen <- function(n) {
  a <- rnorm(n, 0, sqrt(s2a)); w <- rnorm(n, 0, sqrt(s2w))
  mu + outer(a, rep(1, R)) + matrix(bj, n, R, byrow = TRUE) +
    outer(w, delta * camp) + matrix(rnorm(n*R, 0, sqrt(s2e)), n, R)
}
conger_q <- function(m, S) { R<-length(m); nu<-de<-0
  for (i in 1:(R-1)) for (j in (i+1):R){d2<-(m[i]-m[j])^2
    nu<-nu+S[i,i]+S[j,j]-2*S[i,j]+d2; de<-de+S[i,i]+S[j,j]+d2}; 1-nu/de }

Sig_pop <- s2a + (delta^2*s2w)*outer(camp, camp); diag(Sig_pop) <- diag(Sig_pop) + s2e
truth <- conger_q(mu + bj, Sig_pop)
cat(sprintf("Block covariance (within-camp C=%.1f, cross-camp C=%.1f, V=%.1f)\n",
            s2a+delta^2*s2w, s2a-delta^2*s2w, s2a+delta^2*s2w+s2e))
cat(sprintf("truth (quadratic Conger kappa) = %.4f\n\n", truth))

ampute <- function(X, mech, rate = 0.4) {
  n<-nrow(X); R<-ncol(X)
  if (mech=="mcar") X[matrix(runif(n*R)<rate,n,R)] <- NA
  else { a<-log(rate/(1-rate)); pr<-plogis(a + 1.1*scale(X[,1])[,1])  # MAR on anchor rater 1
         for (j in 2:R) X[runif(n)<pr, j] <- NA }
  X
}
est_nt <- function(X){f<-tryCatch(kappa(X,estimator="nt_fiml",weight="quadratic"),error=function(e)NULL)
  if(is.null(f))NA else unname(coef(f)["Conger"])}
est_st <- function(X){R<-ncol(X)
  d<-data.frame(subj=factor(c(row(X))),rat=factor(c(col(X)),levels=seq_len(R)),y=c(X)); d<-d[!is.na(d$y),]
  m<-tryCatch(suppressWarnings(lmer(y~rat+(1|subj),d,REML=FALSE,control=lmerControl(calc.derivs=FALSE))),error=function(e)NULL)
  if(is.null(m))return(NA)
  fe<-fixef(m); mr<-numeric(R); mr[1]<-fe[["(Intercept)"]]; for(j in 2:R) mr[j]<-fe[["(Intercept)"]]+fe[[paste0("rat",j)]]
  s2s<-as.numeric(VarCorr(m)$subj[1]); s2sr<-sigma(m)^2
  Sg<-matrix(s2s,R,R); diag(Sg)<-s2s+s2sr; conger_q(mr,Sg)}

## asymptotic check (large n) — the clean consistency test
cat("ASYMPTOTIC (n=20000), %bias vs truth:\n")
for (mech in c("mcar","anchorMAR")) {
  set.seed(1); X <- ampute(gen(20000), mech)
  nt<-est_nt(X); st<-est_st(X)
  cat(sprintf("  %-10s  nt_fiml=%.4f (%+.1f%%)   struct=%.4f (%+.1f%%)\n",
              mech, nt, 100*(nt-truth)/truth, st, 100*(st-truth)/truth))
}

## finite-n confirm under anchor-MAR
B <- if (exists("B_OVERRIDE")) B_OVERRIDE else 200; S <- 200
M <- t(vapply(seq_len(B), function(k){set.seed(3000+k); X<-ampute(gen(S),"anchorMAR")
  c(nt=est_nt(X), st=est_st(X))}, numeric(2)))
M <- M[stats::complete.cases(M),,drop=FALSE]
cat(sprintf("\nFINITE n=%d anchor-MAR, B=%d:\n", S, nrow(M)))
cat(sprintf("  nt_fiml  bias=%+.1f%%  SD=%.4f  RMSE=%.4f\n",
            100*(mean(M[,"nt"])-truth)/truth, sd(M[,"nt"]), sqrt(mean((M[,"nt"]-truth)^2))))
cat(sprintf("  struct   bias=%+.1f%%  SD=%.4f  RMSE=%.4f\n",
            100*(mean(M[,"st"])-truth)/truth, sd(M[,"st"]), sqrt(mean((M[,"st"]-truth)^2))))
