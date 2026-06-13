# Probe 39 — robustness sweep: does the C++ MLE-RE stay glued to lme4 across
# the conditions ten Hove et al. stress (small S, heavy missingness)?

suppressMessages({ library(Rcpp); library(lme4) })
here <- "experiments/probes/39-ten-hove-icc"
sourceCpp(file.path(here, "icc_fiml.cpp"))

fit_icc <- function(subj, rat, y, S, R) {
  nll <- function(p) { v <- exp(p); -icc_loglik(v[1],v[2],v[3],subj,rat,y,S,R)$loglik }
  v <- exp(optim(log(rep(var(y)/3,3)), nll, method="BFGS",
                 control=list(reltol=1e-12, maxit=500))$par)
  c(s2s=v[1], s2r=v[2], s2sr=v[3], A1=v[1]/sum(v), C1=v[1]/(v[1]+v[3]))
}

sim <- function(seed, S, R, per_subj, s2s=2, s2r=1, s2sr=1.5, mu=3) {
  set.seed(seed)
  a <- rnorm(S,0,sqrt(s2s)); b <- rnorm(R,0,sqrt(s2r))
  do.call(rbind, lapply(seq_len(S), function(s){
    rs <- sort(sample.int(R, per_subj))
    data.frame(subj=s, rat=rs, y=mu+a[s]+b[rs]+rnorm(per_subj,0,sqrt(s2sr)))
  }))
}

grid <- expand.grid(S=c(30,200), R=c(5,12), per=c(2,3), seed=1:5)
grid <- grid[grid$per < grid$R, ]
res <- lapply(seq_len(nrow(grid)), function(i){
  g <- grid[i,]; d <- sim(g$seed, g$S, g$R, g$per)
  mine <- fit_icc(d$subj, d$rat, d$y, g$S, g$R)
  m  <- suppressWarnings(lmer(y~1+(1|subj)+(1|rat), data=d, REML=FALSE))
  vc <- as.data.frame(VarCorr(m))
  l  <- c(vc$vcov[vc$grp=="subj"], vc$vcov[vc$grp=="rat"], vc$vcov[vc$grp=="Residual"])
  l  <- c(l, l[1]/sum(l), l[1]/(l[1]+l[3]))
  data.frame(S=g$S, R=g$R, per=g$per, miss=round(100*(1-g$per/g$R)),
             max_rel = max(abs(mine-l)/pmax(abs(l),1e-6)),
             singular = isSingular(m))
})
res <- do.call(rbind, res)
agg <- aggregate(cbind(max_rel)~S+R+per+miss, res,
                 FUN=function(x) max(x))
agg <- agg[order(agg$miss, agg$S), ]
cat("Worst relative disagreement vs lme4 across 5 seeds per cell:\n")
print(agg, row.names=FALSE, digits=3)
cat(sprintf("\nlme4 singular fits: %d / %d\n", sum(res$singular), nrow(res)))
cat(sprintf("OVERALL worst rel-diff vs lme4 (non-singular): %.2e\n",
            max(res$max_rel[!res$singular])))
