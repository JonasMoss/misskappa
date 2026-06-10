set.seed(2); n <- 40; C <- 4; R <- 4
z <- sample(1:C, n, TRUE)
rate <- function(z) ifelse(runif(n) < 0.7, z, sample(1:C, n, TRUE))
x <- sapply(1:R, function(j) rate(z))
pairs <- t(combn(R, 2))
for (i in 1:n) { p <- pairs[sample(nrow(pairs),1),]; x[i, setdiff(1:R,p)] <- NA }

cells <- as.matrix(expand.grid(rep(list(1:C), R)))
compat <- lapply(1:n, function(i) {
  obs <- !is.na(x[i,])
  which(colSums(t(cells[, obs, drop=FALSE]) == x[i, obs]) == sum(obs))
})

# Constraint matrix A: one row per (pattern, observed-cell) margin, over all
# 6 pair patterns x 16 value combos
A <- NULL
for (k in 1:nrow(pairs)) {
  a <- pairs[k,1]; b <- pairs[k,2]
  for (va in 1:C) for (vb in 1:C)
    A <- rbind(A, as.numeric(cells[,a]==va & cells[,b]==vb))
}
cat("cells:", ncol(A), " margin rows:", nrow(A), " rank(A):", qr(A)$rank,
    " => face dim >=", ncol(A) - qr(A)$rank, "\n")

conger <- function(th) {
  po <- 0; pe <- 0
  m <- sapply(1:R, function(j) sapply(1:C, function(c) sum(th[cells[,j]==c])))
  for (k in 1:nrow(pairs)) {
    po <- po + sum(th[cells[,pairs[k,1]]==cells[,pairs[k,2]]])
    pe <- pe + sum(m[,pairs[k,1]]*m[,pairs[k,2]])
  }
  po <- po/nrow(pairs); pe <- pe/nrow(pairs); (po-pe)/(1-pe)
}

# Smoothed EM: MAP with Dirichlet(1+delta) -> analytic center of MLE face
em <- function(th, delta = 0, iters = 50000, tol = 1e-13) {
  K <- length(th)
  for (it in 1:iters) {
    new <- numeric(K)
    for (i in 1:n) { ci <- compat[[i]]; w <- th[ci]; new[ci] <- new[ci] + w/sum(w) }
    new <- (new + delta) / (n + delta*K)
    if (max(abs(new-th)) < tol) return(new)
    th <- new
  }
  th
}

for (delta in c(0, 1e-4, 1e-3)) {
  ks <- sapply(1:6, function(s) {
    set.seed(100+s); th0 <- rgamma(C^R, runif(1,0.2,3)); th0 <- th0/sum(th0)
    conger(em(th0, delta))
  })
  cat(sprintf("delta=%g: kappa mean=%.6f  spread across starts=%.2e\n",
              delta, mean(ks), diff(range(ks))))
}
