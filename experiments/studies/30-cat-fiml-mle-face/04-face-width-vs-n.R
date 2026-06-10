# Face width of Conger kappa vs n, C=4 R=4 pair-rotation design.
# Grouped E-step: at most 6*16=96 distinct observed records.
C <- 4; R <- 4
cells <- as.matrix(expand.grid(rep(list(1:C), R)))
K <- nrow(cells)
prs <- t(combn(R, 2))

conger <- function(th) {
  po <- 0; pe <- 0
  m <- sapply(1:R, function(j) sapply(1:C, function(c) sum(th[cells[,j]==c])))
  for (k in 1:nrow(prs)) {
    po <- po + sum(th[cells[,prs[k,1]]==cells[,prs[k,2]]])
    pe <- pe + sum(m[,prs[k,1]]*m[,prs[k,2]])
  }
  po <- po/nrow(prs); pe <- pe/nrow(prs); (po-pe)/(1-pe)
}

run_one <- function(n, seed, nstart = 8, delta = 0) {
  set.seed(seed)
  z <- sample(1:C, n, TRUE)
  rate <- function(z) ifelse(runif(n) < 0.7, z, sample(1:C, n, TRUE))
  x <- sapply(1:R, function(j) rate(z))
  for (i in 1:n) { p <- prs[sample(nrow(prs),1),]; x[i, setdiff(1:R,p)] <- NA }
  # group identical observed records
  key <- apply(x, 1, paste, collapse = "/")
  tab <- table(key)
  groups <- lapply(names(tab), function(k) {
    row <- suppressWarnings(as.integer(strsplit(k, "/", fixed=TRUE)[[1]]))
    obs <- !is.na(row)
    which(colSums(t(cells[, obs, drop=FALSE]) == row[obs]) == sum(obs))
  })
  w <- as.numeric(tab)
  em <- function(th, iters = 200000, tol = 1e-13) {
    for (it in 1:iters) {
      new <- numeric(K)
      for (g in seq_along(groups)) {
        ci <- groups[[g]]; p <- th[ci]; new[ci] <- new[ci] + w[g]*p/sum(p)
      }
      new <- (new + delta)/(n + delta*K)
      if (max(abs(new-th)) < tol) return(new)
      th <- new
    }
    th
  }
  ks <- sapply(1:nstart, function(s) {
    set.seed(1000*seed + s)
    th0 <- rgamma(K, runif(1, 0.2, 3)); th0 <- th0/sum(th0)
    conger(em(th0))
  })
  diff(range(ks))
}

for (n in c(40, 160, 640)) {
  widths <- sapply(1:5, function(seed) run_one(n, seed))
  cat(sprintf("n=%4d  median width=%.2e  widths: %s\n",
              n, median(widths), paste(sprintf("%.1e", widths), collapse=" ")))
}
