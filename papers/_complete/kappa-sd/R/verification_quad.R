### ============================================================================
###  Verifying computations of the Fleiss example.
### ============================================================================

#source("disaggregate.R")

#' Quadratic distance function.
d_quadratic <- \(x) {
  mean((x - mean(x))^2)
}

#' Calculate the g-wise Hubert's kappa
#' @param x,g Data and g.
#' @returns g-wise Hubert's kappa
r_quad <- \(x, g = 2) {
  x <- as.matrix(x)
  n <- nrow(x)
  r <- ncol(x)
  size <- length(unique(c(x)))

  indices <- arrangements::combinations(r, g)

  da <- mean(sapply(seq_len(nrow(indices)), \(i) {
    y <- x[, indices[i, ]]
    mean(sapply(seq_len(n), \(i) {
      d_quadratic(y[i, ])
    }))
  }))

  # Calculate df
  uniques <- sort(unique(c(x)))
  pk <- apply(x, 2, function(y) Rfast::Table(c(y, uniques))) - 1
  pk <- pk / (nrow(x))
  i_combs <- as.matrix(expand.grid(lapply(seq(g), \(x) seq(size))))
  dists <- apply(i_combs, 1, d_quadratic)
  r_combs <- arrangements::combinations(r, g)

  probs <- pk[i_combs[, 1], r_combs[, 1]]
  for(i in 2:g) {
    probs <- probs * pk[i_combs[, i], r_combs[, i]]
  }

  df <- mean(probs * dists)

  c(da, 1 - da / (df*nrow(i_combs)))
}

### ============================================================================
###  Making mus
### ============================================================================

x <- as.matrix(agreer::dat.zapf2016)
g <- 3

mus <- \(x, g) {
  r <- ncol(x)
  r_combs <- arrangements::combinations(r, g)
  uniques <- sort(unique(c(x)))
  size <- length(uniques)
  i_combs <- as.matrix(expand.grid(lapply(seq(g - 1), \(x) seq(size))))
  nk <- apply(x, 2, function(y) Rfast::Table(c(y, uniques))) - 1
  dists <- list()

  for(elem in uniques) {
    dists[[elem]] <- sapply(seq(nrow(i_combs)), \(i) d_quadratic(c(elem, i_combs[i, ])))
  }


  nk <- apply(x, 2, function(y) Rfast::Table(c(y, uniques))) - 1
  pk <- nk / nrow(x)
  mus <- rep(0, nrow(x))

  for(i in seq(nrow(x))) {
    ## We need to iterate over all rater combinations while caring about which
    ## rater comes first. The outer loop iterates over ordered unordered
    ## combinations (lexicographically ordered), the inner loop lets us work
    ## with which rater comes first.

    out = 0
    for(j in seq(nrow(r_combs))) {
        for(k in r_combs[j, ]) {
          raters <- setdiff(r_combs[j, ], k)
          ps <- pk[i_combs[, 1], raters[1]]

          if (g > 2) {
            for(l in 2:(g-1)) {
              ps <- ps * pk[i_combs[, l], raters[l]]
            }
          }

          elem <- x[i, k]
          out = out + mean(dists[[elem]] * ps)
        }
    }
    mus[i] = out / (length(r_combs))
  }
  mus * nrow(i_combs)
}







g <- 2
mean(mus(x, g))
r_quad(x, g)[1] / mean(mus(x, g)) / (1 - r_quad(x, g)[2])







### ============================================================================
###  Examples
### ============================================================================

x <- as.matrix(agreer::dat.zapf2016)
cpp_quad<- \(x, g) {
  unname(agreer::kappa(x, type = "cohen", d = "quadratic", g = g)$estimate[1])
}


mean(mus(x, 2))

r_quad(x, 2)
cpp_quad(x, 2)

r_quad(x, 3)
cpp_quad(x, 3)

r_quad(x, 4)
cpp_quad(x, 4)

microbenchmark::microbenchmark(mus(rbind(x, x), 3), cpp_quad(rbind(x, x), 3))