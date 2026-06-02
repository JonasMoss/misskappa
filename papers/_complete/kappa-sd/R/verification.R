### ============================================================================
###  Verifying computations of the Fleiss example.
### ============================================================================

source("disaggregate.R")
x <- read.csv("fleiss1971.csv")

### ============================================================================
###  Functions
### ============================================================================

#' Calculate the g-wise Hubert's kappa
#'
#' Uses the equivalent percent agreement definition of Hubert's kappa.
#'
#' @param x Data.
#' @returns g-wise Hubert's kappa
r_hubert <- \(x, g = 2) {
  n <- nrow(x)
  r <- sum(x[1, ])

  # Disaggregate the data; this is used for grouping.

  x_new <- disaggr(x)

  # Indices for all g-ary combinations.
  indices <- arrangements::combinations(r, g)

  pa <- mean(sapply(seq_len(nrow(indices)), \(i) {
    mean(apply(x_new[, indices[i, ]], 1, \(x) length(unique(x)) == 1))
  }))

  ps <- rowMeans(sapply(seq_len(nrow(indices)), \(i) {
    table(x_new[, indices[i, ]]) / (g * n)
  }))

  pca <- sum(ps^g)

  (pa - pca) / (1 - pca)
}

#' Nominal distance function.
d_nominal <- \(x) {
  tab <- Rfast::Table(x)
  sum(tab) - max(tab)
}

#' Quadratic distance function.
d_quadratic <- \(x) {

}

#' Calculate the g-wise Hubert's kappa
#' @param x,g Data and g.
#' @returns g-wise Hubert's kappa
r_nominal <- \(x, g = 2) {
  n <- nrow(x)
  r <- sum(x[1, ])
  size <- ncol(x)
  # Disaggregate the data; this is used for grouping.

  x_new <- disaggr(x)

  # Indices for all g-ary combinations.
  indices <- arrangements::combinations(r, g)

  da <- mean(sapply(seq_len(nrow(indices)), \(i) {
    y <- x_new[, indices[i, ]]
    mean(sapply(seq_len(n), \(i) {
      d_nominal(y[i, ])
    }))
  }))

  # Calculate df
  counts <- sapply(seq_len(size), \(i) sum(x[, i]))
  ps <- counts / sum(counts)
  combs <- as.matrix(expand.grid(lapply(seq(g), \(x) seq(size))))
  probs <- ps[combs]
  dim(probs) <- c(length(probs) / g, g)
  probs <- apply(probs, 1, prod)
  df <- sum(apply(combs, 1, d_nominal) * probs)


  1 - da / df
}

### ============================================================================
###  Examples
### ============================================================================

#' Calculate Hubert's kappa (c++ way).
#' @param x,g Data and g.
#' @return Hubert's kappa calculate from the agreer package.
cpp_hubert <- \(x, g) {
  unname(agreer::kappa(x, type = "fleiss", d = "hubert", g = g, f = TRUE)$estimate[1])
}

#' Calculate nominal kappa (c++ way).
#' @param x,g Data and g.
#' @return Nominal kappa calculate from the agreer package.
cpp_nominal <- \(x, g) {
  unname(agreer::kappa(x, type = "fleiss", d = "nominal", g = g, f = TRUE)$estimate[1])
}


r_hubert(x, 2) # 0.4302445
cpp_hubert(x, 2) # 0.4302445

r_hubert(x, 3) # 0.3331096
cpp_hubert(x, 3) # 0.3331096

r_hubert(x, 6) # 0.1657997
cpp_hubert(x, 6) # 0.1657997

r_nominal(x, 2) # 0.4302445
cpp_nominal(x, 2) # 0.4302445

r_nominal(x, 3) # 0.4962921
cpp_nominal(x, 3) # 0.4962921

r_nominal(x, 6) # 0.4859784
cpp_nominal(x, 6) # 0.4859784