x = rbind(
  c(1, 1, 2, 1, 1),
  c(1, 2, 3, 2, 2),
  c(2, 1, 1, 1, 1),
  c(2, 3, 4, 4, 5))

agreer::kappa(x, g = 5, disagreement = "absolute")

mus = apply(x, 2, median)
vs = sapply(seq_along(mus), function(i) mean(abs(x[, i] - mus[i])))
