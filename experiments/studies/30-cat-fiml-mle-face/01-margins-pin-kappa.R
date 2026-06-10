# 3 raters, 2 categories, patterns {1,2},{2,3},{1,3} only (no complete cases).
# All pairs co-observed => design guard passes, but the 2x2x2 joint is NOT
# identified: theta + t*(-1)^(i+j+k)*v preserves ALL two-way margins.
# Claim: kappa (Conger/Fleiss/BP) is constant along that face.

idx <- expand.grid(i = 0:1, j = 0:1, k = 0:1)
v <- (-1)^(idx$i + idx$j + idx$k)
theta0 <- c(0.20, 0.08, 0.07, 0.15, 0.10, 0.12, 0.08, 0.20)
stopifnot(abs(sum(theta0) - 1) < 1e-12)

kappas <- function(theta) {
  p <- array(theta, dim = c(2,2,2))
  m12 <- apply(p, c(1,2), sum); m23 <- apply(p, c(2,3), sum); m13 <- apply(p, c(1,3), sum)
  m1 <- apply(p, 1, sum); m2 <- apply(p, 2, sum); m3 <- apply(p, 3, sum)
  po <- (sum(diag(m12)) + sum(diag(m23)) + sum(diag(m13))) / 3
  pec <- (sum(m1*m2) + sum(m2*m3) + sum(m1*m3)) / 3        # Conger
  pbar <- (m1 + m2 + m3) / 3; pef <- sum(pbar^2)            # Fleiss
  c(conger = (po - pec)/(1 - pec), fleiss = (po - pef)/(1 - pef), bp = 2*po - 1)
}

for (t in c(-0.05, 0, 0.03, 0.06)) {
  th <- theta0 + t * v * 0.5
  stopifnot(all(th > 0))
  cat(sprintf("t=%+.2f  theta[1]=%.4f  ", t, th[1])); print(round(kappas(th), 10))
}
