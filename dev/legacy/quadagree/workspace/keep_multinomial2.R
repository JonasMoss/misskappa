x <- as.matrix(dat.fleiss1971)

n <- nrow(x)
c(x = sqrt(compute_kappa(x, raters = 6)$kappa_var),
  y = irrCAC::fleiss.kappa.dist(x)$stderr)

result <- t(sapply(4:30, \(i) {
x <- as.matrix(dat.fleiss1971)[1:i, ]
n <- nrow(x)
c(x = sqrt(compute_kappa(x, raters = 6)$kappa_var),
y = irrCAC::fleiss.kappa.dist(x)$stderr)
}))


plot(result)


x <- as.matrix(dat.fleiss1971)
x[1, ] <- c(0,0,0,5,0)

n <- nrow(x)

# Our result
our_stderr <- sqrt(compute_kappa_final(x, raters = 6)$kappa_var)

# irrCAC result, adjusted for n vs n-1
irrCAC_stderr <- irrCAC::fleiss.kappa.dist(x)$stderr * sqrt(n - 1) / sqrt(n)

print(paste("Our Stderr:    ", round(our_stderr, 6)))
print(paste("irrCAC Stderr: ", round(irrCAC_stderr, 6)))

# Let's see the var_phi matrix on complete data - it should NOT be diagonal
fit <- em_counts_constrained(x, raters = 6)
print("Resulting var_phi (should be non-diagonal):")
print(round(fit$var_phi, 6))
