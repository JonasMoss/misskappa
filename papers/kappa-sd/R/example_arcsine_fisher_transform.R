s = c(0.90, 0.90, 0.90)
true_dist = rep(1,5) / 5
n_sim = 20

kappas = replicate(50000, {
  y = agreeable::simulate_jsm(
    n = n_sim,
    s = s,
    true_dist = true_dist,
    model = "bp")

  agreeable::cohen_kappa(y, weight = "quadratic", skip = TRUE)
})


pdf("figures/example_arcsine_fisher_transform.pdf")
par(mfrow = c(2, 2))
hist(kappas, xlab = expression(paste("Untransformed ", kappa[w])), main = NULL, breaks = 100, freq = FALSE)
hist(asin(kappas), xlab = expression(paste("Arcsine-transformed ", kappa[w])), main = NULL, breaks = 100, freq = FALSE)
hist(atanh(kappas), xlab = expression(paste("Fisher-transformed ", kappa[w])), main = NULL, breaks = 100, freq = FALSE)
dev.off()
