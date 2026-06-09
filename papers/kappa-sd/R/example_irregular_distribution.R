library("agreeable")

s = sqrt(c(0.90, 0.90, 0.90))
true_dist = rep(1,5) / 5
n_sim = 20
n_reps = 200000
kappas = replicate(n_reps, {
  y = agreeable::simulate_jsm(
    n = n_sim,
    s = s,
    true_dist = true_dist,
    model = "bp")

  c(cohen_kappa(y,skip = TRUE, count = TRUE, weight = "unit"),
    cohen_kappa(y,skip = TRUE, count = TRUE, weight = "abs"))
})

z1 = head(Rfast::Table(kappas[1, ]), -1) / n_reps
z2 = head(Rfast::Table(kappas[2, ]), -1) / n_reps

sorted1 = head(sort(unique(kappas[1, ])), -1)
sorted2 = head(sort(unique(kappas[2, ])), -1)

pdf("../figures/example_irregular.pdf", width = 8, heigh = 4)
par(mfrow = c(1, 2))
plot(sorted1, z1, type = "h", xlab = expression(paste("Nominal weights ", kappa[d])), ylab = "")
plot(sorted2, z2, type = "h", xlab = expression(paste("Absolute value weights ", kappa[d])), ylab = "")
dev.off()