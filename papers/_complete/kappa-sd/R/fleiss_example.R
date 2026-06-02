fleiss_table = function(g) {
  ci = round(rbind(agreer::kappa(agreer::dat.fleiss1971, type = "fleiss", disagreement = "nominal", g = g, fleiss_form = TRUE)$conf.int,
                   agreer::kappa(agreer::dat.fleiss1971, type = "fleiss", disagreement = "hubert", g = g, fleiss_form = TRUE)$conf.int), 3)

  est = round(rbind(agreer::kappa(agreer::dat.fleiss1971, type = "fleiss", disagreement = "nominal", g = g, fleiss_form = TRUE)$estimate,
                    agreer::kappa(agreer::dat.fleiss1971, type = "fleiss", disagreement = "hubert", g = g, fleiss_form = TRUE)$estimate), 3)

  cbind(ci, est)[, 1:3]
}

fleiss_table(2)
fleiss_table(3)
fleiss_table(6)