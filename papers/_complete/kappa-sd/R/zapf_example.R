#' Generate confidence intervals and sample estimates for the Zapf example.
#' @param g The g parameter.
#' @param type Type of chance agreement, 'cohen' or 'fleiss'.
#' @return Table of confidence limits and estimates.

zapf_table = function(g, type = "cohen") {
  ci = round(rbind(agreer::kappa(agreer::dat.zapf2016, type = type, disagreement = "nominal", g = g)$conf.int,
                   agreer::kappa(agreer::dat.zapf2016, type = type, disagreement = "absolute", g = g)$conf.int,
                   agreer::kappa(agreer::dat.zapf2016, type = type, disagreement = "quadratic", g = g)$conf.int,
                   agreer::kappa(agreer::dat.zapf2016, type = type, disagreement = "hubert", g = g)$conf.int), 3)
  colnames(ci) = c("0.05", "0.95")

  est = round(rbind(agreer::kappa(agreer::dat.zapf2016, type = type, disagreement = "nominal", g = g)$estimate,
                    agreer::kappa(agreer::dat.zapf2016, type = type, disagreement = "absolute", g = g)$estimate,
                    agreer::kappa(agreer::dat.zapf2016, type = type, disagreement = "quadratic", g = g)$estimate,
                    agreer::kappa(agreer::dat.zapf2016, type = type, disagreement = "hubert", g = g)$estimate), 3)

  cbind(ci, est)[, 1:3]
}


# These are the valus used in the paper.
zapf_table(2)
zapf_table(3)
zapf_table(4)
zapf_table(4, type = "fleiss")