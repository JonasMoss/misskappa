#' Calculate Gwet-style Fleiss' Kappa from Scratch
#'
#' This function replicates the logic used by packages like irrCAC for handling
#' Fleiss' Kappa with missing data (NA values). It is based on the "average of
#' per-subject statistics" philosophy.
#'
#' @param ratings A matrix where rows are subjects and columns are raters.
#'   May contain NA values.
#'
#' @return A list containing three key values:
#'   - d_observed: The observed disagreement proportion (equivalent to 1 - pa).
#'   - d_chance: The chance disagreement proportion (equivalent to 1 - pe).
#'   - kappa: The final chance-corrected agreement coefficient.

calculate_gwet_kappa_from_scratch <- function(ratings) {

  n <- nrow(ratings)

  per_subject_disagreement_proportions <- c()
  subjects_with_pairs <- 0

  for (i in 1:n) {
    subj_ratings <- na.omit(ratings[i, ])
    r_i <- length(subj_ratings)

    if (r_i >= 2) {
      subjects_with_pairs <- subjects_with_pairs + 1

      rater_indices <- 1:r_i
      pairs_idx <- combn(rater_indices, 2)

      total_pairs_for_subj <- ncol(pairs_idx)
      disagreeing_pairs_for_subj <- 0

      for (p in 1:total_pairs_for_subj) {
        rating1 <- subj_ratings[pairs_idx[1, p]]
        rating2 <- subj_ratings[pairs_idx[2, p]]
        if (rating1 != rating2) {
          disagreeing_pairs_for_subj <- disagreeing_pairs_for_subj + 1
        }
      }

      subject_prop <- disagreeing_pairs_for_subj / total_pairs_for_subj
      per_subject_disagreement_proportions <- c(per_subject_disagreement_proportions, subject_prop)
    }
  }

  d_observed <- mean(per_subject_disagreement_proportions)


  categories <- sort(unique(na.omit(as.vector(ratings))))
  q <- length(categories)
  p_ik_matrix <- matrix(0, nrow = n, ncol = q)
  colnames(p_ik_matrix) <- categories

  for (i in 1:n) {
    subj_ratings <- na.omit(ratings[i, ])
    r_i <- length(subj_ratings)

    if (r_i > 0) {
      counts_for_subj <- table(factor(subj_ratings, levels = categories))
      p_ik_matrix[i, ] <- counts_for_subj / r_i
    }
  }

  p_k_vector <- colMeans(p_ik_matrix)
  pe <- sum(p_k_vector^2)
  d_chance <- 1 - pe

  return(list(
    d_observed = d_observed,
    d_chance = d_chance,
    kappa = 1 - (d_observed / d_chance)
  ))
}

calculate_gwet_kappa_from_scratch(as.matrix(dat.klein2018))

irrCAC::fleiss.kappa.raw(dat.klein2018)
