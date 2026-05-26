# Load the tibble library for creating the final data frame
library(tibble)

# --- Data Generation ---
# Each row is created as a vector using c(). We'll then combine them.

# (Part 1: Original and Duplicated Rows)
# These rows are "complete" (sum to 6) and include duplicates to satisfy (b).
complete_and_duplicates <- list(
  c(0, 0, 0, 6, 0),  # Original row 1 from dat.fleiss1971
  c(0, 3, 0, 0, 3),  # Original row 2
  c(0, 0, 0, 0, 6),  # Original row 4
  c(0, 0, 0, 0, 6),  # DUPLICATE of row 4 to satisfy (b)
  c(2, 0, 4, 0, 0),  # Original row 6
  c(2, 0, 3, 1, 0),  # A nice mix from original row 8
  c(5, 1, 0, 0, 0),  # Original row 18
  c(4, 0, 0, 0, 2),  # Original row 27
  c(4, 0, 0, 0, 2)   # DUPLICATE of row 27 to satisfy (b)
)

# (Part 2: Rows with Missing Values)
# This section creates rows with missing values as per your definition (row sum != 6).
# We ensure variety by covering all patterns from 1 to 5 missing values, satisfying (c).

# Rows with 1 missing value (sum = 5)
sum_is_5 <- list(
  c(0, 0, 0, 5, 0),  # Modified from (0,0,0,6,0)
  c(0, 2, 0, 0, 3),  # Modified from (0,3,0,0,3)
  c(1, 0, 4, 0, 0),  # Modified from (2,0,4,0,0)
  c(0, 4, 0, 1, 0),   # Modified from (0,5,0,1,0)
  c(0, 4, 0, 1, 0), # DUP
  c(0, 4, 0, 1, 0) # DUP
)

# Rows with 2 missing values (sum = 4)
sum_is_4 <- list(
  c(0, 2, 0, 2, 0),  # Modified from (0,3,0,3,0)
  c(0, 0, 4, 0, 0),  # Modified from (0,0,4,0,2)
  c(3, 1, 0, 0, 0),  # Modified from (5,1,0,0,0)
  c(1, 0, 2, 1, 0),   # Modified from (2,0,3,1,0)
  c(1, 0, 2, 1, 0)   # Modified from (2,0,3,1,0)
)

# Rows with 3 missing values (sum = 3)
sum_is_3 <- list(
  c(0, 1, 0, 2, 0),  # Modified from (0,1,0,5,0)
  c(1, 0, 1, 1, 0),  # Modified from (2,0,3,1,0)
  c(0, 2, 1, 0, 0)   # Modified from (0,3,3,0,0)
)

# Rows with 4 missing values (sum = 2)
sum_is_2 <- list(
  c(0, 2, 0, 0, 0),  # Modified from (0,3,0,0,3)
  c(2, 0, 0, 0, 0),  # Modified from (2,0,3,1,0)
  c(0, 0, 0, 0, 2)   # Modified from (0,0,0,0,6)
)

# Rows with 5 missing values (sum = 1)
sum_is_1 <- list(
  c(1, 0, 0, 0, 0),  # Modified from (5,1,0,0,0)
  c(0, 0, 1, 0, 0),  # Modified from (1,0,5,0,0)
  c(0, 0, 0, 1, 0)   # Modified from (0,0,0,6,0)
)


# --- Assemble into a Tibble ---
# Combine all the lists of row vectors
all_rows <- c(
  complete_and_duplicates,
  sum_is_5,
  sum_is_4,
  sum_is_3,
  sum_is_2,
  sum_is_1
)

my_data <- as_tibble(do.call(rbind, all_rows))
colnames(my_data) <- c("depression", "personality_disorder", "schizophrenia", "neurosis", "other")

x <- as.matrix(my_data)
storage.mode(x) <- "integer"


get_extensions <- \(out) {
  active_ranks <- out$idx_to_rank

  active_patterns <- unrank_rcpp(active_ranks, r = 6, c = 5)

  colnames(active_patterns) <- colnames(my_data)

  tibble::as_tibble(active_patterns) %>%
    dplyr::mutate(
      rank = active_ranks,
      theta = out$theta,
      .before = 1
    ) %>%
    dplyr::arrange(dplyr::desc(theta))

}

out <- run_em_rcpp(x, 6, max_iter = 10000, alpha = 1e-4)
out$iterations
out$theta
dplyr::arrange(get_extensions(out), rank)
print(my_data, n = 29)


x <- as.matrix(dat.fleiss1971)
storage.mode(x) <- "integer"
r <- 6
em_output <- run_em_rcpp(x, r = r)
final_kappa <- calculate_kappa_rcpp(em_output, r = r, c = 5)
final_kappa$kappa


sapply(6:12, \(r) {
  em_output <- run_em_rcpp(x, r = r, max_iter = 100000)
  final_kappa <- calculate_kappa_rcpp(em_output, r = r, c = 5)
  final_kappa$kappa
})




sqrt(final_kappa$kappa_var)

final_kappa$kappa
irrCAC::fleiss.kappa.dist(x)$coeff

n <- nrow(x)
sqrt(final_kappa$kappa_var)
irrCAC::fleiss.kappa.dist(x)$stderr * sqrt(n - 1) / sqrt(n)
