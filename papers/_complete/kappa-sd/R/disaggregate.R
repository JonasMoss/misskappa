#' Disaggregate agreement data.
#'
#' Take data on aggregated (or "Fleiss") form and construct a possible
#'    disaggregated variant.
#'
#' @param x Fleiss form data such as `dat.fleiss1971`.
#' @returns Output on "wide" form.
disaggr <- \(x) {
  r <- sum(x[1, ])
  size <- ncol(x)
  n <- nrow(x)
  out <- matrix(0, nrow = n, ncol = r)
  for (i in seq(n)) {
    current_r = 1
    current_index = 1
    while (current_index <= size) {
      if (x[i, current_index] == 0) {
        current_index = current_index + 1
      } else {
        out[i, current_r] = current_index
        x[i, current_index] = x[i, current_index] - 1
        current_r = current_r + 1
      }
    }
  }
  out
}

aggr <- \(x, cats = NULL) {
  if(is.null(cats)) cats <- sort(unique(c(as.matrix(x))))
  x <- as.matrix(x)
  results <- matrix(0, nrow(x), length(cats))
  for (i in seq_len(nrow(x))) {
    counts = Rfast::Table(as.matrix(x)[i, ])
    results[i, as.integer(names(counts))] = counts
  }
  results
}


aggr(disaggr(x)) == x
