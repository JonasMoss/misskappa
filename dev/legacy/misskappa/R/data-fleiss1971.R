#' Agreement study from Fleiss, J. L. (1971)
#'
#' Fleiss, J. L. (1971) is best known for introducing Fleiss' kappa. In Table 1
#'    he presents diagnosis data of psychiatric patients. The data is on Fleiss
#'    form, with `n = 30` patients diagnosed by `6` psychiatrists each. It is
#'    not the case that the same set of six psychiatrists diagnosed every
#'    patient.
#'
#' @usage dat.fleiss1971
#'
#' @format The tibble contains the five columns "depression",
#'     "personality disorder", "schizophrenia" "neurosis", and "other".
#'     The content of the `ij`th cell is the number of raters who diagnosed
#'     the `i`th patient as having mental illness `j`.
#'
#' @keywords datasets
#'
#' @references
#' Fleiss, J. L. (1971). Measuring nominal scale agreement among many raters.
#' Psychological Bulletin, 76(5), 378<U+2013>382.
#' https://doi.org/10.1037/h0031619
#'
#' @source
#' The data was scraped by hand from Fleiss (1971).
"dat.fleiss1971"
