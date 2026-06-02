#' Fleiss (1971) psychiatric diagnoses dataset
#'
#' Counts-format rating data from Fleiss (1971), n = 30 subjects, six raters,
#' five categories. Each row is one subject; each column counts the number of
#' raters who assigned that subject to that category.
#'
#' @format A 30 x 5 integer matrix.
#' @source Fleiss, J. L. (1971). Measuring nominal scale agreement among many
#'   raters. \emph{Psychological Bulletin}, 76(5), 378-382.
#' @keywords datasets
"dat.fleiss1971"

#' Gwet (2014) example dataset
#'
#' Example raw rating data carried over from the legacy `misskappa`
#' package; used in Gwet (2014).
#'
#' @keywords datasets
"dat.gwet2014"

#' Klein (2018) example dataset
#'
#' Example raw rating data carried over from the legacy `misskappa` package.
#'
#' @keywords datasets
"dat.klein2018"

#' Zapf (2016) example dataset
#'
#' Example raw rating data carried over from the legacy `misskappa` package.
#'
#' @keywords datasets
"dat.zapf2016"

#' McDuff & Girard (2019) smiling-images agreement study
#'
#' Agreement study from McDuff & Girard (2019) in long format: MTurk judges
#' rate images of people smiling. The design is large and incomplete --- many
#' items and judges, and no judge rates every image --- which makes it a useful
#' missing-data example. There are 273 items, 121 judges, and 1365 rows.
#'
#' @format A 1365 x 4 data frame with one row per (item, judge) rating:
#'   \describe{
#'     \item{rating_positive}{numeric; rated positivity of the image.}
#'     \item{rating_smile}{numeric; rated presence of a smile in the image.}
#'     \item{item}{integer; the image rated.}
#'     \item{judge}{integer; the MTurk judge giving the rating.}
#'   }
#' @source The file \code{04_MTurk.csv} at \url{https://osf.io/n4grd/}.
#' @references McDuff, D., & Girard, J. M. (2019). Democratizing Psychological
#'   Insights from Analysis of Nonverbal Behavior. \emph{2019 8th International
#'   Conference on Affective Computing and Intelligent Interaction (ACII)},
#'   220-226. \doi{10.1109/ACII.2019.8925503}
#' @keywords datasets
"dat.mcduff2019"
