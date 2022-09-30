#' Capron transition onsets (dust)
#'
#' Abrupt warmings observed in the log(Ca2+) record of the NGRIP ice core. Taken from
#' Supplementary data 4 in Capron et al. (2021). This data set is only used to compare
#' our transition dating estimation in "Examples/Paper_results/Paper_results_1.Rmd"
#' and "Examples/Paper_results/Paper_results_2.Rmd".
#'
#' @format A \code{tbl data.frame} with 25 rows and 16 columns. The following columns
#' are relevant to analyses carried out in this package:
#' \describe{
#'   \item{`NGRIP dust`}{Name of the transition.}
#'   \item{`t1 (50\%)`}{The median of the transition onset in years before present
#'   (used as point estimate).}
#'   \item{`d1 (50\%)`}{The median of the transition onset depth in meters.}
#'   }
#' @references Capron, E., Rasmussen, S.O., Popp, T.J. et al. The anatomy of past
#' abrupt warmings recorded in Greenland ice. Nat Commun 12, 2106 (2021).
#' https://doi.org/10.1038/s41467-021-22241-w



"capron_onsets_NGRIP_dust"
