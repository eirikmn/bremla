#' Tie-point probability distribution (Adolphi et al.)
#'
#' Probability distributions for the five tie-points described in Adolphi et al. (2018), and updated in Muscheler et al. (2020).
#' The tiepoints are given in GICC05 ages (yb2k) 11,050, 12,050, 13,050, 22,050 and 42,050, with corresponding NGRIP dephts (m)
#' 1492.50, 1502.8, 1532.7, 1767.25 and 2132.4.
#'
#' @format A \code{data.frame} with 5905 rows and 3 columns (five groups):
#' \describe{
#'   \item{x}{The x-axis (age) for the distribution.}
#'   \item{y}{The y-axis (density) for the distribution.}
#'   \item{group}{Integer indicating to which tie-point each point is associated with.}
#'   }
#' @references Adolphi, F., Bronk Ramsey, C., Erhardt, T., Edwards, R. L., Cheng, H., Turney, C. S. M., Cooper, A., Svensson, A., Rasmussen, S. O., Fischer, H., and Muscheler, R. (2018).
#' Connecting the Greenland ice-core and U∕Th timescales via cosmogenic radionuclides: testing the synchroneity of Dansgaard–Oeschger events,
#' Clim. Past, 14, 1755–1781, https://doi.org/10.5194/cp-14-1755-2018
#' @references Muscheler, R., Adolphi, F., Heaton, T., Bronk Ramsey, C., Svensson, A., Van der Plicht, J., & Reimer, P. (2020).
#' Testing and Improving the IntCal20 Calibration Curve with Independent Records.
#' Radiocarbon, 62(4), 1079-1094. doi:10.1017/RDC.2020.54


"adolphi_tiepoints"
