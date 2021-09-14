#' Optimal data window for DO-events
#'
#' Data windows enclosing DO-events found by optimization scheme to be used in linear ramp model fit. Also includes onset depth and dating given by Rasmussen et al. (2014).
#'
#'
#' @format A data.frame with 29 rows and 10 columns:
#' \describe{
#'   \item{ID}{The row in the \code{GISevents} dataset.}
#'   \item{onsetlabel}{Description of the DO-event.}
#'   \item{NGRIP_depth_m}{Onset of DO-event in NGRIP depth.}
#'   \item{GICC_age.yb2k}{Onset of DO-event in GICC05 age.}
#'   \item{depth_int_lower.m}{Start of data window expressed in terms of NGRIP depth.}
#'   \item{depth_int_upper.m}{End of data window expressed in terms of NGRIP depth.}
#'   \item{age.int.lower}{Start of data window expressed in terms of GICC05 age.}
#'   \item{age.int.upper}{End of data window expressed in terms of GICC05 age.}
#'   \item{index.int.lower}{Start of data window expressed in the number of (5 cm) steps from Holocene.}
#'   \item{index.int.upper}{End of data window expressed in the number of (5 cm) steps from Holocene.}
#'   }
#' @references Sune O. Rasmussen, Matthias Bigler, Simon P. Blockley, Thomas Blunier, Susanne L. Buchardt, Henrik B. Clausen, Ivana Cvijanovic, Dorthe Dahl-Jensen, Sigfus J. Johnsen, Hubertus Fischer, Vasileios Gkinis, Myriam Guillevic, Wim Z. Hoek, J. John Lowe, Joel B. Pedro, Trevor Popp, Inger K. Seierstad, Jørgen Peder Steffensen, Anders M. Svensson, Paul Vallelonga, Bo M. Vinther, Mike J.C. Walker, Joe J. Wheatley, Mai Winstrup, A stratigraphic framework for abrupt climatic changes during the Last Glacial period based on three synchronized Greenland ice-core records: refining and extending the INTIMATE event stratigraphy, Quaternary Science Reviews, Volume 106, 2014, Pages 14-28, https://doi.org/10.1016/j.quascirev.2014.09.007. (https://www.sciencedirect.com/science/article/pii/S0277379114003485)
#' @source \url{https://www.iceandclimate.nbi.ku.dk/data}

"DO_intervals"
