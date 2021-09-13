#' Optimal data window for DO-events
#'
#' Data windows enclosing DO-events found by optimization scheme to be used in linear ramp model fit.
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
#' @references Will be added soon!
#' @source \url{}

"DO_intervals"
