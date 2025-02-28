% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/default_arguments.R
\name{control.sim.default}
\alias{control.sim.default}
\title{Default variables in control.sim}
\usage{
control.sim.default()
}
\value{
Returns a list including default values for all variables in \code{control.sim}.
}
\description{
Sets the default variables in the list \code{control.sim} used to specify the
fitting procedure. The list contains the following arguments:
\describe{
  \item{\code{synchronized} }{Boolean indicating whether or not synchronization is to be performed.
  Default value is \code{TRUE}, but this is set to \code{FALSE} if supplied tie-points are not found.}
  \item{\code{nsims} }{Integer giving how many simulations should be generated. Default value is \code{10000}.}
  \item{\code{ncores} }{Integer which sets the number of cores to be used to sample from
  the joint posterior obtained by INLA. Default value is \code{2}.}
  \item{\code{store.everything} }{Boolean indicating whether or not to store superfluous information
  such as the simulated mean vectors (from fixed effects). Default \code{FALSE}.}
  \item{\code{summary} }{List object with two items
    \describe{
    \item{\code{compute} }{Boolean indicating if summary statistics (posterior marginal means median
    and credible intervals) should be computed. Default \code{TRUE}.}
    \item{\code{CI.type} }{character giving which type of credible intervals to be computed.
    \code{"quantiles"} for quantiles based intervals (Default) or \code{"hpd"} for
    highest posterior density intervals.}}
    }
}
}
\seealso{
\code{\link{bremla},\link{bremla_chronology_simulation},
\link{bremla_synchronized_simulation},\link{set.options}}
}
\author{
Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
}
\keyword{bremla}
\keyword{default}
\keyword{simulation}
