% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/default_arguments.R
\name{synchronization.default}
\alias{synchronization.default}
\title{Default variables in synchronization}
\usage{
synchronization.default()
}
\value{
Returns a list including default values for all variables in
\code{synchronization}.
}
\description{
Sets the default variables in the list \code{synchronization} used to specify the
tie-points used for synchronization. The list contains the following arguments:
\itemize{
  \item \code{locations} Numeric describing the locations of the tie-points. Default value
  is here \code{NULL} meaning that this has to be given by the user. \item \code{locations_unit} Character giving which axis the \code{locations} are given.
  Can be \code{"depth"} (default), \code{"age"} or \code{"index"}.
  \item \code{samples} data.frame giving pre-computed samples. This should be used if you want to provide the sampled tie-points yourself. \code{nrow} must correspond to
  the number of samples and \code{ncol} to the number of tie-points. Default value is \code{NULL}.
  \item \code{method} Character giving which method should be used to generate tie-points. Currently supported: \code{"adolphi"} (default) uses tie-point distributions given by Adolphi et al. (2018) and Muscheler et al. (2020), \code{"normal"} uses a normal distribution with mean and standard deviations supplied by user (TODO). If \code{"adolphi"} is used, then \code{locations}, \code{locations_unit} is filled out accordingly.
  \item \code{x.ref} Numeric giving a reference for the x-axis (given by \code{locations_unit}) that will be included in future \code{plot}.
  \item \code{params} List object containing parameters for parametric tie-point distributions. for \code{method="gauss"} the list has two numeric objects called \code{mean} and \code{sd} representing the mean and standard deviation for each tie point, respectively.
  \item \code{agedisc} List containing arguments relating to the age-discrepancy model. If \code{NULL}, no age-discrepancy model is assigned and the tie-points are treated as realizations of the posterior predictive of the layer model.
  \describe{
  \item{\code{model}}{Age-discrepancy model. Currently, only 2nd order random walk model (\code{rw2}) is supported.}
  \item{\code{method}}{How the age-discrepancy model should be fitted. Currently, only \code{INLA} is supported.}
  \item{\code{options}}{A set of options passed into the fitting procedure. For \code{method=INLA}, this sets the additional arguments for the \code{inla} function.}
  \item{\code{hyperprior}}{List corresponding to the \code{hyper} object in INLA's \code{f} function which specifies the priors of the hyperparameters. Default is \code{NULL} which sets INLA's default priors.}
  \item{\code{stepsizes}}{If INLA struggles to converge, attempt again with different stepsizes according to \code{stepsizes}. default value is 0.005.}
  }
}
}
\seealso{
\code{\link{bremla}, \link{tiepointsimmer},
\link{bremla_synchronized_simulation},\link{set.options}}
}
\author{
Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
}
\keyword{bremla}
\keyword{default}
\keyword{simulation}
