% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpfulfunctions.R
\name{meanmaker}
\alias{meanmaker}
\title{Mean vector from fixed effects}
\usage{
meanmaker(coefs, reg.model, data)
}
\arguments{
\item{coefs}{Fixed effects}

\item{reg.model}{list specifying the structure of the linear regression model.}

\item{data}{data.frame obtained from \code{\link{bremla_prepare}}.}
}
\value{
Returns a vector representing the mean number of layer difference per depth interval.
}
\description{
Gives the mean number of layer differences per depth from the fixed effects.
}
\seealso{
\code{\link{bremla_prepare},\link{bremla_chronology_simulation}}
}
\author{
Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
}
\keyword{bremla}
\keyword{mean}
