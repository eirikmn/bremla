% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bremla_chronology_simulation.R
\name{bremla_chronology_simulation}
\alias{bremla_chronology_simulation}
\title{Bremla chronology simulation}
\usage{
bremla_chronology_simulation(object, control.sim, print.progress = FALSE)
}
\arguments{
\item{object}{List object which is the output of function \code{\link{bremla_modelfitter}}}

\item{control.sim}{List object containing specifications for simulation procedure and
what is to be computed. See \code{\link{control.sim.default}} for details.}

\item{print.progress}{Boolean. If \code{TRUE} progress will be printed to screen}
}
\value{
Returns the same \code{object} list from the input, but appends simulated
chronologies, mean vectors (if \code{store.means=TRUE}) and further information, summary and statistics in
\code{object\$simulation}.
}
\description{
Simulates chronologies based on the fitted regression model.
}
\examples{
\donttest{
if(inlaloader()){
require(stats)
set.seed(1)
n <- 1000
phi <- 0.8
sigma <- 1.2
a_lintrend <- 0.3; a_proxy = 0.8
dy_noise <- as.numeric(arima.sim(model=list(ar=c(phi)),n=n,sd=sqrt(1-phi^2)))
lintrend <- seq(from=10,to=15,length.out=n)

proxy <- as.numeric(arima.sim(model=list(ar=c(0.9)),n=n,sd=sqrt(1-0.9^2)))
dy <- a_lintrend*lintrend + a_proxy*proxy + sigma*dy_noise

y0 = 11700;z0=1200
age = y0+cumsum(dy)
depth = 1200 + 1:n*0.05
depth2 = depth^2/depth[1]^2 #normalize for stability


formula = dy~-1+depth2 + proxy
data = data.frame(age=age,dy=dy,proxy=proxy,depth=depth,depth2=depth2)
data = rbind(c(y0,NA,NA,z0,NA),data) #First row is only used to extract y0 and z0.

events=list(locations=c(1210,1220,1240))
control.fit = list(ncores=2,noise="ar1")
control.sim=list(synchronized=2,
                 summary=list(compute=TRUE))

object = bremla_prepare(formula,data,nsims=5000,reference.label="simulated timescale",
                        events = events,
                        control.fit=control.fit,
                        control.sim=control.sim)
object = bremla_modelfitter(object)
object = bremla_chronology_simulation(object, print.progress=TRUE)
summary(object)
plot(object)
}
}

}
\seealso{
\code{\link{bremla_prepare},\link{bremla_modelfitter}}
}
\author{
Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
}
\keyword{bremla}
\keyword{simulation}
