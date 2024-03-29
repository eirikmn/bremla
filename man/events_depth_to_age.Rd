% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/events_depth_to_age.R
\name{events_depth_to_age}
\alias{events_depth_to_age}
\title{Complete dating uncertainty of abrupt warming transitions}
\usage{
events_depth_to_age(object, control.transition_dating, print.progress = FALSE)
}
\arguments{
\item{object}{List object which is the output of function \code{\link{linrampfitter}}}

\item{control.transition_dating}{List containing speci}

\item{print.progress}{Boolean. If \code{TRUE} progress will be printed to screen}
}
\value{
Returns the same \code{object} list from the input, but appends the simulated transition onset ages along with summary statistics.
}
\description{
Combines linear onset posterior from ramp model fit and simulated chronologies from
the Bayesian regression modeling to estimate complete dating uncertainty of the onset
of abrupt warming transitions events using Monte Carlo simulation.
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
synchronization=list(locations=depth[c(100,400,700)],method="gauss",
                           params=list(mean=c(age[c(100,400,700)]+c(30,-100,50)),
                                       sd=c(50,20,100)
                                       )
                     )
control.sim=list(synchronized=TRUE,
                 summary=list(compute=TRUE))

#simulate transition:
prox = rnorm(n,mean=c(rep(0,400),seq(0,4,length.out=200),rep(4,400)),sd=1)
window = 330:500
control.linramp = list(label="Simulated",proxy=prox,interval=window,interval.unit="index",
    depth.ref=depth[401])
control.transition_dating=list(label="Simulated transition",dating=list(age.ref=age[401]))
object = bremla_prepare(formula,data,nsims=5000,reference.label="simulated timescale",
                        events = events,
                        synchronization=synchronization,
                        control.fit=control.fit,
                        control.sim=control.sim,
                        control.linramp=control.linramp,
                        control.transition_dating=control.transition_dating)
object = bremla_modelfitter(object)
object = tiepointsimmer(object)
object = bremla_synchronized_simulation(object)
object = linrampfitter(object)
object = events_depth_to_age(object,print.progress=TRUE)
summary(object)
plot(object)
}
}
}
\seealso{
\code{\link{bremla_chronology_simulation},\link{linrampfitter}}
}
\author{
Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
}
\keyword{bremla}
\keyword{dating}
\keyword{transition}
