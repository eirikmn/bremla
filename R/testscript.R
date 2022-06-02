if(FALSE){


source("r/bremla_advanced.R")
source("r/bremla_prepare.R")
source("r/bremla_modelfitter.R")
source("r/helpfulfunctions.R")
source("r/default_arguments.R")
source("r/bremla_chronology_simulation.R")
source("r/bremla_synchronized_simulation.R")
source("r/linrampfitter.R")
source("r/events_depth_to_age.R")
source("r/plot.bremla.R")
source("r/summary.bremla.R")
source("r/print.bremla.R")
source("r/rgeneric.uneven.AR1.R")
library(matrixStats)
require(stats)
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


formula = dy~-1+depth2 + proxy
data = data.frame(age=age,dy=dy,proxy=proxy,depth=depth,depth2=depth^2)
dd = rbind(c(y0,NA,NA,z0,NA),data)

events=list(locations=c(1210,1220,1240))
control.fit = list(ncores=2,noise="ar1")
synchronization=list(method="gauss")
control.sim=list(synchronized=2,
                 summary=list(compute=TRUE))
control.bias=NULL

#simulate transition:
prox = rnorm(n,mean=c(rep(0,400),seq(0,4,length.out=20),rep(4,580)),sd=1)
window = 330:500
control.transition_dating=list(label="Simulated transition",
                               linramp=list(proxy=prox,
                                            interval=window,
                                            interval.what="index",
                                            depth.ref=depth[400]),
                               dating=list(age.ref=age[400]))

print.progress=TRUE
stop("stopper her")

object = bremla_advanced(formula,dd,nsims=5000,events=events,
                         synchronization=synchronization,
                         control.fit=control.fit,
                         control.sim=control.sim,
                         control.transition_dating=control.transition_dating,
                         print.progress=TRUE)
summary(object)
plot(object)





object = bremla_prepare(formula, data=dd, ##data must include 'depth' and 'age' in data or input
                        events=events,
                        synchronization=synchronization,
                        control.fit=control.fit,
                        control.sim=control.sim,
                        control.transition_dating=control.transition_dating,
                        control.bias=control.bias)

## se over at initialverdiene som ble funnet i lm passer til inla
object = bremla_modelfitter(object, #set controls to NULL and use via object
                            print.progress=print.progress)

 object = bremla_chronology_simulation(object,print.progress=print.progress)
#
 synchronization=list(method="gauss",nsims=object$.args$control.sim$nsims)
 object = tiepointsimmer(object,synchronization,print.progress=print.progress)
#
 control.sim$synchronized=TRUE
 object = bremla_synchronized_simulation(object,control.sim,print.progress=TRUE)
#


prox = rnorm(n,mean=c(rep(0,400),seq(0,4,length.out=20),rep(4,580)),sd=1)
window = 330:500
plot(prox[window])
control.transition_dating = list(linramp=list(interval=window,proxy=prox,
                                              interval.what="index"))
library(numDeriv)
object = linrampfitter(object,control.transition_dating,
                       print.progress=print.progress)


 #perform Monte Carlo simulations to produce samples for onset age of warming transition
 object = events_depth_to_age(object,control.transition_dating$dating,
                              print.progress=print.progress)


control.bias = list(bias.model="uniform")
object = bremla_biased_chronologies(object,control.bias,print.progress=TRUE)
stop("stop")

plot(object)
summary(object)
print(object)


}
