#' Bayesian regression modeling of layered archives
#'
#' Fits a regression model to the data and produces chronology simulations.
#'
#' @param formula formula describing the predictor used in fitting the data. Partial covariates (psi)
#' can be filled in by specifying the degree and events in \code{events} argument.
#' @param data data.frame containing response and covariates. Must include 'age' and 'depth'.
#' The first row is used for initializing 'age' and 'depth', hence the rest can be set to NA.
#' Covariates related to the events (psi) can be filled in automatically by specifying the degree and events in \code{events} argument.
#' @param reference.label Character label of reference timescale. Used in \code{\link{plot},\link{summary}}.
#' @param x.label Character label for the x-axis (depth).
#' @param y.label Character label for the y-axis (age).
#' @param nsims Number of chronologies to be simulated.
#' @param events List object describing locations, degree and other informations about
#' the specific climatic periods used in the predictor. If used must include at least
#' \code{events\$locations}. See \code{\link{events.default}} for more details.
#' @param synchronization List containing locations and details for sampling tie-points used to synchronize our time scale.
#' @param control.fit List containing specifications for the fitting procedure.
#' See \code{\link{control.fit.default}} for more details.
#' @param control.sim List containing specifications for sampling chronologies (synchronous or otherwise).
#' See \code{\link{control.sim.default}} for more details.
#' @param control.linramp List containing specifications for fitting a linear ramp model to given data.
#' If used, must include \code{control.linramp\$proxy} and \code{control.linramp\$interval}. See
#' \code{\link{control.linramp.default}} for details.
#' @param control.transition_dating List containing specifications for
#' estimating a given transition. See \code{\link{control.transition_dating.default}} for more details.
#' @param control.bias List containing specifications for adding a potential unknown stochastic counting bias.
#' See \code{\link{control.bias.default}} for more details.
#' @param print.progress Boolean. If \code{TRUE} progress will be printed to screen
#'
#' @return Returns an S3 object of class "bremla". This includes output from all functions nested within the bremla function. Including fitted marginals and summary statistics, simulated chronologies, time spent on each step.
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{bremla_modelfitter},\link{bremla_chronology_simulation}}
#' @keywords bremla
#' @examples
#' \donttest{
#' if(inlaloader()){
#' ## Simulation example
#' require(stats)
#' set.seed(1)
#' n <- 1000
#' phi <- 0.8
#' sigma <- 1.2
#' a_lintrend <- 0.3; a_proxy = 0.8
#' dy_noise <- as.numeric(arima.sim(model=list(ar=c(phi)),n=n,sd=sqrt(1-phi^2)))
#' lintrend <- seq(from=10,to=15,length.out=n)
#'
#' proxy <- as.numeric(arima.sim(model=list(ar=c(0.9)),n=n,sd=sqrt(1-0.9^2)))*5
#' dy <- a_lintrend*lintrend + a_proxy*proxy + sigma*dy_noise
#'
#' y0 = 11700;z0=1200
#' age = y0+cumsum(dy)
#' depth = 1200 + 1:n*0.05
#'
#'
#' formula = dy~-1+depth2 + proxy
#' data = data.frame(age=age,dy=dy,proxy=proxy,depth=depth,depth2=depth^2)
#' data = rbind(c(y0,NA,NA,z0,NA),data) #First row is only used to extract y0 and z0.
#'
#' events=list(locations=c(1210,1220,1240))
#' control.fit = list(ncores=2,noise="ar1")
#' synchronization=list(locations=depth[c(100,400,700)],method="gauss",
#'                            params=list(mean=c(age[c(100,400,700)]+c(30,-100,50)),
#'                                        sd=c(50,20,100)
#'                                        )
#'                      )
#' control.sim=list(synchronized=2,
#'                  summary=list(compute=TRUE))
#'
#' #simulate transition:
#' prox = rnorm(n,mean=c(rep(0,400),seq(0,4,length.out=20),rep(4,580)),sd=1)-45
#' window = 330:500
#' control.linramp = list(label="Simulated",proxy=prox,interval=window,interval.unit="index",
#'     depth.ref=depth[401])
#' control.transition_dating=list(label="Simulated transition",dating=list(age.ref=age[401]))
#'
#'
#' object = bremla(formula,data,nsims=5000,reference.label="simulated timescale",
#'                          events=events,
#'                          synchronization=synchronization,
#'                          control.fit=control.fit,
#'                          control.sim=control.sim,
#'                          control.linramp=control.linramp,
#'                          control.transition_dating=control.transition_dating,
#'                          print.progress=TRUE)
#' summary(object)
#' plot(object)
#'
#' }
#' if(inlaloader()){
#' ### Real data example ###
#' require(stringr)
#' data("event_intervals")
#' data("events_rasmussen")
#' data("NGRIP_5cm")
#'
#'
#' age = NGRIP_5cm$age
#' depth = NGRIP_5cm$depth
#' d18O = NGRIP_5cm$d18O
#' proxy=d18O
#' data = data.frame(age=age,dy=c(NA,diff(age)),depth=depth,depth2=depth^2,proxy=proxy)
#' formula = dy~-1+depth2
#'
#'
#' lowerints = which.index(event_intervals$depth_int_lower.m, depth[2:length(depth)])
#' upperints = which.index(event_intervals$depth_int_upper.m, depth[2:length(depth)])
#'
#' eventnumber=13 #number between 1 and 29. specifies which transition to consider
#' transitionlabel = str_sub(event_intervals$onsetlabel[eventnumber],
#'                       str_locate(event_intervals$onsetlabel[eventnumber],"GI")[1])
#' depth.reference = event_intervals$NGRIP_depth_m[eventnumber]
#' age.reference = event_intervals$GICC_age.yb2k[eventnumber]
#' interval = lowerints[eventnumber]:upperints[eventnumber]
#'
#' nsims=5000
#' events=list(locations = events_rasmussen$depth,
#'             locations_unit="depth",degree=1)
#' synchronization = list(method="adolphi")
#' control.fit=list(method="inla")
#' control.sim=list(synchronized = 2)
#'
#' control.linramp=list(proxy=proxy,interval=interval,interval.unit="index",
#'                      depth.reference=depth.reference,
#'                      label=transitionlabel,depth.label="d18O (permil)")
#' control.transition_dating=list(age.ref=age.reference,age.label="years (yb2k)")
#' object = bremla(formula,data,reference.label="GICC05",
#'                    nsims=nsims,
#'                    events=events,
#'                    synchronization=synchronization,
#'                    control.fit=control.fit,
#'                    control.sim=control.sim,
#'                    control.linramp=control.linramp,
#'                    control.transition_dating=control.transition_dating,
#'                    print.progress=TRUE )
#'   summary(object)
#'   plot(object)
#'   }
#' }
#' @export
#' @import matrixStats
bremla <- function(formula,data,reference.label=NULL,
                            x.label=NULL,
                            y.label=NULL,
                            nsims=10000,
                            events=NULL,
                            synchronization=NULL,
                            control.fit=NULL,
                            control.sim=NULL,
                            control.linramp=NULL,
                            control.transition_dating=NULL,
                            control.bias=NULL,
                            print.progress=FALSE
){

  time.start = Sys.time()
  bremla.call = sys.call(which=1)
  if(print.progress) cat("Initiating data formatting...",sep="")


  #prepare bremla object by formatting dataset, writing formulastrings etc
  object = bremla_prepare(formula, data, ##data must include 'depth' and 'age' in data or input
                          reference.label=reference.label,
                          x.label = x.label,
                          y.label = y.label,
                          nsims=nsims,
                          events=events,
                          synchronization=synchronization,
                          control.fit=control.fit,
                          control.sim=control.sim,
                          control.linramp=control.linramp,
                          control.transition_dating=control.transition_dating,
                          control.bias=control.bias)
  if(print.progress) cat(" completed!\n",sep="")

  if(!is.null(control.fit)){
    #fit the data, first by least squares, then by INLA (if specified)

    if( .Platform$OS.type=="windows" && nrow(data)>7000){
      warning("Windows is poorly suited for large data sets. Try a different operating system if one is available.
              Using Windows Subsystem for Linux (WSL) might also provide a solution.")
      return(object)
    }

    object = bremla_modelfitter(object, control.fit, #set controls to NULL and use via object
                                print.progress=print.progress)

  }

  if(!is.null(control.sim)){
    if(nsims>0 && control.sim$synchronized %in% c(FALSE,2)){
      control.sim$nsims=nsims
      #produce samples from the chronologies
      object = bremla_chronology_simulation(object, control.sim=control.sim,
                                            print.progress=print.progress)

      #compute posterior marginal mean, quantiles and other summary statistics
      #object = bremla_simulationsummarizer(object,CI.type=CI.type,sync=FALSE,print.progress=print.progress)

    }
  }

  if(!is.null(synchronization)){
    synchronization$nsims=nsims
    ##format and or simulate tie-points
    object = tiepointsimmer(object, synchronization,
                            print.progress=print.progress)

  }

  if(!is.null(control.sim)){
    if(nsims>0 && control.sim$synchronized %in% c(TRUE,2)){
      control.sim$nsims=nsims
      #produce samples from the chronologies
      object = bremla_synchronized_simulation(object, control.sim=control.sim,
                                              print.progress=print.progress)

      #compute posterior marginal mean, quantiles and other summary statistics
      #object = bremla_simulationsummarizer(object,CI.type=CI.type,sync=TRUE,print.progress=print.progress)
    }
  }



  #if control.transition_dating list object (containing specifications) is included, perform dating estimation
  if(!is.null(control.linramp)){
    #find onset depth posterior by fitting linear ramp model with INLA
    object = linrampfitter(object,control.linramp,
                           print.progress=print.progress)

    if(!is.null(control.transition_dating)){
      control.transition_dating$dating$nsims = nsims
      #perform Monte Carlo simulations to produce samples for onset age of warming transition
      object = events_depth_to_age(object, control.transition_dating$dating,
                                   print.progress=print.progress)
    }

  }
  #if bias list object is included, perform this analysis
  if(!is.null(control.bias)){
    object = bremla_biased_chronologies(object,control.bias,
                                        print.progress=print.progress)
  }
  time.total = difftime(Sys.time(), time.start,units="secs")[[1]]
  object$.args$call = bremla.call
  object$time$total = time.total
  class(object) = "bremla"

  return(object)
}
