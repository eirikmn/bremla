#' Complete dating uncertainty of abrupt warming transitions
#'
#' Combines linear onset posterior from ramp model fit and simulated chronologies from
#' the Bayesian regression modeling to estimate complete dating uncertainty of the onset
#' of abrupt warming transitions events using Monte Carlo simulation.
#'
#' @param object List object which is the output of function \code{\link{linrampfitter}}
#' @param control.transition_dating List containing speci
#' @param print.progress Boolean. If \code{TRUE} progress will be printed to screen
#'
#' @return Returns the same \code{object} list from the input, but appends the simulated transition onset ages along with summary statistics.
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{bremla_chronology_simulation},\link{linrampfitter}}
#' @keywords bremla dating transition
#' @examples
#' \donttest{
#' require(stats)
#' set.seed(1)
#' n <- 1000
#' phi <- 0.8
#' sigma <- 1.2
#' a_lintrend <- 0.3; a_proxy = 0.8
#' dy_noise <- as.numeric(arima.sim(model=list(ar=c(phi)),n=n,sd=sqrt(1-phi^2)))
#' lintrend <- seq(from=10,to=15,length.out=n)
#'
#' proxy <- as.numeric(arima.sim(model=list(ar=c(0.9)),n=n,sd=sqrt(1-0.9^2)))
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
#' control.sim=list(synchronized=TRUE,
#'                  summary=list(compute=TRUE))
#'
#' #simulate transition:
#' prox = rnorm(n,mean=c(rep(0,400),seq(0,4,length.out=20),rep(4,580)),sd=1)
#' window = 330:500
#' control.linramp = list(label="Simulated",proxy=prox,interval=window,interval.unit="index",
#'     depth.ref=depth[401])
#' control.transition_dating=list(label="Simulated transition",dating=list(age.ref=age[401]))
#' object = bremla_prepare(formula,data,nsims=5000,reference.label="simulated timescale",
#'                         events = events,
#'                         synchronization=synchronization,
#'                         control.fit=control.fit,
#'                         control.sim=control.sim,
#'                         control.linramp=control.linramp,
#'                         control.transition_dating=control.transition_dating)
#' object = bremla_modelfitter(object)
#' object = tiepointsimmer(object)
#' object = bremla_synchronized_simulation(object)
#' object = linrampfitter(object)
#' object = events_depth_to_age(object,print.progress=TRUE)
#' summary(object)
#' plot(object)
#'
#' }
#' @export
#' @import INLA
#' @importFrom INLA inla.rgeneric.define inla.emarginal inla.smarginal
#' @importFrom stats optim
#'
events_depth_to_age = function(object, control.transition_dating,print.progress=FALSE){

#nsims = 10000, print.progress=FALSE, label=NULL,age.reference=NULL){
  time.start=Sys.time()


  if(missing(control.transition_dating)){

    if(!is.null(object$.args$control.transition_dating)){
      if(print.progress){
        cat("'control.transition_dating' missing. Importing information from 'object'.",sep="")
      }
      control.transition_dating = object$.args$control.transition_dating
    }else{
      stop("Could not find 'control.transition_dating'. Stopping.")
    }
  }
  #if(!is.null(control.transition_dating))
  control.transition_dating = set.options(control.transition_dating,control.transition_dating.default())

  object$.args$control.transition_dating = control.transition_dating

  ## sample hyperparameters
  if(is.null(object$linramp))
    stop("Linear ramp fit not found. Run 'linrampfitter' first!")


  n=nrow(object$data)

  if(!is.null(control.transition_dating$label) && control.transition_dating$label!=""){
    label=control.transition_dating$label
  }else{
    if(!is.null(object$.args$control.linramp$label) && object$.args$control.linramp$label!=""){
      label=object$.args$control.linramp$label
    }else{
      label=""
    }
  }


  if(print.progress){
    if(!is.null(label) && label !=""){

      cat("Initiating simulation from complete dating uncertainty of event: ",label,"\n",sep="")
    }else{
      cat("Initiating simulation from complete dating uncertainty of unlabeled event...\n",sep="")
    }
  }

  nsims = min(control.transition_dating$nsims,object$.args$control.sim$nsims,
              object$.args$synchronization$nsims)
  if(nsims > ncol(object$simulation$age))
    stop(paste0("Not enough samples from age-depth model (",ncol(object$simulation$age),", ",nsims, "needed)"))
  z.sample = inla.rmarginal(nsims,object$linramp$param$t0$marg.t0)

  if(is.null(object$simulation)){
    warning("Cannot find simulations from age-depth model in 'agedepthmodel'!")
  }else{

    if(control.transition_dating$sync){
      if(is.null(object$simulation$age_sync)){
        warning("Could not find synchronized chronology samples. Using unsynchronized instead.")
        samples = object$simulation$age
      }else{
        samples = object$simulation$age_sync
      }
    }else{
      samples = object$simulation$age
    }

    time.startsim = Sys.time()
    ysamps=numeric(nsims)
    z = object$data$depth
    for(i in 1:nsims){
      if(print.progress==TRUE && i%%1000==0){
        cat("Onset age simulation ",i,"/",nsims,"\n",sep="")
      }
      zgrid = sort(c(z.sample[i],z))
      sampleindex = which(zgrid==z.sample[i] )
      distance = (z.sample[i]-zgrid[sampleindex-1])/(zgrid[sampleindex+1]-zgrid[sampleindex-1])

      ysamps[i] = (1-distance)*samples[sampleindex-1,i] +distance*samples[sampleindex,i]
    }
  }
  time.endsim = Sys.time()
  if(print.progress) cat("Completed in ", difftime(Sys.time(),time.startsim,units="secs")[[1]]," seconds!\n",sep="")



  dens = density(ysamps); Y0marg= cbind(dens$x,dens$y); colnames(Y0marg) = c("x","y"); Y0marg = inla.smarginal(Y0marg)
  z.Y0 = inla.zmarginal(Y0marg,silent=TRUE)

  object$event_dating = list(samples = ysamps, mean = z.Y0$mean,sd=z.Y0$sd,q0.025=z.Y0$quant0.025,q0.5=z.Y0$quant0.5,q0.975=z.Y0$quant0.975)

  object$event_dating$.args = list(nsims = nsims,
                                   label=label,
                                   age.reference=control.transition_dating$age.ref)
  object$time$event_age = list(simulation = difftime(time.endsim,time.startsim,units="secs")[[1]],total = difftime(Sys.time(),time.start,units="secs")[[1]])

  return(object)
}


