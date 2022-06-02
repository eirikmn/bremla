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
#'
#' @examples
#' \donttest{
#' data("event_intervals")
#' data("events_rasmussen")
#' data("NGRIP_5cm")
#'
#' age = NGRIP_5cm$age
#' depth = NGRIP_5cm$depth
#' d18O = NGRIP_5cm$d18O
#' proxy=d18O
#'
#' eventdepths = events_rasmussen$depth
#' eventindexes = c(1,which.index(eventdepths, depth[2:length(depth)]) )
#' eventindexes = unique(eventindexes[!is.na(eventindexes)])
#'
#' nsims = 5000
#' object = bremla_prepare(age,depth,proxy,events=eventdepths,nsims=nsims)
#' object = bremla_modelfitter(object)
#' object = bremla_chronology_simulation(object,nsims=nsims)
#' object = tiepointsimmer(object,nsims=nsims,method="adolphi")
#' object = bremla_synchronized_simulation(object,nsims=nsims)
#' object = bremla_simulationsummarizer(object,CI.type="hpd",sync=TRUE)
#'
#' lowerints = which.index(event_intervals$depth_int_lower.m, depth[2:length(depth)])
#' upperints = which.index(event_intervals$depth_int_upper.m, depth[2:length(depth)])
#'
#' eventnumber=13 #number between 1 and 29. specifies which transition to consider
#' depth.reference = event_intervals$NGRIP_depth_m[eventnumber]
#' age.reference = event_intervals$GICC_age.yb2k[eventnumber]
#' interval = lowerints[eventnumber]:upperints[eventnumber]
#'
#' object = linrampfitter(object,interval,label="GI-11",depth.reference=depth.reference)
#' object = events_depth_to_age(object,nsims=nsims,age.reference=age.reference)
#' plot(object,plot.proxydata=list(age=TRUE,depth=FALSE,xrev=FALSE,label=NULL),
#'     plot.ls = NULL,
#'     plot.inla.posterior = NULL,
#'     plot.inlasims = NULL,
#'     plot.syncsims = NULL,
#'     plot.tiepoints = NULL,
#'     plot.bias = NULL,
#'     plot.linramp = list(depth.reference=NULL,show.t0=TRUE,show.t1=TRUE,xrev=TRUE,label=NULL))
#' }
#'
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


