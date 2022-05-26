#' Complete dating uncertainty of abrupt warming transitions
#'
#' Combines linear ramp model fit and Bayesian regression modeling to estimate complete dating uncertainty of the onset of Dansgaard-Oeschger events using Monte Carlo simulation.
#'
#' @param object List object which is the output of function \code{\link{linrampfitter}}
#' @param nsims Integer denoting the number of simulations of transition onset age. Cannot exceed the number of simulated age chronologies from \code{\link{bremla_chronology_simulation}}.
#' @param print.progress Boolean. If \code{TRUE} progress will be printed to screen
#' @param label character string describing the label designed the transition, e.g. \code{label="GI-11"}.
#' @param age.reference numeric denoting a reference value such as the transition onset age obtained by other means. This appears when using \code{plot}. If \code{depth.reference=NULL} then this is not used.
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
events_depth_to_age = function(object, nsims = 10000, print.progress=FALSE, label=NULL,age.reference=NULL){
  time.start=Sys.time()
  n=length(object$data$y)
  if(is.null(label)) label=object$linramp$.args$label
  if(print.progress){
    if(!is.null(label) && label !=""){
      cat("Initiating simulation from complete dating uncertainty of event: ",label,"\n",sep="")
    }else{
      cat("Initiating simulation from complete dating uncertainty of unlabeled event...\n",sep="")
    }
  }


  if(nsims > dim(object$simulation$age)[2]) stop(paste0("Not enough samples from age-depth model (",dim(object$simulation$age)[2],", ",nsims, "needed)"))
  z.sample = inla.rmarginal(nsims,object$linramp$param$t0$marg.t0)

  if(is.null(object$simulation)){
    warning("Cannot find simulations from age-depth model in 'agedepthmodel'!")
  }else{

    time.startsim = Sys.time()
    ysamps=numeric(nsims)
    z = object$data$z
    for(i in 1:nsims){
      if(print.progress==TRUE && i%%1000==0){
        cat("Onset age simulation ",i,"/",nsims,"\n",sep="")
      }
      zgrid = sort(c(z.sample[i],z))
      sampleindex = which(zgrid==z.sample[i] )
      distance = (z.sample[i]-zgrid[sampleindex-1])/(zgrid[sampleindex+1]-zgrid[sampleindex-1])
      ysamps[i] = (1-distance)*object$simulation$age[sampleindex-1,i] +distance*object$simulation$age[sampleindex,i]
    }
  }
  time.endsim = Sys.time()
  if(print.progress) cat("Completed in ", difftime(Sys.time(),time.startsim,units="secs")[[1]]," seconds!\n",sep="")



  dens = density(ysamps); Y0marg= cbind(dens$x,dens$y); colnames(Y0marg) = c("x","y"); Y0marg = inla.smarginal(Y0marg)
  z.Y0 = inla.zmarginal(Y0marg,silent=TRUE)

  object$event_dating = list(samples = ysamps, mean = z.Y0$mean,sd=z.Y0$sd,q0.025=z.Y0$quant0.025,q0.5=z.Y0$quant0.5,q0.975=z.Y0$quant0.975)

  object$event_dating$.args = list(nsims = nsims,label=label,age.reference=age.reference)
  object$time$event_age = list(simulation = difftime(time.endsim,time.startsim,units="secs")[[1]],total = difftime(Sys.time(),time.start,units="secs")[[1]])

  return(object)
}


