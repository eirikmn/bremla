#' Complete dating uncertainty of DO-onsets
#'
#' Combines linear ramp model fit and Bayesian regression modeling to estimate complete dating uncertainty of the onset of Dansgaard-Oeschger events using Monte Carlo simulation.
#'
#' @param object List object which is the output of function \code{\link{linrampfitter}}
#' @param nsims Integer denoting the number of simulations of DO onset age. Cannot exceed the number of simulated age chronologies from \code{\link{bremla_chronology_simulation}}.
#' @param print.progress Boolean. If \code{TRUE} progress will be printed to screen
#' @param label character string describing the label designed the transition, e.g. \code{label="GI-11"}.
#' @param age.reference numeric denoting a reference value such as the transition onset age obtained by other means. This appears when using \code{plot}. If \code{depth.reference=NULL} then this is not used.
#'
#' @return Returns the same \code{object} list from the input, but appends the simulated DO onset ages along with summary statistics.
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{bremla_chronology_simulation},\link{linrampfitter}}
#' @keywords bremla dating DO
#'
#' @examples
#'
#' @export
#' @import INLA
#' @importFrom INLA inla.rgeneric.define inla.emarginal inla.smarginal
#' @importFrom stats optim
#'
DO_depth_to_age = function(object, nsims = 10000, print.progress=FALSE, label=NULL,age.reference=NULL){
  time.start=Sys.time()
  n=length(object$data$y)
  if(is.null(label)) label=object$linramp$.args$label
  if(print.progress){
    if(!is.null(label) && label !=""){
      cat("Initiating simulation from complete dating uncertainty of DO-event: ",label,"\n",sep="")
    }else{
      cat("Initiating simulation from complete dating uncertainty of unlabeled DO-event...\n",sep="")
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

    }
  }
  time.endsim = Sys.time()
  if(print.progress) cat("Completed in ", difftime(Sys.time(),time.startsim,units="secs")[[1]]," seconds!\n",sep="")



  dens = density(ysamps); Y0marg= cbind(dens$x,dens$y); colnames(Y0marg) = c("x","y"); Y0marg = inla.smarginal(Y0marg)
  z.Y0 = inla.zmarginal(Y0marg,silent=TRUE)

  object$DO_dating = list(samples = ysamps, mean = z.Y0$mean,sd=z.Y0$sd,q0.025=z.Y0$quant0.025,q0.5=z.Y0$quant0.5,q0.975=z.Y0$quant0.975)

  object$DO_dating$.args = list(nsims = nsims,label=label,age.reference=age.reference)
  object$time$DO_age = list(simulation = difftime(time.endsim,time.startsim,units="secs")[[1]],total = difftime(Sys.time(),time.start,units="secs")[[1]])

  return(object)
}


