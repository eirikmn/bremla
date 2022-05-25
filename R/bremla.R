#' Bayesian regression modeling of layered archives
#'
#' Fits a regression model to the data and produces chronology simulations.
#'
#' @param age Vector of observed ages y_0,...,y_n
#' @param depth Vector of depths z_0,...,z_n corresponding to \code{age}
#' @param proxy Vector of observed proxy (e.g. del-18O) at depths \code{depth}
#' @param events Vector describing locations (in dimension specified by \code{eventmeasure}) of climate transitions for use in linear regression model
#' @param nsims Number of chronologies to be simulated
#' @param eventmeasure Character describing in which dimension the climate events are located ("age", "depth", "index")
#' @param reference.label String denoting the label used for age data/reference
#' @param reg.model List of booleans specifies which effects to be included in the linear regression model: intercept (\code{const}), linear response to depth (\code{depth1}), quadratic response to depth (\code{depth2}), linear response to proxy (\code{proxy}), individual intercepts for each climate period (\code{psi0}), individual linear responses to depth for each climate period (\code{psi1})
#' @param proxy.type String denoting which proxy is used (for plot-labels).
#' @param transform String denoting which transformation should be applied to the residuals. Default "identity", "log" is also supported, but experimental.
#' @param noise Character specifying the noise model: independent identically distributed (\code{iid}), first order autoregressive (\code{ar1}) or second order autoregressive (\code{ar2})
#' @param method Character specifying how the model is fitted. Currently only least squares (\code{ls}) and INLA (\code{inla}) are supported, and least squares is always run.
#' @param event.estimation List specifying the settings for the linear ramp fit and subsequent age simulation of a specified transition. See example.
#' @param bias List specifying the settings for the simulation of chronologies subjected to unknown counting bias. See example for implementation.
#' @param store.everything Boolean describing whether optional information (e.g. simulation of age vectors) should be stored.
#' @param CI.type Character describing which type of uncertainty intervals to be used. \code{CI.type="quantiles"} is computed from the quantiles, \code{CI.type="hpd"} computes highest posterior density intervals (takes a little longer to compute). If distribution is unimodal and symmetric these are equivalent.
#' @param synchronization List containing information used for simulating tie-points and synchronizing age-depth model to them.
#' @param print.progress Boolean. If \code{TRUE} progress will be printed to screen
#'
#' @return Returns an S3 object of class "bremla". This includes output from all functions nested within the bremla function. Including fitted marginals and summary statistics, simulated chronologies, time spent on each step.
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{bremla_modelfitter},\link{bremla_chronology_simulation}}
#' @keywords bremla
#'
#' @examples
#'
#' @export
#' @import matrixStats
bremla = function(age,depth,proxy, events=NULL,nsims=10000, eventmeasure = "depth",
                  reference.label=NULL, proxy.type="d18O",transform="identity",
                  reg.model = list(
                    const=FALSE,depth1=FALSE,depth2=TRUE,proxy=TRUE,psi0=TRUE,psi1=TRUE),
  noise="ar1", method="inla", CI.type="quantiles",
  synchronization = NULL,#list(locations=c(11050,12050,13050,22050,42050),locations.type="age",method="adolphi",samples=NULL),
  event.estimation = NULL, bias = NULL,store.everything=FALSE,print.progress=FALSE){

  time.start = Sys.time()
  bremla.call = sys.call(which=1)
  if(print.progress) cat("Initiating data formatting...",sep="")

  #prepare bremla object by formatting dataset, writing formulastrings etc
  object = bremla_prepare(age,depth,proxy, events=events,nsims=nsims, eventmeasure = eventmeasure,
                          reg.model = reg.model,noise=noise, method=method,
                          reference.label=reference.label, transform=transform,
                          proxy.type=proxy.type)
  if(print.progress) cat(" completed!\n",sep="")

  #fit the data, first by least squares, then by INLA (if specified)
  object = bremla_modelfitter(object, method=method,print.progress=print.progress)

  #produce samples from the chronologies
  object = bremla_chronology_simulation(object, nsims=nsims, method=method,store.means=store.everything,print.progress=print.progress)

  #compute posterior marginal mean, quantiles and other summary statistics
  object = bremla_simulationsummarizer(object,CI.type=CI.type,sync=FALSE,print.progress=print.progress)

  if(!is.null(synchronization)){
    if(tolower(noise) %in% c("iid","ar2","2","ar(2)") ){
      warning("Synchronization is currently only supported for AR(1) processes. Skipping this part...")
    }else{
      ##format and or simulate tie-points
      object = tiepointsimmer(object, nsims=nsims, locations=synchronization$locations,
                              locations.type=synchronization$locations.type,
                              method=synchronization$method,samples=synchronization$samples)
      ##simulate synchronized chronologies
      object = bremla_synchronized_simulation(object,nsims=nsims) #rest of input arguments are collected from 'object'
    }
  }


  #compute posterior marginal mean, quantiles and other summary statistics for synchronous time scale
  object = bremla_simulationsummarizer(object,CI.type=CI.type,sync=TRUE,print.progress=print.progress)

  #if event.estimation list object (containing specifications) is included, perform dating estimation
  if(!is.null(event.estimation)){
    #find onset depth posterior by fitting linear ramp model with INLA
    object = linrampfitter(object,interval=event.estimation$interval,
                           h=event.estimation$h,t1.sims=event.estimation$t1.sims,
                           rampsims=event.estimation$rampsims,
                           label=event.estimation$label,
                           depth.reference = event.estimation$depth.reference,
                           print.progress=print.progress)

    #perform Monte Carlo simulations to produce samples for onset age of warming transition
    object = events_depth_to_age(object, nsims = nsims, print.progress=print.progress,
                                label=event.estimation$label,
                                age.reference = event.estimation$age.reference)
  }
  #if bias list object is included, perform this analysis
  if(!is.null(bias)){
    object = bremla_biased_chronologies(object,bias.model=bias$bias.model,
                                        biasparams = bias$biasparams,nsims=nsims,
                                        store.samples=store.everything)
  }
  time.total = difftime(Sys.time(), time.start,units="secs")[[1]]
  object$.args$call = bremla.call
  object$time$total = time.total
  class(object) = "bremla"

  return(object)
}
