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
#' @param reg.model List of booleans specifies which effects to be included in the linear regression model: intercept (\code{const}), linear response to depth (\code{depth1}), quadratic response to depth (\code{depth2}), linear response to proxy (\code{proxy}), individual intercepts for each climate period (\code{psi0}), individual linear responses to depth for each climate period (\code{psi1})
#' @param noise Character specifying the noise model: independent identically distributed (\code{iid}), first order autoregressive (\code{ar1}) or second order autoregressive (\code{ar2})
#' @param method Character specifying how the model is fitted. Currently only least squares (\code{ls}) and INLA (\code{inla}) are supported, and least squares is always run.
#' @param DO.estimation List specifying the settings for the linear ramp fit and subsequent age simulation of a specified DO-event. See example.
#' @param bias List specifying the settings for the simulation of chronologies subjected to unknown counting bias. See example for implementation.
#' @param store.everything Boolean describing whether optional information (e.g. simulation of age vectors) should be stored.
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
bremla = function(age,depth,proxy, events=NULL,nsims=10000, eventmeasure = "depth",reg.model = list(
  const=FALSE,depth1=FALSE,depth2=TRUE,proxy=TRUE,psi0=TRUE,psi1=TRUE), noise="ar1", method="inla",
  DO.estimation = NULL, bias = NULL,store.everything=FALSE,print.progress=FALSE){
  #DO.estimation = list(interval=...,h=0.1,t1.sims=50000,rampsims=50000,label="GI-11,depth.reference=NULL,age.reference=NULL)
  #bias = list(bias.model="uniform",biasparams=cbind(c(1,1),c(0.98,1.02),c(0.96,1.04)),store.samples=FALSE)

  time.start = Sys.time()
  if(print.progress) cat("Initiating data formatting...",sep="")

  object = bremla_prepare(age,depth,proxy, events=events,nsims=nsims, eventmeasure = eventmeasure,reg.model = reg.model,
                         noise=noise, method=method)
  if(print.progress) cat(" completed!\n",sep="")
  object = bremla_modelfitter(object, method=method,print.progress=print.progress)

  object = bremla_chronology_simulation(object, nsims=nsims, method=method,store.means=store.everything,print.progress=print.progress)
  object = bremla_simulationsummarizer(object,print.progress=print.progress)
  if(!is.null(DO.estimation)){
    object = linrampfitter(object,interval=DO.estimation$interval,h=DO.estimation$h,t1.sims=DO.estimation$t1.sims,
               rampsims=DO.estimation$rampsims,label=DO.estimation$label,print.progress=print.progress)

    object = DO_depth_to_age(object, nsims = nsims, print.progress=print.progress,label=DO.estimation$label)
  }
  time.total = Sys.time() - time.start
  object$time = list(total=time.total)
  class(object) = "bremla"

  return(object)
}
