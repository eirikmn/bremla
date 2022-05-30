bremla_advanced <- function(formula,data,nsims=10000,events=NULL,
                            synchronization=NULL,
                            control.fit=NULL,
                            control.sim=NULL,
                            control.transition_dating=NULL,
                            control.bias=NULL
                            ){

  time.start = Sys.time()
  bremla.call = sys.call(which=1)
  if(print.progress) cat("Initiating data formatting...",sep="")

  #prepare bremla object by formatting dataset, writing formulastrings etc
  object = bremla_prepare(data, events=events, formula, control.fit=control.fit,
                          control.sim=control.sim,
                          control.transition_dating=control.transition_dating,
                          control.bias=control.bias)
  if(print.progress) cat(" completed!\n",sep="")

  if(!is.null(control.fit)){
    #fit the data, first by least squares, then by INLA (if specified)
    object = bremla_modelfitter(object, control.fit,print.progress=print.progress)

  }


  if(nsims>0 && control.sim$synchronized %in% c(FALSE,2)){
    #produce samples from the chronologies
    object = bremla_chronology_simulation(object, control.sim=control.sim,
                                          print.progress=print.progress)

    #compute posterior marginal mean, quantiles and other summary statistics
    #object = bremla_simulationsummarizer(object,CI.type=CI.type,sync=FALSE,print.progress=print.progress)

  }
  if(!is.null(synchronization)){
    ##format and or simulate tie-points
    object = tiepointsimmer(object, nsims=nsims, locations=synchronization$locations,
                            locations.type=synchronization$locations.type,
                            method=synchronization$method,samples=synchronization$samples)

  }


  if(nsims>0 && control.sim$synchronized %in% c(FALSE,2)){
    if(is.null(synchronization)){
      warning("Synchronization is currently only supported for AR(1) processes. Skipping this part...")
    }else{
      #produce samples from the chronologies
      object = bremla_synchronized_simulation(object, control.sim=control.sim,
                                              print.progress=print.progress)

      #compute posterior marginal mean, quantiles and other summary statistics
      #object = bremla_simulationsummarizer(object,CI.type=CI.type,sync=TRUE,print.progress=print.progress)
    }


  }


  #if event.estimation list object (containing specifications) is included, perform dating estimation
  if(!is.null(event.estimation)){
    #find onset depth posterior by fitting linear ramp model with INLA
    object = linrampfitter(object,control.transition_dating$linramp,
                           print.progress=print.progress)

    #perform Monte Carlo simulations to produce samples for onset age of warming transition
    object = events_depth_to_age(object, nsims = nsims, print.progress=print.progress,
                                 control.transition_dating$dating)
  }
  #if bias list object is included, perform this analysis
  if(!is.null(control.bias)){
    object = bremla_biased_chronologies(object,control.bias)
  }
  time.total = difftime(Sys.time(), time.start,units="secs")[[1]]
  object$.args$call = bremla.call
  object$time$total = time.total
  class(object) = "bremla"

}
