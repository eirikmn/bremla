#' Mean vector from fixed effects
#'
#' Gives the mean number of layer differences per depth from the fixed effects.
#'
#' @param coefs Fixed effects
#' @param reg.model list specifying the structure of the linear regression model.
#' @param nevents The number of climate transitions
#' @param data data.frame obtained from \code{\link{bremla_prepare}}.
#'
#' @return Returns a vector representing the mean number of layer difference per depth interval.
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{bremla_prepare},\link{bremla_chronology_simulation}}
#' @keywords bremla mean
#'
#' @examples
#' @export
meanmaker = function(coefs,reg.model,nevents=69,data){
  ## Computes mean vector from given fixed effects 'coefs'.
  ## Requires specification of which effects to include ('reg.model'), the number of climate transitions ('nevents') and a data.frame with covariates ('data')
  coefcounter=1
  fitted=numeric(dim(data)[1])
  if(reg.model$const){
    fitted=coefs[1]
    coefcounter=coefcounter+1
  }
  if(reg.model$depth1){
    fitted = fitted + coefs[coefcounter]*data$z
    coefcounter=coefcounter+1
  }
  if(reg.model$depth2){
    fitted = fitted + coefs[coefcounter]*data$z2
    coefcounter=coefcounter+1
  }
  if(reg.model$proxy){
    fitted=fitted + coefs[coefcounter]*data$x
    coefcounter=coefcounter+1
  }
  if(nevents>0){
    for(i in 2:nevents){
      if(reg.model$psi0){
        #aa = coefs[coefcounter]
        fitted = fitted + coefs[coefcounter]*data[[paste0("a",i-1)]]
        coefcounter=coefcounter+1
      }
      if(reg.model$psi1){
        #cc = coefs[coefcounter]
        fitted = fitted + coefs[coefcounter]*data[[paste0("c",i-1)]]
        coefcounter=coefcounter+1
      }
    }
  }
  return(fitted)
}
#' Linear ramp function
#'
#' Computes the (noiseless) linear ramp function.
#'
#' @param t Vector describing the x-axis (here, depth)
#' @param t0 the onset of the transition
#' @param dt The transition duration
#' @param y0 The initial level of the observed values (here, proxy values)
#' @param dy The increase in the observed values
#'
#' @return Returns a vector y-axis corresponding to \code{t}
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{linrampfitter}}
#' @keywords linramp
#'
#' @examples
#' @export
linramp = function(t,t0=0,dt=1,y0=0,dy=1){
  y = numeric(length(t))
  y = y0 + dy*(t-t0)/dt
  y[t<t0]=y0
  y[t>t0+dt]=y0+dy
  return(y)
}



linramprev = function(t,t0=0,dt=1,y0=0,dy=1){
  y = numeric(length(t))
  y = y0 + dy*(t-t0)/dt

  y[t>t0]=y0
  y[t<t0+dt]=y0+dy
  return(y)
}


#' Corresponding index
#'
#' Finds the indices where events best match values in a given vector.
#'
#' @param events Vector describing the values of interest (here: either depth or age of events)
#' @param record Vector we want to search (here: either full depth or age record)
#'
#' @return Returns a vector of indices for which elements of \code{record} best matches the elements of \code{events}-
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{bremla}}
#' @keywords indices
#'
#' @examples
#' @export
which.index = function(events, record){ ## Finds which indices of 'record' that are located closest to the 'event's
  eventindexes = numeric(length(events))
  for(i in 1:length(events)){
    if(events[i] < min(record) || events[i] > max(record)){ #Gives NA if located outside range of 'record'
      warning(paste0("Event ",i,", located at ",events[i]," is outside the interval covered by 'record' (",min(record),", ",max(record),"). The event will be omitted!"))
      eventindexes[i] = NA
    }else{
      eventindexes[i] = which(abs(events[i]-record) == min(abs(events[i]-record)))
    }

  }
  #eventindexes = unique(c(1,eventindexes[!is.na(eventindexes)])) #Placing transition at the start of record. Removing NA and duplicates
  return(eventindexes)
}

#' Simulation-summarizer
#'
#' Computes posterior marginal mean and uncertainty intervals from simulations.
#'
#' @param object List object which is the output of function \code{\link{bremla_chronology_simulation}}
#' @param interval character describing which uncertainty intervals should be used. Highest posterior density is used by default
#' @param print.progress Boolean describing whether or not progress should be printed on screen.
#'
#' @return Returns the \code{object} list from the input and appends additional summary statistics: posterior marginal mean and uncertainty interval-
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{bremla_chronology_simulation}}
#' @keywords indices
#'
#' @examples
#' @export
#' @importFrom matrixStats rowMeans2 rowVars
#' @importFrom INLA inla.hpdmarginal
bremla_simulationsummarizer = function(object,interval="hpd",print.progress=FALSE){
  if(print.progress) cat("Computing posterior marginal mean and 95% hpd intervals from chronology samples...\n",sep="")
  time.start = Sys.time()
  n = dim(object$simulation$age)[1]
  nsims = dim(object$simulation$age)[2]
  hpdlower = numeric(n); hpdupper = numeric(n)
  meanvek = rowMeans2(object$simulation$age)
  sdvek = sqrt(rowVars(object$simulation$age))
  for(i in 1:n){
    dens = density(object$simulation$age[i,])
    hpdlower[i] = inla.hpdmarginal(0.95,dens)[1]
    hpdupper[i] = inla.hpdmarginal(0.95,dens)[2]
  }
  time.summary = Sys.time()
  object$simulation$summary = list(mean=meanvek,sd=sdvek,hpd0.025=hpdlower,hpd0.975=hpdupper, sim.sum.time=time.summary,.args=list(interval=interval,print.progress=print.progress))
  if(print.progress) cat(" completed in ",difftime(time.summary,time.start,units="secs")[[1]],"\n",sep="")
  object$time$samplesummary = list(total=difftime(time.summary,time.start,units="secs")[[1]])
  return(object)
}
