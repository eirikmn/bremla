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
  eventindexes = unique(c(1,eventindexes[!is.na(eventindexes)])) #Placing transition at the start of record. Removing NA and duplicates
  return(eventindexes)
}

bremla_simulationsummarizer = function(object,print.progress=FALSE){
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
  if(print.progress) cat(" completed in ",difftime(time.summary,time.start,units="secs")[[1]],"\n",sep="")
  object$simulation$summary = list(mean=meanvek,sd=sdvek,hpd0.025=hpdlower,hpd0.975=hpdupper, sim.sum.time=time.summary)
  object$time$samplesummary = list(total=difftime(time.summary,time.start,units="secs")[[1]])
  return(object)
}
