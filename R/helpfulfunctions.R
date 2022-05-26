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
#' @export
linramp = function(t,t0=0,dt=1,y0=0,dy=1){
  y = numeric(length(t))
  y = y0 + dy*(t-t0)/dt
  y[t<t0]=y0
  y[t>t0+dt]=y0+dy
  return(y)
}


#' Linear ramp function (reverse)
#'
#' Computes the (noiseless) linear ramp function, but in reverse.
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
#' @export
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
#' @param CI.type character describing which uncertainty intervals should be used (\code{"quantiles" or "hpd"}). Highest posterior density is used by default
#' @param sync boolean. Set \code{TRUE} if summary statistics should be computed for syncrhonized samples (if any). \code{FALSE} for unsynchronized.
#' @param print.progress Boolean describing whether or not progress should be printed on screen.
#'
#' @return Returns the \code{object} list from the input and appends additional summary statistics: posterior marginal mean and uncertainty interval-
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{bremla_chronology_simulation}}
#' @keywords indices
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
#' object = bremla_prepare(age,depth,proxy,events=eventdepths,nsims=0)
#' object = bremla_modelfitter(object)
#' object = bremla_chronology_simulation(object,nsims=1000)
#' object = bremla_simulationsummarizer(object,CI.type="quant",sync=FALSE)
#' plot(object)
#' }
#' @export
#' @importFrom matrixStats rowMeans2 rowSds rowMedians
#' @importFrom INLA inla.hpdmarginal
#' @importFrom stats median
bremla_simulationsummarizer = function(object,CI.type="hpd",sync=TRUE,print.progress=FALSE){
  if(print.progress) cat("Computing posterior marginal mean and 95% credible intervals from chronology samples...\n",sep="")
  time.start = Sys.time()
  if(sync && !is.null(object$simulation$age_sync)){
    if(is.null(object$simulation$age_sync)){
      warning("Could not find synchronous simulations. Try again with sync=FALSE")
      return(object)
    }else{
      samples = object$simulation$age_sync
    }
  }else{
    samples = object$simulation$age
  }
  n = nrow(samples)
  nsims = ncol(samples)

  meanvek = rowMeans2(samples)
  sdvek = rowSds(samples)
  medianvek = rowMedians(samples)
  lower = numeric(n); upper = numeric(n)
  if(CI.type=="hpd"){

    for(i in 1:n){
      if(sync && !is.null(object$simulation$age_sync) && length(unique(samples[i,]))==1){ #if the tie-points are fixed

        lower[i] = samples[i,1]
        upper[i] = samples[i,1]
      }else{
        dens = density(samples[i,])


        lower[i] = inla.hpdmarginal(0.95,dens)[1]
        upper[i] = inla.hpdmarginal(0.95,dens)[2]
      }
    }
  }else{
    lower = meanvek-1.96*sdvek
    upper = meanvek+1.96*sdvek
  }

  time.summary = Sys.time()
  if(sync){
    object$simulation$summary_sync = list(mean=meanvek,median=medianvek,sd=sdvek,lower=lower,upper=upper,
                                     .args=list(print.progress=print.progress,CI.type=CI.type))
  }else{
    object$simulation$summary = list(mean=meanvek,median=medianvek,sd=sdvek,lower=lower,upper=upper,
                                     .args=list(print.progress=print.progress,CI.type=CI.type))
  }

  #if(CI.type=="hpd") object$simulation$summary$mode = modevek
  if(print.progress) cat(" completed in ",difftime(time.summary,time.start,units="secs")[[1]],"seconds.\n",sep="")

  object$time$samplesummary = list(total=difftime(time.summary,time.start,units="secs")[[1]])

  if(sync){
    object$time$samplesyncsummary = list(total=difftime(time.summary,time.start,units="secs")[[1]])
    object$simulation$summary_sync$sim.sum.time = difftime(time.summary,time.start,units="secs")[[1]]
  }else{
    object$time$samplesummary = list(total=difftime(time.summary,time.start,units="secs")[[1]])
    object$simulation$summary$sim.sum.time = difftime(time.summary,time.start,units="secs")[[1]]
  }


  return(object)
}



control.fixed.priors = function(reg.model, fit, nevents){

  my.control.fixed = list(mean=list(  ))

  if(reg.model$depth1) my.control.fixed$mean[["z1"]] = fit$coefficients[["z1"]]
  if(reg.model$depth2) my.control.fixed$mean[["z2"]] = fit$coefficients[["z2"]]
  if(reg.model$proxy) my.control.fixed$mean[["x"]] = fit$coefficients[["x"]]

  if(reg.model$psi0 || reg.model$psi1){
    for(i in 1:(nevents-1)){
      if(reg.model$psi0){
        my.control.fixed$mean[[paste0("a",i)]] = fit$coefficients[[paste0("a",i)]]
      }
      if(reg.model$psi0){
        my.control.fixed$mean[[paste0("c",i)]] = fit$coefficients[[paste0("c",i)]]
      }
    }
  }

  return(my.control.fixed)
}

#' Skew-quasi-Gaussian sampler
#' Produces simulations from the merged Gaussian process
#' @param n The number of samples to produce.
#' @param mode Mode where the two Gaussians are merged.
#' @param sdL lower standard deviation.
#' @param sdU upper standard deviation.
#' @param log boolean describing whether or not the log value should be returned.
#' @param plothist list object with arguments describing how (and if) a histogram should be produced.
#'
#' @return returns a numeric object with \code{n} samples.
#'
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @keywords skew sample tiepoint
#' @export
skewsampler = function(n,mode=0,sdL=1,sdU=1,log=FALSE,plothist=list(compute=FALSE,breaks=50,col="orange")){
  uppers = runif(n) < sdU/(sdL+sdU)
  samples = numeric(n)
  for(i in 1:n){
    if(uppers[i]){
      sample = abs(rnorm(1))
      samples[i] = sdU*sample+mean
    }else{
      sample = -abs(rnorm(1))
      samples[i] = sdL*sample+mean
    }
  }
  if(log) samples = log(samples)

  return(samples)
}


#' Load Adolphi tie-point distributions
#' Returns the probability distributions for the tie-points given by Adolphi et al. (2018) and Muscheler et al. (2020)
#' @param tieshifts numeric which gives the amount each tie-point should be shifted, e.g. to go from describing offset from GICC05 to years before present.
#' @param plotdens boolean. If \code{TRUE} plot pdfs.
#' @param x.ref numeric. Gives reference value on the x-axis when distributions are plotted.
#' @return returns a list containing the pdfs for each tie-point.
#' @examples
#' adolphipdfs = adolphiloader(tieshifts=c(11050,12050,13050,22050,42050))
#' plot(adolphipdfs$tie1,type="l",xlab="Time (yb2k)",ylab="Density",main="Tie-point #1")
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @references Adolphi, F., Bronk Ramsey, C., Erhardt, T., Edwards, R. L., Cheng, H., Turney, C. S. M., Cooper, A., Svensson, A., Rasmussen, S. O., Fischer, H., and Muscheler, R. (2018).
#' Connecting the Greenland ice-core and U∕Th timescales via cosmogenic radionuclides: testing the synchroneity of Dansgaard–Oeschger events,
#' Clim. Past, 14, 1755–1781, https://doi.org/10.5194/cp-14-1755-2018
#' @references Muscheler, R., Adolphi, F., Heaton, T., Bronk Ramsey, C., Svensson, A., Van der Plicht, J., & Reimer, P. (2020).
#' Testing and Improving the IntCal20 Calibration Curve with Independent Records.
#' Radiocarbon, 62(4), 1079-1094. doi:10.1017/RDC.2020.54
#' @keywords adolphi tiepoint simulation
#' @export
#' @importFrom utils data
adolphiloader = function(tieshifts=numeric(5), plotdens=FALSE, x.ref=NULL){
  #data("adolphi_tiepoints",package = "bremla",envir=environment())


  #adolphi_tiepoints[adolphi_tiepoints[,3]==1,1]
  tie1 = adolphi_tiepoints[adolphi_tiepoints[,3]==1,1:2]
  tie2 = adolphi_tiepoints[adolphi_tiepoints[,3]==2,1:2]
  tie3 = adolphi_tiepoints[adolphi_tiepoints[,3]==3,1:2]
  tie4 = adolphi_tiepoints[adolphi_tiepoints[,3]==4,1:2]
  tie5 = adolphi_tiepoints[adolphi_tiepoints[,3]==5,1:2]

  returlist = list(tie1=tie1,tie2=tie2,tie3=tie3,tie4=tie4,tie5=tie5)
  for(i in 1:5){
    returlist[[i]][,1] = returlist[[i]][,1]+tieshifts[i]
  }

  if(plotdens){
    la=layout(mat=matrix(c(1,4,1,4,2,4,2,5,3,5,3,5) ,nrow=2))
    if(sum(tieshifts)==0){
      xlab = "x - GICC05 (years)"
    }else{
      xlab = "Age (y b2k)"
    }
    for(i in 1:5){
      plot(returlist[[i]],type="l",main=paste0("Tie-point ",i),xlab=xlab,ylab="Density")
      zmarg = inla.zmarginal(returlist[[i]],silent=TRUE)
      abline(v=zmarg$mean)
      abline(v=zmarg$quant0.025,col="gray")
      abline(v=zmarg$quant0.975,col="gray")
      if(!is.null(x.ref)){
        abline(v=x.ref[i],col="blue")
      }
    }
  }
  return(returlist)
}


#' Simulate from Adolphi tie-point distributions
#' Returns the probability distributions for the tie-points given by Adolphi et al. (2018) and Muscheler et al. (2020)
#' @param nsims Integer. The number of samples to be produced.
#' @param plotdens boolean. If \code{TRUE} plot pdfs.
#' @param tieshifts numeric which gives the amount each tie-point should be shifted, e.g. to go from describing offset from GICC05 to years before present.
#' @param plotdens boolean that is passed on to adolphiloader. Set \code{TRUE} if you want to plot distributions
#' @param x.ref numeric. Gives reference value on the x-axis when distributions are plotted.
#' @param plothist list object describing if and how histogram of samples should be plotted
#' @return returns a list containing the pdfs for each tie-point.
#' @examples
#' tieshifts= c(11050,12050,13050,22050,42050)
#' samples = adolphi_tiepoint_simmer(nsims=2000,tieshifts=tieshifts)
#' hist(samples[,3],col="orange",freq=0,breaks=50,main="Tie-point 3",xlab="Age (yb2k)")
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @keywords adolphi tiepoint pdf
#' @export
#' @importFrom graphics layout
adolphi_tiepoint_simmer = function(nsims=10000,tieshifts = numeric(5), plotdens=FALSE,x.ref=NULL,
                                   plothist=list(compute=FALSE,breaks=50,col="orange")){

  tie_pdfs = adolphiloader(tieshifts = tieshifts,plotdens=plotdens,x.ref=x.ref)


  tie_samples = matrix(NA,nrow=nsims,ncol=5)

  for(i in 1:5){
    sims = inla.rmarginal(nsims,tie_pdfs[[i]])
    tie_samples[,i] = sims
  }

  if(sum(tieshifts)==0){
    xlab = "x - GICC05 (years)"
  }else{
    xlab = "Age (y b2k)"
  }

  if(!is.null(plothist) && plothist$compute==TRUE){
    la=layout(mat=matrix(c(1,4,1,4,2,4,2,5,3,5,3,5) ,nrow=2))
    hist(tie_samples[,1],freq=0,col=plothist$col,breaks=plothist$breaks,main="Tie-point 1",xlab=xlab)
    hist(tie_samples[,2],freq=0,col=plothist$col,breaks=plothist$breaks,main="Tie-point 2",xlab=xlab)
    hist(tie_samples[,3],freq=0,col=plothist$col,breaks=plothist$breaks,main="Tie-point 3",xlab=xlab)
    hist(tie_samples[,4],freq=0,col=plothist$col,breaks=plothist$breaks,main="Tie-point 4",xlab=xlab)
    hist(tie_samples[,5],freq=0,col=plothist$col,breaks=plothist$breaks,main="Tie-point 5",xlab=xlab)
    par(mfrow=c(1,1))
  }


  #par(mfrow=c(1,1))
  return(tie_samples)
}


#' Tie-point simulation
#' Wrapper function which calls specified functions to produce samples from tie-point distributions and places them in \code{bremla} object.
#' @param object Bremla list object to store samples in
#' @param nsims Integer. The number of samples to be produced.
#' @param locations numeric describing location of each tie-point
#' @param locations.type string. what unit is the location given in? Can be "depth", "age" and "index".
#' @param method string. Set \code{method="adolphi"} to use adolphi tie-points. Other options will be made available soon
#' @param x.ref numeric of length \code{length(locations)}, adds reference locations (x-axis) to which the tie-points can be compared to later.
#' @param samples data.frame including premade tie-point samples. \code{ncol} must correspond to the number of tie-points.
#' @param ... Further arguments to be passed down to the selected tie-point simulation function.
#'
#' @return returns the same \code{object} from input but appends tie-point samples and information.
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
#' object = bremla_prepare(age,depth,proxy,events=eventdepths,nsims=0)
#' object = tiepointsimmer(object,nsims=10000,method="adolphi")
#' hist(object$tie_points$samples[,2],col="orange",freq=0,breaks=40,
#'     main="Tie-point 3",xlab="Age (yb2k)")
#' }
#'
#'
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @keywords tiepoint sample
#' @export
tiepointsimmer = function(object, nsims=10000, locations=NULL, locations.type="depth",
                          method="adolphi", x.ref=NULL, samples=NULL,...){

  time.start = Sys.time()

  if(!is.null(samples)){
    nsims = ncol(samples)
    if(ncol(samples) != length(locations)) warning("number of columns in 'samples' must correspond to length of 'locations'!")

    object$tie_points = list(samples=samples,locations=locations,locations.type=locations.type,
                             method=method,nsims=nsims)
  }

  if(method=="adolphi"){
    locations = c(11050,12050,13050,22050,42050)
    locations.type="age"
    tieshifts = locations
    x.ref=tieshifts

    samples = adolphi_tiepoint_simmer(nsims=nsims,tieshifts=tieshifts,...)
  }


  object$tie_points = list(samples=samples, locations=locations,locations.type=locations.type,
                           method=method,nsims=nsims,x.ref=x.ref)

  time.total = difftime(Sys.time(), time.start,units="secs")[[1]]

  object$time$tiepoints = time.total

  class(object) = "bremla"
  return(object)
}


#' Efficient Gaussian simulation
#' Simulate from multivariate Gaussian distribution efficiently by using sparse precision matrix
#' @param nsims The number of samples to be produced.
#' @param Q The precision matrix, should be sparse.
#' @param muvek The mean vector.
#'
#' @return Returns matrix with nsims samples of length equal to \code{nrow(Q)}.
#'
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @keywords simulation GMRF sparse
Qsimmer = function(nsims, Q, muvek){
  nn=dim(Q)[1]
  samples = matrix(NA,nrow=nn,ncol=nsims)

  La = chol(Q)

  for(i in 1:nsims){
    samples_z = rnorm(nn)
    v = solve(La,samples_z)[,1]
    samples[,i] = muvek[(length(muvek)-nn+1):length(muvek)] + v
  }
  return(samples)
}

#' Precision matrix of cumulative AR(1) process
#' Produces the sparse precision matrix of a cumulative AR(1) process
#' @param n The dimension of the multivariate process.
#' @param sigma The standard deviation (scaling parameter).
#' @param rho The lag-one correlation coefficient of the AR(1) process.
#'
#' @return Returns sparse precision matrix of dimension \code{n x n}.
#'
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @keywords GMRF sparse precision matrix
Qmaker_ar1cum = function(n,sigma,rho){

  noise = sigma^2*(1-rho^2)
  ii = c(1L, n, 2L:(n - 1L), 1L:(n - 1L),1L:(n-2L)); jj = c(1L, n, 2L:(n - 1L), 2L:n,3L:n)
  if(n==3){
    ii = c(1L, n, n-1L,  1L:(n - 2L), n-1L, 1L:(n-2L))
    jj = c(1L, n, n-1L,  2L:(n - 1L), n    ,3L:n)
    values = 1/(noise)*c(2+2*rho+1*rho^2, 1, rep(2+2*rho+1*rho^2,n-2L), rep(-(1+2*rho+rho^2),n-2),-(1+rho), rep(rho,n-2))
  }else{
    ii = c(1L, n, n-1L, 2L:(n-2L), 1L:(n-2L), n-1L, 1L:(n-2L))
    jj = c(1L, n, n-1L, 2L:(n-2L), 2L:(n-1),  n,    3L:n)
    values = 1/(noise)*c( 2+2*rho+rho^2, 1, rho^2+2*rho+2,rep(2+2*rho+2*rho^2,n-3L), rep(-(1+2*rho+rho^2),n-2),-(1+rho), rep(rho,n-2)  )
  }
  return(sparseMatrix(i=ii, j=jj, x=values, giveCsparse = TRUE, symmetric = TRUE))
}

#' Precision matrix of cumulative Gaussian process
#' Produces the sparse precision matrix of a cumulative Gaussian process
#' @param Qx The precision matrix of the Gaussian process to be transformed.
#'
#' @return Returns sparse precision matrix of dimension \code{n x n}.
#'
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @keywords GMRF sparse precision matrix
Qymaker = function(Qx){
  nn = dim(Qx)[1]
  ii = c(1:nn,2:nn); jj = c(1:nn,1:(nn-1)); xx = c(rep(1,nn),rep(-1,nn-1))
  S = sparseMatrix(i=ii,j=jj,x=xx)
  return(t(S)%*%Qx%*%S)
}
