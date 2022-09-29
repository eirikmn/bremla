<<<<<<< HEAD
#' Mean vector from fixed effects
#'
#' Gives the mean number of layer differences per depth from the fixed effects.
#'
#' @param coefs Fixed effects
#' @param reg.model list specifying the structure of the linear regression model.
#' @param data data.frame obtained from \code{\link{bremla_prepare}}.
#'
#' @return Returns a vector representing the mean number of layer difference per depth interval.
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{bremla_prepare},\link{bremla_chronology_simulation}}
#' @keywords bremla mean
#'
#' @export
meanmaker = function(coefs,reg.model,data){
  ## Computes mean vector from given fixed effects 'coefs'.
  ## Requires specification of which effects to include ('reg.model'), the number of climate transitions ('nevents') and a data.frame with covariates ('data')

  n=nrow(data)
  fitted = numeric(n)
  for (i in 1:length(coefs)){
    name = names(reg.model)[i]
    if(name == "`(Intercept)`"){
      fitted = fitted + coefs[i]
    }else{
      fitted = fitted + coefs[i]*data[[name]]
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
#' @param suppress Boolean. Set to \code{TRUE} to suppress warnings that events fall outside designated
#' interval.
#'
#' @return Returns a vector of indices for which elements of \code{record} best matches the elements of \code{events}-
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{bremla}}
#' @keywords indices
#'
#' @export
which.index = function(events, record,suppress=TRUE){ ## Finds which indices of 'record' that are located closest to the 'event's
  eventindexes = numeric(length(events))
  for(i in 1:length(events)){
    if(events[i] < min(record) || events[i] > max(record)){ #Gives NA if located outside range of 'record'
      if(!suppress){
        warning(paste0("Event ",i,", located at ",events[i]," is outside the interval covered by 'record' (",min(record),", ",max(record),"). The event will be omitted!\n"))
      }

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
#' control.sim=list(synchronized=2,
#'                  summary=list(compute=TRUE))
#'
#' object = bremla_prepare(formula,data,nsims=5000,reference.label="simulated timescale",
#'                         events = events,
#'                         control.fit=control.fit,
#'                         control.sim=control.sim)
#' object = bremla_modelfitter(object)
#' object = bremla_chronology_simulation(object)
#' object = bremla_simulationsummarizer(object,sync=FALSE,print.progress=TRUE)
#' summary(object)
#' plot(object)
#'
#' }
#' @export
#' @importFrom matrixStats rowMeans2 rowSds rowMedians
#' @importFrom INLA inla.hpdmarginal
#' @importFrom stats median
bremla_simulationsummarizer = function(object,sync=TRUE,print.progress=FALSE){

  CI.type=object$.args$control.sim$summary$CI.type


  #sync=object$.args$control.sim$synchronized
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
        dens = density(samples[i,]); dens = data.frame(x=dens$x,y=dens$y)

        hpds = suppressWarnings( inla.hpdmarginal(0.95,inla.smarginal(dens)) )

        lower[i] = hpds[1]
        upper[i] = hpds[2]
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

#' @importFrom stringr str_replace_all
cleanstring <- function(string){
  str=str_replace_all(string," ","")
  return(str)
}

#' @importFrom stringr str_locate str_detect str_locate_all str_sub
lat.selector <- function(formulastring){
  formulastring=cleanstring(formulastring)
  lat.selection=list()
  tildepos = str_locate(formulastring,"~")[1]
  if(str_detect(formulastring,"-1+")){
    explicit_int= TRUE
  }else{
    lat.selection$`(Intercept)`=1
    if(str_detect(formulastring,"~1+")){
      explicit_int= TRUE
    }else{
      explicit_int= FALSE
    }
  }
  pluspos = (str_locate_all(formulastring,"\\+")[[1]])[,1]


  if(explicit_int){ ##omit first position (intercept) if stated explicitly in formula
    separators = pluspos
  }else{
    separators = c(tildepos,pluspos)
  }
  if(length(separators)>1){
    for(i in 1:(length(separators)-1) ){
      cov.name = str_sub(formulastring,start=separators[i]+1,end=separators[i+1]-1)
      lat.selection[[cov.name]]=1
    }
  }
  if(length(separators)>=1){
    cov.name = str_sub(formulastring,start=tail(separators,1)+1,end=-1)
    lat.selection[[cov.name]]=1
  }

  return(lat.selection)
}




control.fixed.priors = function(reg.model, fit, nevents){

  #see '?control.fixed' to see what I am doing
  my.control.fixed = list(mean=list(  ))

  startindex=1
  if(!is.null(reg.model$`(Intercept)`) && reg.model$`(Intercept)` == 1){
    my.control.fixed$mean.intercept = fit$coefficients[["(Intercept)"]]
    startindex=startindex+1
  }
  for (i in startindex:length(reg.model)){
    varname = names(reg.model)[i]
    my.control.fixed$mean[[ varname ]] = fit$coefficients[[varname]]
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
#' Connecting the Greenland ice-core and U/Th timescales via cosmogenic radionuclides: testing the synchroneity of Dansgaard-Oeschger events,
#' Clim. Past, 14, 1755-1781, https://doi.org/10.5194/cp-14-1755-2018
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
#'
#' Wrapper function which calls specified functions to produce samples from tie-point distributions and places them in \code{bremla} object.
#'
#' @param object Bremla list object to store samples in
#' @param synchronization List containing specifications for generating tie-points.
#' Must include \code{synchronization\$locations}. Includes different methods for
#' producing tie-points specified by \code{synchronization\$method}, and can also
#' import pre-computed samples in \code{synchronization\$samples}. See
#' \code{\link{synchronization.default}} for more details.
#' @param print.progress Boolean. If \code{TRUE} progress will be printed to screen.
#' @param ... Further arguments to be passed down to the selected tie-point simulation
#' function.
#'
#' @return returns the same \code{object} from input but appends tie-point samples
#' and information in \code{object\$tie_points}.
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
#' synchronization=list(locations=depth[c(100,400,700)],method="gauss",locations_unit="depth",
#'         params=list(mean=age[c(100,400,700)]+c(20,5,-20),sd=c(20,10,50)))
#' object = bremla_prepare(formula,data,nsims=5000,
#'                         reference.label="simulated timescale",
#'                         events = events,
#'                         synchronization=synchronization)
#' object = tiepointsimmer(object)
#' summary(object)
#' plot(object)
#' }
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @keywords tiepoint sample
#' @export
tiepointsimmer = function(object, synchronization,print.progress=FALSE,...){

  time.start = Sys.time()
  if(missing(synchronization)){

    if(!is.null(object$.args$synchronization)){
      if(print.progress){
        cat("'synchronization' missing. Importing information from 'object'.",sep="")
      }
      synchronization = object$.args$synchronization
    }else{
      stop("Could not find 'synchronization'. Stopping.")
    }
  }
  #if(!is.null(synchronization))
  synchronization = set.options(synchronization,synchronization.default())

  object$.args$synchronization = synchronization


  if(!is.null(synchronization$samples)){
    nsims = ncol(synchronization$samples)
    if(ncol(synchronization$samples) != length(synchronization$locations))
      warning("number of columns in 'samples' must correspond to length of 'locations'!")
    samples = synchronization$samples


  }else{
    if(synchronization$method=="adolphi"){
      synchronization$locations = c(11050,12050,13050,22050,42050)
      synchronization$locations_unit="age"
      locations_indexes = which.index(synchronization$locations,object$data$age)
      tieshifts = synchronization$locations
      synchronization$x.ref=tieshifts

      samples = adolphi_tiepoint_simmer(nsims=synchronization$nsims,
                                        tieshifts=tieshifts,...)
    }else if(synchronization$method %in% c("normal","gaussian","gauss")){
      #synchronization$locations = object$data$depth[c(100,400,700)]

      if(tolower(synchronization$locations_unit) %in% c("age","y","time")){
        locations_indexes = which.index(synchronization$locations,object$data$age)
      }else if(tolower(synchronization$locations_unit) %in% c("depth","z")){
        locations_indexes = which.index(synchronization$locations,object$data$depth)
      }else{
        locations_indexes = synchronization$locations
      }

      ntie = length(synchronization$locations)
      if(is.null(synchronization$params)){
        no_offset = object$original.chron$age[locations_indexes]
        synchronization$params$mean = rep(no_offset,ntie)
        synchronization$params$sd = rep(1,ntie)
      }
      means = synchronization$params$mean
      sds = synchronization$params$sd
      samples = matrix(NA,nrow=synchronization$nsims,ncol=ntie)

      #loc.ind = which.index(synchronization$locations,object$data$depth)

      for(i in 1:ntie){
        samples[,i] = rnorm(synchronization$nsims,mean=means[i],sd=sds[i])
      }
      ## temporary

      # meanvek = object$data$age[loc.ind]+c(1,-100,300)
      # sdvek = object$data$age[loc.ind]/meanvek[1]*(1:3)*50
      # for(i in 1:ntie){
      #   samples[,i] = rnorm(synchronization$nsims,mean=meanvek[i],sd=sdvek[i])
      # }
      #locations_indexes=loc.ind
    }



  }
  if(synchronization$locations_unit %in% c("depth","z")){
    locations_indexes = which.index(synchronization$locations,object$data$depth)
  }else if(synchronization$locations_unit %in% c("age","y")){
    locations_indexes = which.index(synchronization$locations,object$data$age)
  }else if(synchronization$locations_unit %in% c("indexes","index","i")){
    locations_indexes=synchronization$locations
  }

  object$tie_points = list(samples=samples,
                           locations=synchronization$locations,
                           locations_unit=synchronization$locations_unit,
                           method=synchronization$method,
                           nsims=synchronization$nsims,
                           x.ref=synchronization$x.ref,
                           locations_indexes=locations_indexes,
                           tie_n=length(synchronization$locations))

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
#' @examples
#' n=100
#' sigma=1
#' rho=0.8
#' Q = Qmaker_ar1cum(n,sigma,rho)
#'
#' nsims = 1000
#' muvek = seq(1,5,length.out=n)
#' samples = Qsimmer(nsims,Q,muvek)
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @keywords simulation GMRF sparse
#' @export
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
#' @examples
#' n=100
#' sigma=1
#' rho=0.8
#' Q = Qmaker_ar1cum(n,sigma,rho)
#'
#' nsims = 1000
#' muvek = seq(1,5,length.out=n)
#' samples = Qsimmer(nsims,Q,muvek)
#'
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @keywords GMRF sparse precision matrix
#' @seealso \code{\link{Qsimmer}, \link{Qymaker}}
#' @export
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
#' @examples
#' n=100
#' sigma=1
#' rho=0.8
#' Q_method1 = Qmaker_ar1cum(n,sigma,rho)
#' Q_ar1 = Matrix::sparseMatrix(
#'        i = c(1,n,2:(n-1),1:(n-1)),
#'        j = c(1,n,2:(n-1),2:n),
#'        x = 1/(1-rho^2)*1/sigma^2*c(1,1,rep(1+rho^2,n-2),rep(-rho,n-1)),
#'        symmetric = TRUE
#' )
#' Q_method2 = Qymaker(Q_ar1)
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @keywords GMRF sparse precision matrix
#' @seealso \code{\link{Qsimmer}, \link{Qmaker_ar1cum}}
#' @export
Qymaker = function(Qx){
  nn = dim(Qx)[1]
  ii = c(1:nn,2:nn); jj = c(1:nn,1:(nn-1)); xx = c(rep(1,nn),rep(-1,nn-1))
  S = sparseMatrix(i=ii,j=jj,x=xx)
  return(t(S)%*%Qx%*%S)
}
=======
#' Mean vector from fixed effects
#'
#' Gives the mean number of layer differences per depth from the fixed effects.
#'
#' @param coefs Fixed effects
#' @param reg.model list specifying the structure of the linear regression model.
#' @param data data.frame obtained from \code{\link{bremla_prepare}}.
#'
#' @return Returns a vector representing the mean number of layer difference per depth interval.
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{bremla_prepare},\link{bremla_chronology_simulation}}
#' @keywords bremla mean
#'
#' @export
meanmaker = function(coefs,reg.model,data){
  ## Computes mean vector from given fixed effects 'coefs'.
  ## Requires specification of which effects to include ('reg.model'), the number of climate transitions ('nevents') and a data.frame with covariates ('data')

  n=nrow(data)
  fitted = numeric(n)
  for (i in 1:length(coefs)){
    name = names(reg.model)[i]
    if(name == "`(Intercept)`"){
      fitted = fitted + coefs[i]
    }else{
      fitted = fitted + coefs[i]*data[[name]]
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
#' control.sim=list(synchronized=2,
#'                  summary=list(compute=TRUE))
#'
#' object = bremla_prepare(formula,data,nsims=5000,reference.label="simulated timescale",
#'                         events = events,
#'                         control.fit=control.fit,
#'                         control.sim=control.sim)
#' object = bremla_modelfitter(object)
#' object = bremla_chronology_simulation(object)
#' object = bremla_simulationsummarizer(object,sync=FALSE,print.progress=TRUE)
#' summary(object)
#' plot(object)
#'
#' }
#' @export
#' @importFrom matrixStats rowMeans2 rowSds rowMedians
#' @importFrom INLA inla.hpdmarginal
#' @importFrom stats median
bremla_simulationsummarizer = function(object,sync=TRUE,print.progress=FALSE){

  CI.type=object$.args$control.sim$summary$CI.type


  #sync=object$.args$control.sim$synchronized
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

#' @importFrom stringr str_replace_all
cleanstring <- function(string){
  str=str_replace_all(string," ","")
  return(str)
}

#' @importFrom stringr str_locate str_detect str_locate_all str_sub
lat.selector <- function(formulastring){
  formulastring=cleanstring(formulastring)
  lat.selection=list()
  tildepos = str_locate(formulastring,"~")[1]
  if(str_detect(formulastring,"-1+")){
    explicit_int= TRUE
  }else{
    lat.selection$`(Intercept)`=1
    if(str_detect(formulastring,"~1+")){
      explicit_int= TRUE
    }else{
      explicit_int= FALSE
    }
  }
  pluspos = (str_locate_all(formulastring,"\\+")[[1]])[,1]


  if(explicit_int){ ##omit first position (intercept) if stated explicitly in formula
    separators = pluspos
  }else{
    separators = c(tildepos,pluspos)
  }
  if(length(separators)>1){
    for(i in 1:(length(separators)-1) ){
      cov.name = str_sub(formulastring,start=separators[i]+1,end=separators[i+1]-1)
      lat.selection[[cov.name]]=1
    }
  }
  if(length(separators)>=1){
    cov.name = str_sub(formulastring,start=tail(separators,1)+1,end=-1)
    lat.selection[[cov.name]]=1
  }

  return(lat.selection)
}




control.fixed.priors = function(reg.model, fit, nevents){

  #see '?control.fixed' to see what I am doing
  my.control.fixed = list(mean=list(  ))

  startindex=1
  if(!is.null(reg.model$`(Intercept)`) && reg.model$`(Intercept)` == 1){
    my.control.fixed$mean.intercept = fit$coefficients[["(Intercept)"]]
    startindex=startindex+1
  }
  for (i in startindex:length(reg.model)){
    varname = names(reg.model)[i]
    my.control.fixed$mean[[ varname ]] = fit$coefficients[[varname]]
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
#' Connecting the Greenland ice-core and U/Th timescales via cosmogenic radionuclides: testing the synchroneity of Dansgaard-Oeschger events,
#' Clim. Past, 14, 1755-1781, https://doi.org/10.5194/cp-14-1755-2018
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
#'
#' Wrapper function which calls specified functions to produce samples from tie-point distributions and places them in \code{bremla} object.
#'
#' @param object Bremla list object to store samples in
#' @param synchronization List containing specifications for generating tie-points.
#' Must include \code{synchronization\$locations}. Includes different methods for
#' producing tie-points specified by \code{synchronization\$method}, and can also
#' import pre-computed samples in \code{synchronization\$samples}. See
#' \code{\link{synchronization.default}} for more details.
#' @param print.progress Boolean. If \code{TRUE} progress will be printed to screen.
#' @param ... Further arguments to be passed down to the selected tie-point simulation
#' function.
#'
#' @return returns the same \code{object} from input but appends tie-point samples
#' and information in \code{object\$tie_points}.
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
#' synchronization=list(locations=depth[c(100,400,700)],method="gauss")
#' object = bremla_prepare(formula,data,nsims=5000,
#'                         reference.label="simulated timescale",
#'                         events = events,
#'                         synchronization=synchronization)
#' object = tiepointsimmer(object)
#' summary(object)
#' plot(object)
#' }
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @keywords tiepoint sample
#' @export
tiepointsimmer = function(object, synchronization,print.progress=FALSE,...){

  time.start = Sys.time()
  if(missing(synchronization)){

    if(!is.null(object$.args$synchronization)){
      if(print.progress){
        cat("'synchronization' missing. Importing information from 'object'.",sep="")
      }
      synchronization = object$.args$synchronization
    }else{
      stop("Could not find 'synchronization'. Stopping.")
    }
  }
  #if(!is.null(synchronization))
  synchronization = set.options(synchronization,synchronization.default())

  object$.args$synchronization = synchronization


  if(!is.null(synchronization$samples)){
    nsims = ncol(synchronization$samples)
    if(ncol(synchronization$samples) != length(synchronization$locations))
      warning("number of columns in 'samples' must correspond to length of 'locations'!")
    samples = synchronization$samples
  }else{
    if(synchronization$method=="adolphi"){
      synchronization$locations = c(11050,12050,13050,22050,42050)
      synchronization$locations_unit="age"
      locations_indexes = which.index(synchronization$locations,object$data$age)
      tieshifts = synchronization$locations
      synchronization$x.ref=tieshifts

      samples = adolphi_tiepoint_simmer(nsims=synchronization$nsims,
                                        tieshifts=tieshifts,...)
    }else if(synchronization$method %in% c("normal","gaussian","gauss")){
      synchronization$locations = object$data$depth[c(100,400,700)]

      ntie = length(synchronization$locations)
      ## temporary
      samples = matrix(NA,nrow=synchronization$nsims,ncol=ntie)
      loc.ind = which.index(synchronization$locations,object$data$depth)
      meanvek = object$data$age[loc.ind]+c(1,-100,300)
      sdvek = object$data$age[loc.ind]/meanvek[1]*(1:3)*50
      for(i in 1:ntie){
        samples[,i] = rnorm(synchronization$nsims,mean=meanvek[i],sd=sdvek[i])
      }
      locations_indexes=loc.ind
    }



  }

  object$tie_points = list(samples=samples,
                           locations=synchronization$locations,
                           locations_unit=synchronization$locations_unit,
                           method=synchronization$method,
                           nsims=synchronization$nsims,
                           x.ref=synchronization$x.ref,
                           locations_indexes=locations_indexes,
                           tie_n=length(synchronization$locations))

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
#' @examples
#' n=100
#' sigma=1
#' rho=0.8
#' Q = Qmaker_ar1cum(n,sigma,rho)
#'
#' nsims = 1000
#' muvek = seq(1,5,length.out=n)
#' samples = Qsimmer(nsims,Q,muvek)
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @keywords simulation GMRF sparse
#' @export
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
#' @examples
#' n=100
#' sigma=1
#' rho=0.8
#' Q = Qmaker_ar1cum(n,sigma,rho)
#'
#' nsims = 1000
#' muvek = seq(1,5,length.out=n)
#' samples = Qsimmer(nsims,Q,muvek)
#'
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @keywords GMRF sparse precision matrix
#' @seealso \code{\link{Qsimmer}, \link{Qymaker}}
#' @export
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
#' @examples
#' n=100
#' sigma=1
#' rho=0.8
#' Q_method1 = Qmaker_ar1cum(n,sigma,rho)
#' Q_ar1 = Matrix::sparseMatrix(
#'        i = c(1,n,2:(n-1),1:(n-1)),
#'        j = c(1,n,2:(n-1),2:n),
#'        x = 1/(1-rho^2)*1/sigma^2*c(1,1,rep(1+rho^2,n-2),rep(-rho,n-1)),
#'        symmetric = TRUE
#' )
#' Q_method2 = Qymaker(Q_ar1)
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @keywords GMRF sparse precision matrix
#' @seealso \code{\link{Qsimmer}, \link{Qmaker_ar1cum}}
#' @export
Qymaker = function(Qx){
  nn = dim(Qx)[1]
  ii = c(1:nn,2:nn); jj = c(1:nn,1:(nn-1)); xx = c(rep(1,nn),rep(-1,nn-1))
  S = sparseMatrix(i=ii,j=jj,x=xx)
  return(t(S)%*%Qx%*%S)
}
>>>>>>> parent of 19db1dc (Added some more features and Rmd files to generate the plots and results from both papers)
