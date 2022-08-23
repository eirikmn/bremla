#' Bremla model fitting
#'
#' Fits the bremla linear regression model to observations.
#'
#' @param object List object which is the output of function \code{\link{bremla_prepare}}
#' @param control.fit List containing specifications for fitting procedure. See
#' \code{\link{control.fit.default}} for details.
#' @param print.progress Boolean. If \code{TRUE} progress will be printed to screen
#'
#' @return Returns the same \code{object} list from the input, but appends information
#' from the fitting procedure in \code{object\$fitting}. Including fitted values, residuals, posterior distributions
#' and summary statistics for both fixed and random effects.
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{bremla_prepare},\link{bremla_chronology_simulation}}
#' @keywords bremla fitting
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
#' object = bremla_prepare(formula,data,nsims=5000,reference.label="simulated timescale",
#'                         events = events,
#'                         control.fit=control.fit)
#' object = bremla_modelfitter(object, print.progress=TRUE)
#' summary(object)
#' plot(object)
#' }
#' @export
#' @importFrom INLA inla inla.tmarginal inla.zmarginal inla.ar.pacf2phi
#' @importFrom stats acf arima density lm sd

bremla_modelfitter = function(object, control.fit,
                              print.progress=FALSE){

  if( .Platform$OS.type=="windows" && nrow(object$data)>7000){
    warning("Windows is poorly suited for large data sets. Try a different operating system if one is available.
              Using Windows Subsystem for Linux (WSL) might also provide a solution.")
    return(object)
  }

  if(missing(control.fit)){
    if(!is.null(object$.args$control.fit)){
      if(print.progress){
        cat("'control.fit' missing. Importing information from 'object'.",sep="")
      }
      control.fit = object$.args$control.fit
    }else{
        stop("Could not find 'control.fit'. Stopping.")
    }
  }



  #if(!is.null(control.fit))
  control.fit = set.options(control.fit,control.fit.default())
  object$.args$control.fit = control.fit




  method = control.fit$method
  noise = control.fit$noise

  time.start = Sys.time()

  formula = object$formula
  if(print.progress) cat("Performing least squares fit...",sep="")

  fit = lm(object$.internal$formula.ls,object$data) #fits model using least squares


  if(sum(is.na(fit$coefficients))){
    stop(paste0("least squares yields NA estimates for ",sum(is.na(fit$coefficients)),
                " coefficients. Make sure 'formula' does not include linearly dependent variables.
                If control.fit$degree==2, try omitting intercept (add -1 to 'formula'),
                depth and depth from 'formula'."))
  }

  dmean = fit$fitted.value
  resi = fit$residuals

  ## extract parameter estimates
  object$fitting = list(LS = list(fit=fit, params=list(meanvector=as.numeric(dmean))))

  if(tolower(object$.args$control.fit$noise) %in% c("iid","independent","ar(0)")){
    sigma = sd(resi)
    object$fitting$LS$params[["sigma"]] = sigma
  }else if(tolower(object$.args$control.fit$noise) %in% c("ar1",1,"ar(1)")){
    noisefit = arima(resi,order = c(1,0,0))
    phi = noisefit$coef[1]
    sigma = sd(resi)
    object$fitting$LS$params[["sigma"]] = sigma
    object$fitting$LS$params[["phi"]] = phi
    object$fitting$LS$noisefit=noisefit
  }else if(tolower(object$.args$control.fit$noise) %in% c("ar2",2,"ar(2)")){
    noisefit = arima(resi,order = c(2,0,0))
    phi1 = noisefit$coef[1]
    phi2 = noisefit$coef[2]
    sigma = sd(resi)
    object$fitting$LS$params[["sigma"]] = sigma
    object$fitting$LS$params[["phi1"]] = phi1
    object$fitting$LS$params[["phi2"]] = phi2
    object$fitting$LS$noisefit=noisefit
  }else{
    stop(paste0("Noise model is set to ",noise,". Only iid, ar1 and ar2 are currently supported!"))
  }
  time.ls=Sys.time()
  if(print.progress) cat(" completed!\n",sep="")
  #}



  if(tolower(control.fit$method) == "inla"){

    ## will use results from least squares fit as starting point in INLA optimization. Requires proper parametrization
    if(print.progress) cat("Performing INLA fit...",sep="")

    #set initial values for fixed parameters based on least squares 'fit'
    my.control.fixed = control.fixed.priors(object$.internal$lat.selection, fit,
                                            object$.args$events$nevents)
    object$.internal$initial_fixed = my.control.fixed
    resi = object$fitting$LS$fit$residuals

    initialmodes = log(1/object$fitting$LS$params$sigma^2)

    if(tolower(object$.args$control.fit$noise) %in% c(1,"ar1","ar(1)")){
      phi.ls = object$fitting$LS$params$phi
      initialmodes = c(initialmodes, log( (1+phi)/(1-phi) ))

    }else if(tolower(object$.args$control.fit$noise) %in% c(2,"ar2","ar(2)")){
      phi1.ls = object$fitting$LS$params$phi1
      phi2.ls = object$fitting$LS$params$phi2

      initialmodes=c(initialmodes, log( (1+phi1/(1-phi2))/(1-phi1/(1-phi2)) ), log( (1+phi2)/(1-phi2) ) )
    }
    object$.internal$initial_hyper=initialmodes

    num.threads=object$.args$control.fit$ncores #rgeneric can sometimes be more stable ifonly one core is used

    object$data$idy=1:nrow(object$data) #create covariate for random effect in INLA

    my.control.fixed$prec=1

    ## fit using INLA
    inlafit = inla(object$formula, family="gaussian",data=object$data,
                   control.family=list(hyper=list(prec=list(initial=12, fixed=TRUE))) ,
                   #control.fixed=my.control.fixed,
                   num.threads = object$.args$control.fit$ncores,
                   control.compute=list(config=TRUE),
                   #verbose=TRUE,
                   verbose=object$.args$control.fit$verbose,
                   control.inla=list(restart=TRUE,h=0.005,dz=0.75,
                                     int.strategy="auto"),
                   control.mode=list(theta=initialmodes,restart=TRUE)
                   )

    if(print.progress){
      cat(" completed.\n",sep="")
    }
    object$fitting$inla = list(fit=inlafit)


    if(print.progress){
      cat("Computing remaining posteriors using Monte Carlo simulation...\n",sep="")
    }

    ## extract posteriors for hyperparameters
    posterior_sigma = inla.tmarginal(function(x)1/sqrt(x),inlafit$marginals.hyperpar$`Precision for idy`); zmarg_sigma=inla.zmarginal(posterior_sigma,silent=TRUE)
    object$fitting$inla$hyperparameters = list(posteriors=list(sigma_epsilon=posterior_sigma))
    object$fitting$inla$hyperparameters$results$sigma_epsilon = zmarg_sigma
    if(tolower(object$.args$control.fit$noise)%in% c(1,"ar1","ar(1)") ){
      posterior_phi = inlafit$marginals.hyperpar$`Rho for idy`; zmarg_phi = inla.zmarginal(posterior_phi,silent=TRUE)
      object$fitting$inla$hyperparameters$posteriors$phi = posterior_phi
      object$fitting$inla$hyperparameters$results$phi = zmarg_phi

    }else if(object$.args$control.fit$noise %in% c(2,"ar2","ar(2)") ){
      hypersamples = inla.hyperpar.sample(50000,inlafit)

      p=2
      phii = hypersamples[, 2L:(2L+(p-1L))]
      phis = apply(phii, 1L, inla.ar.pacf2phi)
      posterior_phi1 = cbind(density(phis[1,])$x,density(phis[1,])$y);colnames(posterior_phi1)=c("x","y"); zmarg_phi1 = inla.zmarginal(posterior_phi1,silent=TRUE)
      posterior_phi2 = cbind(density(phis[2,])$x,density(phis[2,])$y);colnames(posterior_phi2)=c("x","y"); zmarg_phi2=inla.zmarginal(posterior_phi2,silent=TRUE)
      object$fitting$inla$hyperparameters$posteriors$phi1 = posterior_phi1
      object$fitting$inla$hyperparameters$results$phi1 = zmarg_phi1
      object$fitting$inla$hyperparameters$posteriors$phi2 = posterior_phi2
      object$fitting$inla$hyperparameters$results$phi2 = zmarg_phi2
    }


    time.inla=Sys.time()
    elapsed.inla = difftime(time.inla,time.ls,units="secs")[[1]]
    if(print.progress) cat("INLA fit completed in ",elapsed.inla, " seconds!\n",sep="")

    object$time$fit$inla = elapsed.inla
    object$time$fit$ls = difftime(time.ls,time.start,units="secs")[[1]]
    object$time$fit$total = time.inla-time.start
  }

  return(object)
}

