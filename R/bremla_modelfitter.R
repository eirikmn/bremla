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
#' summary(object)
#' }
#'
#'
#' @export
#' @importFrom INLA inla inla.tmarginal inla.zmarginal inla.ar.pacf2phi
#' @importFrom stats acf arima density lm sd

bremla_modelfitter = function(object, control.fit,
                              print.progress=FALSE){

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
    if(print.progress) cat("Performing INLA fit...\n",sep="")

    #set initial values for fixed parameters based on least squares 'fit'
    my.control.fixed = control.fixed.priors(object$.internal$lat.selection, fit,
                                            object$.args$events$nevents)
    object$.internal$initial_fixed = my.control.fixed
    resi = object$fitting$LS$fit$residuals

    initialmodes = log(1/object$fitting$LS$params$sigma^2)

    if(tolower(object$.args$control.fit$noise) %in% c(1,"ar1","ar(1)")){
      phi.ls = object$fitting$LS$params$phi
      initialmodes = c(initialmodes, log( (1+phi)/(1-phi) ))

    }else if(tolower(object$.args$noise) %in% c(2,"ar2","ar(2)")){
      phi1.ls = object$fitting$LS$params$phi1
      phi2.ls = object$fitting$LS$params$phi2

      initialmodes=c(initialmodes, log( (1+phi1/(1-phi2))/(1-phi1/(1-phi2)) ), log( (1+phi2)/(1-phi2) ) )
    }
    object$.internal$initial_hyper=initialmodes

    num.threads=object$.args$control.fit$ncores #rgeneric can sometimes be more stable ifonly one core is used

    object$data$idy=1:nrow(object$data) #create covariate for random effect in INLA

    ## fit using INLA
    inlafit = inla(object$formula, family="gaussian",data=object$data,
                   control.family=list(hyper=list(prec=list(initial=12, fixed=TRUE))) ,
                   #control.fixed=my.control.fixed,
                   num.threads = object$.args$control.fit$ncores,
                   control.compute=list(config=TRUE),
                   verbose=object$.args$control.fit$verbose,
                   control.inla=list(restart=TRUE,h=0.1),
                   #control.mode=list(theta=initialmodes,restart=TRUE)
                   )


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

