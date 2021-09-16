#' Bremla model fitting
#'
#' Fits the bremla linear regression model to observations.
#'
#' @param object List object which is the output of function \code{\link{bremla_prepare}}
#' @param method Character specifying which method of inference to be used. Currently only \code{ls} and \code{inla} are supported, and if \code{inla} is chosen then both methods are performed. \code{age}
#' @param print.progress Boolean. If \code{TRUE} progress will be printed to screen
#' @param verbose Boolean. If \code{TRUE} details of the inla optimization procedure will be printed to the screen.
#'
#' @return Returns the same \code{object} list from the input, but appends information from the fitting procedure. Including fitted values, residuals and summary statistics for the least squares method. If \code{method="inla"} then the \code{inla} type object from the INLA package are included, as well as posterior marginal distributions and summary statistics for the hyperparameters.
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{bremla_prepare},\link{bremla_chronology_simulation}}
#' @keywords bremla fitting
#'
#' @examples
#'
#' @export
#' @importFrom INLA inla inla.tmarginal inla.zmarginal inla.ar.pacf2phi
#' @importFrom stats acf arima density lm sd

bremla_modelfitter = function(object, method="inla",print.progress=FALSE,verbose=FALSE){

  time.start = Sys.time()
  formula = object$formula
  if(print.progress) cat("Performing least squares fit...",sep="")
  #if(tolower(method) != "inla"){
  fit = lm(object$ls.formula,object$data)
  #plot(data$z,data$dy,typ="l",col="gray"); lines(data$z,fit$fitted.values,col="red")
  dmean = fit$fitted.value
  mean = object$data$y[1]+cumsum(dmean)
  resi = fit$residuals

  object$LS.fitting = list(fit=fit, params=list(meanvector=as.numeric(mean)))
  #object$params = list(mean=mean)
  if(tolower(object$.args$noise) %in% c("iid","independent","ar(0)")){
    sigma = sd(resi)
    object$LS.fitting$params[["sigma"]] = sigma
  }else if(tolower(object$.args$noise) %in% c("ar1",1,"ar(1)")){
    noisefit = arima(resi,order = c(1,0,0))
    phi = noisefit$coef[1]
    sigma = sd(resi)
    object$LS.fitting$params[["sigma"]] = sigma
    object$LS.fitting$params[["phi"]] = phi
    object$LS.fitting$noisefit=noisefit
  }else if(tolower(object$.args$noise) %in% c("ar2",2,"ar(2)")){
    noisefit = arima(resi,order = c(2,0,0))
    phi1 = noisefit$coef[1]
    phi2 = noisefit$coef[2]
    sigma = sd(resi)
    object$LS.fitting$params[["sigma"]] = sigma
    object$LS.fitting$params[["phi1"]] = phi1
    object$LS.fitting$params[["phi2"]] = phi2
    object$LS.fitting$noisefit=noisefit
  }
  time.ls=Sys.time()
  if(print.progress) cat(" completed!\n",sep="")
  #}

  ### move to separate function

  if(tolower(method) == "inla"){

    ## will use results from least squares fit as starting point in INLA optimization. Requires proper parametrization
    if(print.progress) cat("Performing INLA fit...\n",sep="")

    initialmodes = log(1/object$LS.fitting$params$sigma^2)
    if(tolower(object$.args$noise) %in% c(1,"ar1","ar(1)")){
      phi.ls = object$LS.fitting$params$phi
      initialmodes = c(initialmodes, log( (1+phi)/(1-phi) ))
    }else if(tolower(object$.args$noise) %in% c(2,"ar2","ar(2)")){
      phi1.ls = object$LS.fitting$params$phi1
      phi2.ls = object$LS.fitting$params$phi2
      initialmodes = c(initialmodes, log( (1+phi1/(1-phi2))/(1-phi1/(1-phi2)) ), log( (1+phi2)/(1-phi2) ) )
    }

    my.control.fixed = control.fixed.priors(reg.model, fit, nevents)

    inlafit = inla(object$formula, family="gaussian",data=object$data, control.family=list(hyper=list(prec=list(initial=12, fixed=TRUE))),
                   control.fixed=my.control.fixed,
                   control.compute=list(config=TRUE),verbose=verbose,control.inla=list(restart=TRUE,h=0.1), control.mode=list(theta=initialmodes)  )

    #control.fixed=list(mean.intercept=df$layers[1],prec.intercept=0.001)
    object$fitting = list(fit=inlafit)
    posterior_sigma = inla.tmarginal(function(x)1/sqrt(x),inlafit$marginals.hyperpar$`Precision for idy`); zmarg_sigma=inla.zmarginal(posterior_sigma,silent=TRUE)
    object$fitting$hyperparameters = list(posteriors=list(sigma_epsilon=posterior_sigma))
    object$fitting$hyperparameters$results = list(sigma_epsilon = zmarg_sigma)


    if(tolower(object$.args$noise)%in% c(1,"ar1","ar(1)") ){
      posterior_phi = inlafit$marginals.hyperpar$`Rho for idy`; zmarg_phi = inla.zmarginal(posterior_phi,silent=TRUE)
      object$fitting$hyperparameters$posteriors$phi = posterior_phi
      object$fitting$hyperparameters$results$phi = zmarg_phi

    }else if(object$.args$noise %in% c(2,"ar2","ar(2)") ){
      hypersamples = inla.hyperpar.sample(50000,inlafit)

      p=2
      phii = hypersamples[, 2L:(2L+(p-1L))]
      phis = apply(phii, 1L, inla.ar.pacf2phi)
      posterior_phi1 = cbind(density(phis[1,])$x,density(phis[1,])$y);colnames(posterior_phi1)=c("x","y"); zmarg_phi1 = inla.zmarginal(posterior_phi1,silent=TRUE)
      posterior_phi2 = cbind(density(phis[2,])$x,density(phis[2,])$y);colnames(posterior_phi2)=c("x","y"); zmarg_phi2=inla.zmarginal(posterior_phi2,silent=TRUE)
      object$fitting$hyperparameters$posteriors$phi1 = posterior_phi1
      object$fitting$hyperparameters$results$phi1 = zmarg_phi1
      object$fitting$hyperparameters$posteriors$phi2 = posterior_phi2
      object$fitting$hyperparameters$results$phi2 = zmarg_phi2
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
