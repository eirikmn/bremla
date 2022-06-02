#' Bremla chronology simulation
#'
#' Simulates chronologies based on the fitted regression model.
#'
#' @param object List object which is the output of function \code{\link{bremla_modelfitter}}
#' @param control.sim List object containing specifications for simulation procedure and
#' what is to be computed. See \code{\link{control.sim.default}} for details.
#' @param print.progress Boolean. If \code{TRUE} progress will be printed to screen
#'
#' @return Returns the same \code{object} list from the input, but appends simulated
#' chronologies, mean vectors (if \code{store.means=TRUE}) and further information, summary and statistics in
#' \code{object\$simulation}.
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{bremla_prepare},\link{bremla_modelfitter}}
#' @keywords bremla simulation
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
#' nsims = 5000
#' object = bremla_prepare(age,depth,proxy,events=eventdepths,nsims=nsims)
#' object = bremla_modelfitter(object)
#' object = bremla_chronology_simulation(object,nsims=nsims)
#' plot(object,plot.inlasims = list(nsims=30,legend=NULL,xrev=FALSE,label=NULL))
#' }
#'
#'
#' @export
#' @import INLA
#' @importFrom INLA inla.hyperpar.sample inla.tmarginal inla.zmarginal inla.ar.pacf2phi
#' @importFrom stats acf arima arima.sim as.formula rnorm
#' @importFrom parallel detectCores
bremla_chronology_simulation = function(object, control.sim,print.progress=FALSE){

  if(missing(control.sim)){

    if(!is.null(object$.args$control.sim)){
      if(print.progress){
        cat("'control.sim' missing. Importing information from 'object'.",sep="")
      }
      control.sim = object$.args$control.sim
    }else{
      stop("Could not find 'control.sim'. Stopping.")
    }
  }
  #if(!is.null(control.sim))
   control.sim = set.options(control.sim,control.sim.default())

  object$.args$control.sim = control.sim

  ## sample hyperparameters
  if(is.null(object$fitting)) stop("Fitting results not found. Run 'bremla_modelfitter' first.")

  nsims = object$.args$control.sim$nsims
  method = object$.args$control.fit$method
  noise = object$.args$control.fit$noise

  time.start = Sys.time()

  if(tolower(method) == "inla"){ #if INLA is used
    reg.model = object$.internal$lat.selection

    if(print.progress) cat("Simulating ",nsims, " hyperparameters from INLA posterior...",sep="")


    if(tolower(noise) %in% c(0,"ar(0)","ar0","iid","independent")){
      hypersamples = inla.hyperpar.sample(nsims,object$fitting$inla$fit)
      object$simulation = list(sigma = 1/sqrt(hypersamples[,1]))
    }else if (tolower(noise) %in% c(1,"ar1","ar(1)")){
      hypersamples = inla.hyperpar.sample(nsims,object$fitting$inla$fit)
      object$simulation = list(sigma = 1/sqrt(hypersamples[,1]), phi=hypersamples[,2])

    }else if (tolower(noise) %in% c(2,"ar2","ar(2)")){
      hypersamples = inla.hyperpar.sample(nsims,object$fitting$inla$fit)
      p=2
      hypersamplesar2 = inla.hyperpar.sample(nsims,object$fitting$inla$fit)
      phii = hypersamplesar2[, 2L:(2L+(p-1L))]
      phis = apply(phii, 1L, inla.ar.pacf2phi)
      object$simulation = list(sigma = 1/sqrt(hypersamples[,1]),phi1=phis[1,],phi2=phis[2,])
    }


    if(print.progress) cat(" completed!\n",sep="")
  }


  ## sample mean ("fixed") vector and ("stochastic") noise component

  if(tolower(method) == "inla"){
    if(print.progress) cat("Simulating mean vector from fitted coefficients...",sep="")
    time.startmean = Sys.time()

    ##sample fixed parameters first
    latentselection = object$.internal$lat.selection
    # latentselection = list()
    # if(reg.model$const) latentselection$`(Intercept)`=1
    # if(reg.model$depth1) latentselection$z=1
    # if(reg.model$depth2) latentselection$z2=1
    # if(reg.model$proxy) latentselection$x = 1

    # for(i in 2:object$.args$nevents){
    #   if(reg.model$psi0) latentselection[[paste0("a",i-1)]] = 1
    #   if(reg.model$psi1)latentselection[[paste0("c",i-1)]] = 1
    # }
    #latentsamples = inla.posterior.sample(nsims,object$fitting$fit,selection=latentselection,verbose=FALSE,add.names=FALSE)
    ncores_postsamp = max(1,object$.args$control.sim$ncores)
    latentsamples = inla.posterior.sample(nsims,object$fitting$inla$fit,
                                          selection=latentselection,verbose=FALSE,
                                          add.names=FALSE,num.threads = ncores_postsamp)

    n=nrow(object$data)
    if(object$.args$control.sim$store.everything) object$simulation$dmean = matrix(NA,nrow=n,ncol=nsims)
    time.endmean = Sys.time()
    if(print.progress) cat(" completed in ",difftime(time.endmean,time.startmean,units="secs")[[1]]," seconds!\n",sep="")

    object$simulation$age = matrix(NA,nrow=n,ncol=nsims)

    if(print.progress) cat("Simulating chronologies...\n",sep="")
    time.startage=Sys.time()



    for(i in 1:nsims){
      if(i %% 1000 == 0 && print.progress){
        cat("Age simulation ",i,"/",nsims,". Elapsed time: ",difftime(Sys.time(),time.startage,units="secs")[[1]]," seconds...","\n",sep="")
      }
      ## from fixed parameters compute fixed model component using 'meanmaker' function
      coefs = latentsamples[[i]]$latent

      dmeansim = meanmaker( coefs, reg.model, data = object$data )



      if(control.sim$store.everything) object$simulation$dmean[,i] = dmeansim #store mean if we want

      ##sample noise component
      if(tolower(noise) %in% c(0,"iid","independent","ar0","ar(0)")){
        noisesim = rnorm(n,mean=0,sd=object$simulation$sigma[i])

      }else if(tolower(noise) %in% c(1,"ar1","ar(1)")){
        noisesim = arima.sim(n=n,list(ar=c(object$simulation$phi[i])),
                             sd = object$simulation$sigma[i]*sqrt(1-object$simulation$phi[i]^2))

      }else if(tolower(noise) %in% c(2,"ar2","ar(2)")){
        gamma0 = (1-phis[2,i])/((1+phis[2,i])*(1-phis[1,i]-phis[2,i])*(1+phis[1,i]-phis[2,i]))

        noisesim = arima.sim(n = n, list(ar = c( phis[1,i],phis[2,i])),
                             sd = object$simulation$sigma[i]*sqrt(1/gamma0))
      }


      ## Take cumulatives. If log transformation is used, transform back first
      if(object$.args$control.fit$transform == "log"){
        object$simulation$age[,i] = object$initials$age + cumsum(exp(dmeansim+noisesim))
      }else{
        object$simulation$age[,i] = object$initials$age + cumsum(dmeansim+noisesim)
      }


    }
  }

  time.endage = Sys.time()
  if(print.progress) cat("Completed in ",difftime(time.endage,time.startage,units="secs")[[1]]," seconds!\n",sep="")
  object$time$simulation$mean = difftime(time.endmean,time.startmean,"secs")[[1]]
  object$time$simulation$age = difftime(time.endage,time.startage,"secs")[[1]]
  object$time$simulation$total = difftime(time.endage,time.start,units="secs")[[1]]

  if(object$.args$control.sim$summary$compute){
    object = bremla_simulationsummarizer(object,sync=FALSE,print.progress=print.progress)
  }

  return(object)
}
