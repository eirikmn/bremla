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
#' if(inlaloader()){
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
#' depth2 = depth^2/depth[1]^2 #normalize for stability
#'
#'
#' formula = dy~-1+depth2 + proxy
#' data = data.frame(age=age,dy=dy,proxy=proxy,depth=depth,depth2=depth2)
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
#' object = bremla_chronology_simulation(object, print.progress=TRUE)
#' summary(object)
#' plot(object)
#' }
#'}
#'
#' @export
#' @importFrom stats acf arima arima.sim as.formula rnorm
#' @importFrom parallel detectCores
bremla_chronology_simulation = function(object, control.sim,print.progress=FALSE){

  if(missing(control.sim)){

    if(!is.null(object$.args$control.sim)){
      if(print.progress){
        cat("'control.sim' missing. Importing information from 'object'.\n",sep="")
      }
      control.sim = object$.args$control.sim
    }else{
      stop("Could not find 'control.sim'. Stopping.")
    }
  }
  #if(!is.null(control.sim))
   control.sim = set.options(control.sim,control.sim.default())

   if(control.sim$synchronized==TRUE) control.sim$synchronized=2
  object$.args$control.sim = control.sim

  ## sample hyperparameters
  if(is.null(object$fitting)) stop("Fitting results not found. Run 'bremla_modelfitter' first.\n")

  nsims = object$.args$control.sim$nsims
  method = object$.args$control.fit$method
  noise = object$.args$control.fit$noise

  time.start = Sys.time()


  ## sample mean ("fixed") vector and ("stochastic") noise component

  if(tolower(method) == "inla"){
    if(print.progress) cat("Simulating mean vector from fitted coefficients...",sep="")
    time.startmean = Sys.time()

    ##sample fixed parameters first
    latentselection = object$.internal$lat.selection
    reg.model = object$.internal$lat.selection

    ncores_postsamp = max(1,object$.args$control.sim$ncores)



    latentsamples = INLA::inla.posterior.sample(nsims,object$fitting$inla$fit,
                                          selection=latentselection,verbose=FALSE,
                                          add.names=FALSE,num.threads = ncores_postsamp)

    int_hyper = data.frame(matrix(NA, nrow=nsims,ncol=length(latentsamples[[1]]$hyperpar))) #get hyperparameters from joint distribution

    colnames(int_hyper) = names(latentsamples[[1]]$hyperpar)

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
      int_hyper = latentsamples[[i]]$hyperpar

      dmeansim = meanmaker( coefs, reg.model, data = object$data )



      if(control.sim$store.everything) object$simulation$dmean[,i] = dmeansim #store mean if we want

      ##sample noise component
      if(tolower(noise) %in% c("rgeneric","custom")){

        theta = int_hyper
        param.names = object$.args$control.fit$rgeneric$param.names
        hypersamples = int_hyper
        for(i in 1:length(object$.args$control.fit$rgeneric$from.theta)){
          paramsamp = object$.args$control.fit$rgeneric$from.theta[[i]](hypersamples[i])
          if(is.null(param.names[i]) || is.na(param.names[i])){
            tempname = paste0("hyperparameter",i)
            object$simulation$params[[tempname]] = paramsamp
          }else{
            object$simulation$params[[param.names[i] ]] = paramsamp
          }
        }


        Q = object$.args$control.fit$rgeneric$Q(theta,n,ntheta=length(theta))
        muvek = numeric(n)
        noisesim = Qsimmer(1, Q, muvek)
        #hypersamples = int_hyper
      }else{
        hypersamples = int_hyper

        if(tolower(noise) %in% c(0,"iid","independent","ar0","ar(0)")){
          object$simulation$params$sigma[i] = 1/sqrt(hypersamples[1])
          #object$simulation$params = list(sigma = 1/sqrt(hypersamples[,1]))

          noisesim = rnorm(n,mean=0,sd=object$simulation$params$sigma[i])


        }else if(tolower(noise) %in% c(1,"ar1","ar(1)")){
          hypersamples = int_hyper
          object$simulation$params$sigma[i] = 1/sqrt(hypersamples[1])
          object$simulation$params$phi[i] = hypersamples[2]

          noisesim = arima.sim(n=n,list(ar=c(object$simulation$params$phi[i])),
                               sd = object$simulation$params$sigma[i]*sqrt(1-object$simulation$params$phi[i]^2))


        }else if(tolower(noise) %in% c(2,"ar2","ar(2)")){
          p=2
          hypersamplesar2 = INLA::inla.hyperpar.sample(nsims,object$fitting$inla$fit)
          phii = hypersamplesar2[, 2L:(2L+(p-1L))]
          phis = apply(phii, 1L, INLA::inla.ar.pacf2phi)
          object$simulation$params$sigma[i] = 1/sqrt(hypersamples[1])
          object$simulation$params$phi1[i] = phis[1]
          object$simulation$params$phi2[i] = phis[2]


          gamma0 = (1-phis[2,i])/((1+phis[2,i])*(1-phis[1,i]-phis[2,i])*(1+phis[1,i]-phis[2,i]))

          noisesim = arima.sim(n = n, list(ar = c( phis[1,i],phis[2,i])),
                               sd = object$simulation$params$sigma[i]*sqrt(1/gamma0))


          #object$simulation$params = list(sigma = 1/sqrt(hypersamples[,1]),phi1=phis[1,],phi2=phis[2,])
        }

      }



      ## Take cumulatives. If log transformation is used, transform back first
      if(object$.args$control.fit$transform == "log"){
        object$simulation$age[,i] = object$initials$age + cumsum(exp(dmeansim+noisesim))
      }else{
        object$simulation$age[,i] = object$initials$age + cumsum(dmeansim+noisesim)
      }



    #if(print.progress) cat("Simulating ",nsims, " hyperparameters from INLA posterior...",sep="")

    # if(tolower(noise) %in% c("rgeneric","custom")){
      #hypersamples = INLA::inla.hyperpar.sample(nsims,object$fitting$inla$fit)




    # }else if(tolower(noise) %in% c(0,"ar(0)","ar0","iid","independent")){
    #   hypersamples = INLA::inla.hyperpar.sample(nsims,object$fitting$inla$fit)
    #   object$simulation = list(sigma = 1/sqrt(hypersamples[,1]))
    # }else if (tolower(noise) %in% c(1,"ar1","ar(1)")){
    #   hypersamples = INLA::inla.hyperpar.sample(nsims,object$fitting$inla$fit)
    #   object$simulation = list(sigma = 1/sqrt(hypersamples[,1]), phi=hypersamples[,2])
    #
    # }else if (tolower(noise) %in% c(2,"ar2","ar(2)")){
    #   hypersamples = INLA::inla.hyperpar.sample(nsims,object$fitting$inla$fit)
    #   p=2
    #   hypersamplesar2 = INLA::inla.hyperpar.sample(nsims,object$fitting$inla$fit)
    #   phii = hypersamplesar2[, 2L:(2L+(p-1L))]
    #   phis = apply(phii, 1L, INLA::inla.ar.pacf2phi)
    #   object$simulation = list(sigma = 1/sqrt(hypersamples[,1]),phi1=phis[1,],phi2=phis[2,])
    # }
    }
  }

  object$.args$sim = list(unsync=TRUE)
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
