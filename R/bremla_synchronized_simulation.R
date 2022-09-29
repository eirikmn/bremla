#' Simulate synchronized chronologies
#'
#' Produces simulations from the synchronized time scale given tie-point samples.
#'
#' @param object Output of function \code{\link{tiepointsimmer}}.
#' @param control.sim List containing specifications for the simulation procedure,
#' including the number of samples to be generated and whether or not the chronologies
#' should be synchronized (for this function this should be set to \code{TRUE}).
#' See \code{\link{control.sim.default}} for details.
#' @param print.progress Boolean. If \code{TRUE} progress will be printed to screen.
#'
#' @return Returns the \code{object} list from the input and appends a list which
#' includes all synchronized simulations, summary statistics and other information
#' computed or related to the simulation procedure in \code{object\$simulation}.
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
#' synchronization=list(locations=depth[c(100,400,700)],locations_unit="depth",method="gauss",
#'         params=list(mean=age[c(100,400,700)]+c(20,5,-20),sd=c(20,10,50)))
#' control.sim=list(synchronized=TRUE,
#'                  summary=list(compute=TRUE))
#'
#' object = bremla_prepare(formula,data,nsims=5000,reference.label="simulated timescale",
#'                         events = events,
#'                         synchronization=synchronization,
#'                         control.fit=control.fit,
#'                         control.sim=control.sim)
#' object = bremla_modelfitter(object)
#' object = tiepointsimmer(object)
#' object = bremla_synchronized_simulation(object, print.progress=TRUE)
#' summary(object)
#' plot(object)
#' }
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{bremla_chronology_simulation}} \code{\link{bremla}}
#' @keywords simulation synchronization tiepoint
#'
#' @export
bremla_synchronized_simulation = function(object,control.sim,print.progress=FALSE){

time.start = Sys.time()
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
  if(control.sim$synchronized==FALSE) control.sim$synchronized=2

  ## sample hyperparameters
  if(is.null(object$fitting)) stop("Fitting results not found. Run 'bremla_modelfitter' first.")
  if(is.null(object$tie_points)) stop("Tie-points samples not found. Run 'tiepointsimmer' first.")

  nsims = control.sim$nsims
  if(!is.null(object$tie_points)){
    tie_locations=object$tie_points$locations
    tie_locations_unit=object$tie_points$locations_unit
    tiepointsims = object$tie_points$samples
    if(object$tie_points$nsims != nsims) stop("Number of samples given does not correspond to 'nsims'...")
  }else{
    object$tie_points = list(samples=tiepointsims,locations=tie_locations,
                             locations.type=tie_locations_unit,method="precomputed",nsims=nsims)
  }

  if(is.null(tiepointsims)&& is.null(object$tie_points)) stop("No tie-points specified")



  n = nrow(object$data)
  noise = object$.args$control.fit$noise

  if(tie_locations_unit == "age"){
    tie_indexes = which.index(tie_locations,object$data$age)
  }else if(tie_locations_unit == "depth"){
    tie_indexes = which.index(tie_locations,object$data$depth)
  }else if(tie_locations_unit == "index"){
    tie_indexes = tie_locations
  }else{
    stop("Improper 'tie_locations_unit' specified")
  }
  tiepointsims = tiepointsims[,!is.na(tie_indexes)]
  tie_indexes = unique(tie_indexes[!is.na(tie_indexes)])
  m = length(tie_indexes)

  free_n = n-m
  free_indexes = (1:n)[-tie_indexes]

  if(tolower(object$.args$control.fit$method) == "inla" && is.null(object$simulation)){ #if INLA is used
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

  if(print.progress) cat("Sampling fixed coefficients...",sep="")
  latentselection = object$.internal$lat.selection
  # latentselection = list()
  # reg.model = object$.args$reg.model
  # if(reg.model$const) latentselection$`(Intercept)`=1
  # if(reg.model$depth1) latentselection$z=1
  # if(reg.model$depth2) latentselection$z2=1
  # if(reg.model$proxy) latentselection$x = 1
  #
  # for(i in 2:object$.args$nevents){
  #   if(reg.model$psi0) latentselection[[paste0("a",i-1)]] = 1
  #   if(reg.model$psi1)latentselection[[paste0("c",i-1)]] = 1
  # }
  timecoef.start = Sys.time()
  ncores_postsamp = max(1,object$.args$control.sim$ncores)
  latentsamples = inla.posterior.sample(nsims,object$fitting$inla$fit,
                                        selection=latentselection,verbose=FALSE,
                                        add.names=FALSE,num.threads = ncores_postsamp)

  timecoef.end = Sys.time()
  if(print.progress) cat(" completed in ",difftime(timecoef.end,timecoef.start,units="secs")[[1]]," seconds.\n",sep="")
  samples = matrix(NA,nrow=n,ncol=nsims)


  if(print.progress) cat("Simulating synchronized chronologies...\n",sep="")
  timeage.start = Sys.time()
  for(r in 1:nsims){
    sigma_sample = object$simulation$sigma[r]
    phi_sample = object$simulation$phi[r]

    if(object$.args$control.fit$noise=="ar1"){
      Qfull = Qmaker_ar1cum(n,sigma_sample,phi_sample)
    }else{
      # 'noise' is the precision matrix of the layer differences Q_x
      Qfull = Qymaker(object$.args$control.fit$noise)
    }

    Qa = Qfull[-tie_indexes,-tie_indexes]
    Qab = Qfull[-tie_indexes,tie_indexes]
    if(m==1){
      Qab = as.matrix(Qab,ncol=1)
    }


    coefs = latentsamples[[r]]$latent
    dmeansim = meanmaker( coefs, latentselection,data = object$data )

    y_mu = cumsum(c(object$initials$age,dmeansim))[2:(n+1)]
    free_mu = y_mu[free_indexes]
    tie_mu = y_mu[tie_indexes]


    La = t(chol(Qa))
    b_temp = (-Qab%*%(tiepointsims[r,]-tie_mu))[,1]
    w = solve(La,b_temp)
    mu_temp = solve(t(La),w)[,1]
    z0 = rnorm(free_n)
    v = solve(t(La),z0)[,1]
    samples[free_indexes,r] = mu_temp + free_mu + v
    samples[tie_indexes,r] = tiepointsims[r,]

    #mu_amidb = (free_mu - solve(Qa)%*%Qab%*%(skewsamples[r]-tie_mu))[,1]

    #samples[free_indexes,r] = Qsimmer(1,Qa,mu_amidb)
    if(print.progress && (r %% 1000 == 0)){
      cat("Synchronous age simulation ",r,"/",nsims,". Elapsed time: ",difftime(Sys.time(),timeage.start,units="secs")[[1]]," seconds...","\n",sep="")
    }
  }

  timeage.end = Sys.time()

  object$simulation$age_sync = samples
  object$tie_points$free_indexes=free_indexes
  object$tie_points$tie_indexes=tie_indexes
  object$tie_points$free_n=length(free_indexes)
  object$tie_points$tie_n=length(tie_indexes)

  #object$time$tiepoints = time.end
  if(object$.args$control.sim$synchronized == 2){
    object$time$simulation$total = object$time$simulation$total+difftime(timeage.end,time.start,units="secs")[[1]]
  }else{
    object$time$simulation$total = difftime(timeage.end,time.start,units="secs")[[1]]
  }

  object$.args$sim = list(sync=TRUE)
  if(object$.args$control.sim$summary$compute){
    object = bremla_simulationsummarizer(object,sync=TRUE,print.progress=print.progress)
  }
  return(object)
}
