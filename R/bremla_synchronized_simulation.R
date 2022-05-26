#' Simulate synchronized chronologies
#' Produces simulations from the synchronized time scale given tie-point samples.
#'
#' @param object List object which is the output of function \code{\link{bremla_chronology_simulation}}
#' @param nsims Integer denoting how many simulations should be produced.
#' @param tie_locations Numeric that provides the fixed locations for which the tie-points are located.
#' @param tie_locations_type string denoting which type of locations. Can be "depth", "age" or "index".
#' @param tiepointsims data.frame which gives the tiepoint-samples to which the models should be synchronized.
#' \code{nrow} should correspond to \code{nsims} and \code{ncol} should correspond to the number of tie-points.
#' For fixed tie-points you can just repeat each tie-point \code{nsims} times.
#' @param ncores The number of cores to use in simulating from inla.posterior.sample
#'
#' @return Returns the \code{object} list from the input and appends a list which includes all synchronized simulations
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
#' object = tiepointsimmer(object,nsims=nsims,method="adolphi")
#' object = bremla_synchronized_simulation(object,nsims=nsims)
#' object = bremla_simulationsummarizer(object,CI.type="hpd",sync=TRUE)
#' plot(object)
#' }
#'
#'
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{bremla_chronology_simulation}} \code{\link{bremla}}
#' @keywords simulation synchronization tiepoint
#'
#' @export
bremla_synchronized_simulation = function(object,nsims=10000, tie_locations=NULL,
                                          tie_locations_type="depth",tiepointsims=NULL,
                                          ncores=2){

  if(!is.null(object$tie_points)){
    tie_locations=object$tie_points$locations
    tie_locations_type=object$tie_points$locations.type
    tiepointsims = object$tie_points$samples
    if(object$tie_points$nsims != nsims) stop("Number of samples given does not correspond to 'nsims'...")
  }else{
    object$tie_points = list(samples=tiepointsims,locations=tie_locations,
                             locations.type=tie_locations_type,method="precomputed",nsims=nsims)
  }

  if(is.null(tiepointsims)&& is.null(object$tie_points)) stop("No tie-points specified")



  n = length(object$data$y)

  if(tie_locations_type == "age"){
    tie_indexes = which.index(tie_locations,object$data$y)
  }else if(tie_locations_type == "depth"){
    tie_indexes = which.index(tie_locations,object$data$z)
  }else if(tie_locations_type == "index"){
    tie_indexes = tie_locations
  }else{
    stop("Improper 'tie_locations_type' specified")
  }
  tiepointsims = tiepointsims[,!is.na(tie_indexes)]
  tie_indexes = unique(tie_indexes[!is.na(tie_indexes)])
  m = length(tie_indexes)

  free_n = n-m
  free_indexes = (1:n)[-tie_indexes]

  latentselection = list()
  reg.model = object$.args$reg.model
  if(reg.model$const) latentselection$`(Intercept)`=1
  if(reg.model$depth1) latentselection$z=1
  if(reg.model$depth2) latentselection$z2=1
  if(reg.model$proxy) latentselection$x = 1

  for(i in 2:object$.args$nevents){
    if(reg.model$psi0) latentselection[[paste0("a",i-1)]] = 1
    if(reg.model$psi1)latentselection[[paste0("c",i-1)]] = 1
  }
  ncores_postsamp = max(1,ncores)
  latentsamples = inla.posterior.sample(nsims,object$fitting$fit,
                                        selection=latentselection,verbose=FALSE,
                                        add.names=FALSE,num.threads = ncores_postsamp)

  samples = matrix(NA,nrow=n,ncol=nsims)


  time.start = Sys.time()
  for(r in 1:nsims){
    sigma_sample = object$simulation$sigma[r]
    phi_sample = object$simulation$phi[r]

    if(object$.args$noise=="ar1"){
      Qfull = Qmaker_ar1cum(n,sigma_sample,phi_sample)
    }else{
      # 'noise' is the precision matrix of the layer differences Q_x
      Qfull = Qymaker(object$.args$noise)
    }

    Qa = Qfull[-tie_indexes,-tie_indexes]
    Qab = Qfull[-tie_indexes,tie_indexes]
    if(m==1){
      Qab = as.matrix(Qab,ncol=1)
    }


    coefs = latentsamples[[r]]$latent
    dmeansim = meanmaker( coefs, reg.model, nevents=object$.args$nevents,data = object$data )

    y_mu = cumsum(c(object$preceeding$y0,dmeansim))[2:(n+1)]
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
    if((r %% 1000) == 0){
      cat("Synchronous age simulation ",r,"/",nsims,". Elapsed time: ",difftime(Sys.time(),time.start,units="secs")[[1]]," seconds...","\n",sep="")
    }
  }

  time.end = Sys.time()-time.start

  object$simulation$age_sync = samples
  object$tie_points$free_indexes=free_indexes
  object$tie_points$tie_indexes=tie_indexes
  object$tie_points$free_n=length(free_indexes)
  object$tie_points$tie_n=length(tie_indexes)

  object$time$tiepoints = time.end

  return(object)
}
