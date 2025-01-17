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
#' formula = dy~-1+depth2 + proxy
#' data = data.frame(age=age,dy=dy,proxy=proxy,depth=depth,depth2=depth2)
#' data = rbind(c(y0,NA,NA,z0,NA),data) #First row is only used to extract y0 and z0.
#'
#' events=list(locations=c(1210,1220,1240))
#' control.fit = list(ncores=2,noise="ar1")
#' synchronization=list(locations=depth[c(100,400,700)],method="gauss",
#'                            params=list(mean=c(age[c(100,400,700)]+c(30,-100,50)),
#'                                        sd=c(50,20,100)
#'                                        )
#'                      )
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
#' }
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{bremla_chronology_simulation}} \code{\link{bremla}}
#' @keywords simulation synchronization tiepoint
#'
#' @export
#' @importFrom utils write.table
bremla_synchronized_simulation = function(object,control.sim,print.progress=FALSE){

time.start = Sys.time()
  if(missing(control.sim)){

    if(!is.null(object$.args$control.sim)){
      if(print.progress){
        # cat("'control.sim' missing. Importing information from 'object'.",sep="")
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
    if(object$tie_points$nsims < nsims){
      stop("Number of samples to be produced (nsims) exceeds the number of tie-point samples available. Run 'tiepointsimmer' again...")
    }else{
      tiepointsims = tiepointsims[1:nsims,]
    }
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

  if(tolower(object$.args$control.fit$method) == "inla"){ # && is.null(object$simulation)){ #if INLA is used


    ## Collect unsynchronized samples



    if(print.progress){
      cat("Preparing synchronization...\n",sep="")
    }
    if(is.null(object$.args$synchronization$agedisc)){

      ## Old
      if(print.progress){
        cat("Assuming no structural biases...\n",sep="")
      }

      reg.model = object$.internal$lat.selection
      latentselection = object$.internal$lat.selection


      timecoef.start = Sys.time()
      ncores_postsamp = max(1,object$.args$control.sim$ncores)
      latentsamples = INLA::inla.posterior.sample(nsims,object$fitting$inla$fit,
                                                  selection=latentselection,verbose=FALSE,
                                                  add.names=FALSE,num.threads = ncores_postsamp)
      #int_hyper = latentsamples[[i]]

      #int_hyper = data.frame(matrix(NA, nrow=nsims,ncol=length(latentsamples[[1]]$hyperpar))) #get hyperparameters from joint distribution

      #colnames(int_hyper) = names(latentsamples[[1]]$hyperpar)
      timecoef.end = Sys.time()
      if(print.progress) cat(" completed in ",difftime(timecoef.end,timecoef.start,units="secs")[[1]]," seconds.\n",sep="")
      samples = matrix(NA,nrow=n,ncol=nsims)




      if(print.progress) cat("Sampling fixed coefficients...",sep="")


      if(print.progress) cat("\nSimulating synchronized chronologies...\n",sep="")
      timeage.start = Sys.time()
      for(r in 1:nsims){
        hypersamples = latentsamples[[r]]$hyperpar
        int_hyper = hypersamples

        if(tolower(noise) %in% c("rgeneric","custom")){

          param.names = object$.args$control.fit$rgeneric$param.names
          for(i in 1:length(object$.args$control.fit$rgeneric$from.theta)){
            paramsamp = object$.args$control.fit$rgeneric$from.theta[[i]](hypersamples[i])
            if(is.null(param.names[i]) || is.na(param.names[i])){
              tempname = paste0("hyperparameter",i)
              object$simulation$params[[tempname]][i] = paramsamp
            }else{
              object$simulation$params[[param.names[i] ]][i] = paramsamp
            }
          }

          theta = hypersamples

          Qx = object$.args$control.fit$rgeneric$Q(theta,n,ntheta=length(theta))

          Qfull = Qymaker(Qx)

        }else if( tolower(noise) %in% c(1,"ar1","ar(1)") ){

          hypersamples = int_hyper
          #hypersamples = INLA::inla.hyperpar.sample(nsims,object$fitting$inla$fit)
          object$simulation$params$sigma[r] = 1/sqrt(hypersamples[1])
          object$simulation$params$phi[r] = hypersamples[2]

          sigma_sample = object$simulation$params$sigma[r]
          phi_sample = object$simulation$params$phi[r]

          if(object$.args$control.fit$noise=="ar1"){
            Qfull = Qmaker_ar1cum(n,sigma_sample,phi_sample)
          }else{
            # 'noise' is the precision matrix of the layer differences Q_x
            Qfull = Qymaker(object$.args$control.fit$noise)
          }


        }else{
          stop("Currently, only the 'ar1' noise model is implemented. Other models can be specified via the 'rgeneric' model, see '?rgeneric.fitting' for details.")
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
        b_temp = (-Qab%*%(tiepointsims[r,]-tie_mu))
        if(!is.null(dim(b_temp))){
          b_temp = b_temp[,1]
        }
        w = solve(La,b_temp)
        mu_temp = solve(t(La),w)
        if(!is.null(dim(mu_temp))){
          mu_temp = mu_temp[,1]
        }
        z0 = rnorm(free_n)
        v = solve(t(La),z0)
        if(!is.null(dim(v))){
          v = v[,1]
        }

        samples[free_indexes,r] = mu_temp + free_mu + v
        samples[tie_indexes,r] = tiepointsims[r,]

        #mu_amidb = (free_mu - solve(Qa)%*%Qab%*%(skewsamples[r]-tie_mu))[,1]

        #samples[free_indexes,r] = Qsimmer(1,Qa,mu_amidb)
        if(print.progress && (r %% 1000) == 0){
          cat("Synchronous age simulation ",r,"/",nsims,". Elapsed time: ",difftime(Sys.time(),timeage.start,units="secs")[[1]]," seconds...","\n",sep="")
        }
      }

      object$simulation$age_sync = samples







    }else{ #define age discrepancy method







      # Set up tie-points
      n_tie = length(object$tie_points$locations_indexes)

      mt = matrix(NA,nrow=n,ncol=nsims)
      mt[object$tie_points$tie_indexes,] = t(object$tie_points$samples)

      y2mat = mt - object$simulation$age #discrepancies



      if(!is.null(object$.args$synchronization$agedisc)){
        if(is.null(object$.args$synchronization$agedisc$model)){
          object$.args$synchronization$agedisc$model = "rw2"
        }
        if(is.null(object$.args$synchronization$agedisc$method)){
          object$.args$synchronization$agedisc$method = "INLA"
        }
        if(is.null(object$.args$synchronization$agedisc$options)){
          object$.args$synchronization$agedisc$options$stepsizes = c(0.005,0.001,0.01)
          object$.args$synchronization$agedisc$options$restart.fromlast = FALSE
          if(is.null(object$.args$synchronization$agedisc$options$inla.options)){
            object$.args$synchronization$agedisc$options$inla.options = list(
              control.compute=list(config=TRUE),
              silent=1L,
              num.threads=object$.args$control.fit$ncores,
              control.family = list(hyper=list(prec=list(initial=12, fixed=TRUE)))
            )
          }
        }

        if(is.null(object$.args$synchronization$agedisc$hyperprior)){
          if(object$.args$synchronization$agedisc$model %in% c("rw2", "randomwalk2")){
            #object$.args$synchronization$agedisc$hyperprior=NULL
          }else{
            stop("No hyperprior for age discrepancy model given, and no default value available!")
          }
        }
      }




      if(tolower(object$.args$synchronization$agedisc$method) != "inla"){
        stop("Currently, only INLA method is available for age discrepancies")
      }
      # Assigning default values


      disc.model=object$.args$synchronization$agedisc$model

      if(print.progress){
        cat("Assuming age discrepancy is described by a ", disc.model, " model.\n",sep="")
      }
      if(tolower(disc.model) %in% c("rw1","randomwalk1")){
        if(is.null(object$.args$synchronization$agedisc$hyperprior)){
          object$.args$synchronization$agedisc$hyperprior=list(prec=list(param=c(10,5e-05)))
        }
        idx2 = seq(1,1000,length.out=n)
        control.mode =list(restart=TRUE, theta=-2)
        formula2 = y2 ~ -1 + f(idx2,model="rw1", values=idx2, constr=FALSE,
                               hyper=object$.args$synchronization$agedisc$hyperprior)
      }else if(tolower(disc.model) %in% c("rw2","randomwalk2")){
        idx2 = seq(1,1000,length.out=n)
        if(is.null(object$.args$synchronization$agedisc$hyperprior)){
          object$.args$synchronization$agedisc$hyperprior=list(prec=list(param=c(10,5e-05)))
        }
        #control.mode =list(restart=TRUE, theta=-2)
        control.mode =NULL
        formula2 = y2 ~ -1 + f(idx2,model="rw2", values=idx2, constr=FALSE,
                               hyper=object$.args$synchronization$agedisc$hyperprior)
      }else{
        idx2 = 1:n
        control.mode =list(restart=TRUE, theta=c(2,2))
        if(is.null(object$.args$synchronization$agedisc$hyperprior)){
          object$.args$synchronization$agedisc$hyperprior =
            list(prec=list(param=c(8, 5e-05),initial=2),
                 rho=list(param=c(2,0.15), initial=2))
        }
        formula2 = y2 ~ -1 + f(idx2,model="ar1",
                               hyper=object$.args$synchronization$agedisc$hyperprior
                                )
      }

      object$age_discrepancy = object$.args$synchronization$agedisc
      object$age_discrepancy$formula=formula2
      object$age_discrepancy$options$inla.options$control.mode=control.mode

      x2sims = matrix(NA,nrow=n,ncol=nsims)



      time.agedisc.start = Sys.time()

      for(i in 1:nsims){


        elapsed.time = Sys.time()-time.agedisc.start
        if(i==1 && print.progress){
            cat("Starting fitting the model with INLA and simulating from posterior predictor.\n",sep="")
        }
        #if( (i %% round(0.1*nsims) == 0) && print.progress ){
        if( print.progress ){
        #if(i %% 1000 == 0 && print.progress){
          cat("Fitting age discrepancy sample #",i," / ",nsims, " - elapsed time: ",round(difftime(Sys.time(),time.agedisc.start,units="secs")[[1]], digits=1)," seconds...\n",sep="")
        }


        if(i>1 && object$.args$synchronization$agedisc$options$restart.fromlast){
          #object$age_discrepancy$options$control.mode
          control.mode$result=r2
          control.mode$restart=TRUE
        }else{
          control.mode$result=NULL
          control.mode = object$age_discrepancy$options$inla.options$control.mode
        }
        inla.options = object$age_discrepancy$options$inla.options
        inla.options$control.mode = control.mode


        safe_function <- function(stepHs, inla.options){
          success = FALSE
          #stepHs = c(0.005, 0.001, 0.01,0.002,0.003,0.008)

          attempt = 1
          while(!success){
            tryCatch({
              if(attempt > length(stepHs)){
                stepH = rnorm(1,mean=0.005,sd=0.02)
              }else{
                stepH = stepHs[attempt]
              }

              #cat("attempt: ",attempt,"\n",sep="")
              if(attempt >7){
                success <- TRUE #break to avoid infinite loop
                return(NA)
              }
              inla.options$control.inla$h = stepH
              inla.options$control.compute$config=TRUE #this needs to be set to TRUE



              r2 = do.call(INLA::inla, c(list(formula = formula2,
                                              data = data.frame(y2=y2mat[,i], idx2=idx2)),
                                         inla.options) )

              success <- TRUE
              if(attempt>1){
                cat("Function succeeded after ",attempt," attemps.\n",sep="")
              }

              return(r2)
            },error = function(e){
              #message("error",e$message)
              #res = inla.climate(yw,fw, model=model, inla.options=inla.options, stepLength=0.001)
              message("Attempt ", attempt, " failed with error: ", e$message)
              attempt <<- attempt + 1

              #return(res)
            })
          }

        }

        # Doing some stuff as recommended by INLA, to prevent errors
        if(tolower(disc.model) %in% c("rw2", "randomwalk2")){
          invisible(INLA::inla.models())
          m = get("inla.models", INLA::inla.get.inlaEnv())
          m.old = m
          m$latent$rw2$min.diff = NULL
          assign("inla.models", m, INLA::inla.get.inlaEnv())
        }

        r2 = safe_function(object$.args$synchronization$agedisc$options$stepsizes, inla.options = inla.options)

        # Reverting
        if(tolower(disc.model) %in% c("rw2", "randomwalk2")){
          m$latent$rw2$min.diff = m.old$latent$crw2$min.diff
          assign("inla.models", m, INLA::inla.get.inlaEnv())
        }


        r2s = INLA::inla.posterior.sample(1,r2, selection=list(Predictor=1:n))

        x2sims[,i] = r2s[[1]]$latent

        if(is.character(object$.args$synchronization$agedisc$options$backupfilename)){
          # backup if filename is specified
           write.table(x2sims + object$simulation$age, file=object$.args$synchronization$agedisc$options$backupfilename)
        }

      }
      if(print.progress){

      }


      # chronsims = object$simulation$age + x2sims


      object$simulation$age_sync = object$simulation$age + x2sims


      #compute mean and credible intervals

      # object = bremla_simulationsummarizer(object,sync=TRUE,print.progress=print.progress)
#
#       mean=numeric(n)
#       lower=numeric(n)
#       upper=numeric(n)
#       for(i in 1:n){
#         dens=density(chronsims[i,]); dens=data.frame(x=dens$x,y=dens$y)
#         zm = inla.zmarginal(dens,silent=TRUE)
#         mean[i] = zm$mean
#         lower[i] = zm$quant0.025
#         upper[i] = zm$quant0.975
#       }


    }

#     if(FALSE){
#   plot(object$simulation$summary_sync$mean-object$data$age,
#        ylim=range(object$simulation$summary_sync$lower-object$data$age,
#                   object$simulation$summary_sync$upper-object$data$age),type="l")
#   lines(object$simulation$summary_sync$lower-object$data$age,col="red")
#   lines(object$simulation$summary_sync$upper-object$data$age,col="red")
#
#   for(ii in 1:4){
#     quants = quantile(tiepointsims[,ii], probs = c(0.025,0.975)) - object$data$age[object$tie_points$locations_indexes[ii]]
#     segments(object$tie_points$locations_indexes[ii], quants[1],
#              object$tie_points$locations_indexes[ii], quants[2])
#   }
#
#   object$tie_points$tie_indexes
#
# }





  }

  timeage.end = Sys.time()

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
