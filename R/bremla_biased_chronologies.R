#' Linear ramp fit
#'
#' Fits the linear ramp described by Erhardt et al. (2019) to proxy values using INLA.
#'
#' @param object List object which is the output of function \code{\link{bremla_chronology_simulation}}
#' @param control.bias List containing specifications for how the unknown random bias
#' should be implemented. See \code{\link{control.bias.default}} for details.
#' @param print.progress Boolean. If \code{TRUE} progress will be printed to the screen.
#'
#' @return Returns the same \code{object} list from the input, but appends another
#' list of the summary statistics for each analysis (and samples if \code{store.samples=TRUE}),
#' and inputs (\code{settings}). These are stored in \code{object\$biases.}
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{bremla_chronology_simulation}}
#' @keywords bremla bias
#' @examples
#' \donttest{
#' data("event_intervals")
#' data("events_rasmussen")
#' data("NGRIP_5cm")
#' age = NGRIP_5cm$age
#' depth = NGRIP_5cm$depth
#' d18O = NGRIP_5cm$d18O
#' proxy=d18O
#' eventdepths = events_rasmussen$depth
#' eventindexes = c(1,which.index(eventdepths, depth[2:length(depth)]) )
#' eventindexes = unique(eventindexes[!is.na(eventindexes)])
#' nsims=5000
#' object = bremla_prepare(age,depth,proxy,events=eventdepths,nsims=nsims)
#' object = bremla_modelfitter(object)
#' object = bremla_chronology_simulation(object,nsims=nsims)
#' object = bremla_biased_chronologies(object,bias.model="uniform",nsims=nsims,
#'            biasparams=cbind( c(1,1),c(0.98,1.02),c(0.96,1.04) ) )
#' plot(object)
#' }
#' @export
#' @importFrom stats runif
bremla_biased_chronologies = function(object,control.bias,print.progress=FALSE){
#bias.model="uniform",biasparams = c(0.99,1.01),nsims=10000,store.samples=FALSE){


  time.start = Sys.time()



  if(missing(control.bias)){

    if(!is.null(object$.args$control.bias)){
      if(print.progress){
        cat("'control.bias' missing. Importing information from 'object'.",sep="")
      }
      control.bias = object$.args$control.bias
    }else{
      stop("Could not find 'control.bias'. Stopping.")
    }
  }
  #if(!is.null(control.bias))
  control.bias = set.options(control.bias,control.bias.default())

  object$.args$control.bias = control.bias

  ## sample hyperparameters
  if(is.null(object$linramp))
    stop("Linear ramp fit not found. Run 'linrampfitter' first!")

  nsims = control.bias$nsims
  if(nsims > ncol(object$simulation$age))
    stop("Number of simulated biases exceeds number of simulated chronologies! Shutting down...")

  n = nrow(object$simulation$age)
  biasparams = control.bias$biasparams

  if(class(biasparams) == "numeric"){
    biasparams = matrix(biasparams,ncol=1)
  }
  m = ncol(biasparams)
  for(iter in 1:m){
    biasparam = biasparams[,iter]
    biases = runif(nsims,min=biasparam[1],max=biasparam[2])

    if(control.bias$store.samples){
      biasedages = matrix(NA,nrow=n,ncol=nsims)
    }
    bias.x1=numeric(n)
    bias.x2=numeric(n)
    for(i in 1:nsims){
      sample = biases[i]*object$simulation$age[,i]
      if(control.bias$store.samples){
        biasedages[,i] = sample
      }
      bias.x1 = bias.x1 + sample
      bias.x2 = bias.x2 + sample^2
    }
    biasmean = bias.x1/nsims
    biassd = sqrt(  1/(nsims-1) * (bias.x2 - 2*bias.x1*biasmean + nsims*biasmean**2)   )

    listr = paste0("bias",iter)
    object$biases[[listr]] = list(mean = biasmean, sd = biassd, quant0.025 = biasmean-1.96*biassd,quant0.975=biasmean+1.96*biassd
    )
    if(control.bias$store.samples){
      object$biases[[listr]]$simulations = biasedages
    }

  }

  object$biases$.args = list(bias.model=control.bias$bias.model,
                             biasparam=biasparams,
                             nsims=nsims,
                             store.samples=control.bias$store.samples,
                             nbiases = m)
  time.end = Sys.time()
  time.full = difftime(time.end,time.start,units="secs")[[1]]
  object$biases$time = time.full
  object$time$biases = time.full
  return(object)
}
