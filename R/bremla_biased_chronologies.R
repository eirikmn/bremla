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
#' require(stats)
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
#' object = bremla_biased_chronologies(object, control.bias = list(bias.model="uniform"),
#'              print.progress=TRUE)
#' summary(object)
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

  nsims = control.bias$nsims
  if(nsims > ncol(object$simulation$age)){
    warning("Number of simulated biases exceeds number of simulated chronologies! Uses lowest number.")
    nsims = min(nsims,ncol(object$simulation$age))
  }



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
    object$biases[[listr]] = list(mean = biasmean, sd = biassd,
                                  quant0.025 = biasmean-1.96*biassd,
                                  quant0.975 = biasmean+1.96*biassd
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
