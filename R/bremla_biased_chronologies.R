#' Linear ramp fit
#'
#' Fits the linear ramp described by Erhardt et al. (2019) to proxy values using INLA.
#'
#' @param object List object which is the output of function \code{\link{bremla_chronology_simulation}}
#' @param bias.model character string describing from which latent distribution the unknown (stochastic) bias factors are drawn. Currently, only \code{bias.model="uniform"} is supported.
#' @param biasparams matrix describing the parameters of the latent distribution. For the uniform distribution only 2 parameters must be specified, hence the number of rows are 2. If the number of columns are larger than 1 then repeated analyses are performed with different sets of parameters.
#' @param nsims the number of samples to be generated. Cannot exceed the number of simulated chronologies from \code{\link{bremla_chronology_simulation}}.
#' @param store.samples Boolean. If \code{store.samples=TRUE} then samples are stored. If \code{FALSE} they are not.
#'
#' @return Returns the same \code{object} list from the input, but appends another list of the summary statistics for each analysis (and samples if \code{store.samples=TRUE}), and inputs (\code{settings})
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{bremla_chronology_simulation}}
#' @keywords bremla bias
#'
#' @examples
#'
#' @export
#' @importFrom stats runif
bremla_biased_chronologies = function(object,bias.model="uniform",biasparams = c(0.99,1.01),nsims=10000,store.samples=FALSE){
  if(nsims > dim(object$simulation$age)[2]) stop("Number of simulated biases exceeds number of simulated chronologies! Shutting down...")
  n = dim(object$simulation$age)[1]
  m = dim(biasparams)[2]
  for(iter in 1:m){
    biasparam = biasparams[,iter]
    biases = runif(nsims,min=biasparam[1],max=biasparam[2])

    if(store.samples){
      biasedages = matrix(NA,nrow=n,ncol=nsims)
    }
    bias.x1=numeric(n)
    bias.x2=numeric(n)
    for(i in 1:nsims){
      sample = biases[i]*object$simulation$age[,i]
      if(store.samples){
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
    if(store.samples){
      object$biases[[listr]]$simulations = biasedages
    }
    object$biases[[listr]]$.args = list(bias.model=bias.model,biasparam=biasparam,nsims=nsims,store.samples=store.samples)
  }

  return(object)
}
