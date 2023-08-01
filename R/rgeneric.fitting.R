#' rgeneric model specification for fitting custom models
#'
#' Non-standard models in INLA requires specification using the rgeneric modeling framework.
#' This includes key functions that provide information on precision matrix, mean vector, priors, graph etc.
#' The example belows shows how this can be used, and how it is supported within the 'bremla' framework.
#' The necessary arguments for use in the bremla framework is specified under \code{rgeneric}
#' in \code{\link{control.fit.default}}.
#' See example below on how rgeneric  can be used.
#'
#' @param cmd Vector containing list of function names necessary for the rgeneric model.
#' @param theta Vector describing the hyperparameters in internal scaling.
#'
#' @return When used an an input argument for \code{inla.rgeneric.define} this will return a model eligible to be used within the R-INLA framework (see example).
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{linrampfitter}}
#' @keywords bremla rgeneric inla
#'
#' @examples
#' \donttest{
#' if(inlaloader()){
#' require(INLA)
#' set.seed(1)
#'
#' # Gaussian white noise model using rgeneric and bremla
#'
#' # The from.theta function is used to transform hyperparameters from its internal scaling 'theta'
#' #   to the more natural scaling. This is used to compute posteriors and is not used in
#' #   the INLA fitting procedure. This must be specified as a list containing functions for each
#' #   hyperparameter, i.e. length(from.theta) must correspond to length(theta). This function
#' #   cannot be used by the rgeneric functions 'Q' and 'log.prior'.
#' from.theta = list(function(x)1/sqrt( exp(x) )) # convert from log(precision) -> standard deviation
#'
#' # The Q function computes the precision matrix. This is used in fitting the rgeneric model
#' #   with INLA as well as for computing simulations via precision matrix. Must be a function
#' #   of n, theta and ntheta event if they are not used in the function.
#' Q = function(theta,n,ntheta){# iid precision matrix
#'   sigma = 1/sqrt(theta[1])
#'   ii = 1:n
#'   jj = 1:n
#'   xx = 1/sigma^2*rep(1,n)
#'   return(Matrix::sparseMatrix(
#'     i=ii, j=jj, x=xx, symmetric=TRUE
#'   ))
#' }
#' # The log.prior function computes the logarithm of the joint prior of the hyperparameters in
#' #   the internal parameterization (theta). Must be a function of n, theta and ntheta even if
#' #   they are not used in the function.
#' log.prior = function(theta,n,ntheta){# gamma prior on precision (inverse variance)
#'   kappa = exp(theta[1])
#'   lprior = dgamma(kappa,shape=1.0,rate=0.02,log=TRUE) + log(kappa)
#'   return(lprior)
#' }
#'
#'
#' # R-INLA framework example
#' n = 300
#' depth = 1:n
#' ntheta = 1
#' dy = rnorm(n)+depth/40
#' model.rgeneric = inla.rgeneric.define(rgeneric.fitting, n=n,ntheta = 1,
#'                                       Q = Q,
#'                                       log.prior=log.prior)
#'
#' formula = dy ~ -1+depth+ f(idy, model=model.rgeneric)
#'
#' result = inla(formula,family="gaussian", data=data.frame(dy=dy,depth=depth,idy=1:n),
#'               num.threads = 2,
#'               control.family = list(hyper = list(prec = list(initial = 12, fixed=TRUE))) )
#'
#'
#'
#' ## bremla framework example
#'
#' y = cumsum(dy)
#'
#' formula = dy ~ -1+depth
#' data = data.frame(age=y,dy=dy,depth=depth)
#' data = rbind(c(0,0,0),data) #First row is only used to extract y0 and z0.
#'
#' control.fit = list(
#'     noise = "rgeneric", ncores=2,
#'     rgeneric = list(
#'          param.names = "std",
#'          from.theta=from.theta,
#'          Q = Q,
#'          log.prior = log.prior
#'          )
#'  )
#'  control.sim = list(synchronized=TRUE,nsims=5000)
#'
#' synchronization=list(locations=depth[c(20,50,150)],method="gauss",
#'     params=list(mean=c(y[c(20,50,150)]+c(3,-10,5)),
#'                 sd=c(5,2,10) )
#'     )
#'
#'
#' object = bremla(formula,data,nsims=5000,
#'                      control.fit=control.fit,
#'                      control.sim=control.sim,
#'                      synchronization=synchronization)
#' summary(object)
#' plot(object)
#'
#'
#'}
#'}
#'
#'
#' @export
#' @import Matrix
#' @importFrom stats dgamma dnorm
#' @importFrom Matrix sparseMatrix
rgeneric.fitting = function( #specifies necessary functions for INLA to define the linear ramp model
  cmd = c("graph", "Q","mu", "initial", "log.norm.const", "log.prior", "quit"),
  theta = NULL)
{

  # envir = environment(sys.call()[[1]])
  envir = parent.env(environment())

  if(!is.null(envir)){
    n=get("n",envir)
    ntheta=get("ntheta",envir)
  }

  if(!is.null(envir)){
    Q=get("Q",envir)
    log.prior=get("log.prior",envir)

  }
  mu = function(theta,n,ntheta){
    return(numeric(0))
  }

  graph = function(theta,n,ntheta){ #graphs of conditional dependence structure. 1 where Q[i,j] != 0, 0 where Q[i,j] = 0
    if(!is.null(envir)){
      n=get("n",envir)
      ntheta=get("ntheta",envir)
    }
    G = Q(theta,n,ntheta)
    G[G != 0] = 1
    return (G)
  }

  log.norm.const = function(theta,n,ntheta){ #INLA computes this automatically
    return(numeric(0))
  }


  initial = function(theta,n,ntheta){
    if(!is.null(envir)){
      ntheta=get("ntheta",envir)

    }
    ini = numeric(ntheta)
    return (ini)
  }

  quit = function(theta,n,ntheta){
    return ()
  }


  if (!length(theta)) theta = initial(theta,n,ntheta)

  cmd = match.arg(cmd)
  val = do.call(cmd, args = list(theta,n,ntheta))
  return (val)
}
