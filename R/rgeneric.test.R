
rgeneric.test = function( #specifies necessary functions for INLA to define the linear ramp model
  cmd = c("graph", "Q","mu", "initial", "log.norm.const", "log.prior", "quit"),
  theta = NULL)
{

  #envir = environment(sys.call()[[1]])
  envir = parent.env(environment())
  #TAU = exp(12)

  interpret.theta = function() { #helpful function to transform back from internal parametrization
    if(!is.null(envir)){
      dMCEE=get("dMCE",envir)
    }

    kappa = exp(theta[1])
    sigma = 1/sqrt(kappa)
    sigmas = sigma*dMCEE/2 #+1/TAU
    kappas = 1/sigmas^2
cat("sigma:",sigma,"\n")
cat("rangesigmas:",range(sigmas),"\n")
    return(list(sigma=sigma,sigmas=sigmas,kappa=kappa,kappas=kappas))
  }

  mu = function() { #mean vector defined as a linear ramp

    return(numeric(0))
  }


  graph = function(){ #graphs of conditional dependence structure. 1 where Q[i,j] != 0, 0 where Q[i,j] = 0
    G = Q()
    G[G != 0] = 1
    return (G)
  }

  Q = function(){ #inverse covariance matrix: AR(1) that allow for unequal spacing
    if(!is.null(envir)){
      #n=get("n",envir)
      #timepoints=get("timepoints",envir)
    }
    param=interpret.theta()


    kappas = param$kappas
    n=length(kappas)
    ii = c(1:(n-1), n, 1:(n-1))
    jj = c(1:(n-1), n, 2:n)

    xx = c(kappas[1:(n-1)]+kappas[2:n], kappas[n], -kappas[2:n])

    Q = Matrix::sparseMatrix(
      i = ii,
      j = jj,
      x = xx,
      symmetric = TRUE
    )
    cat("Qrange:",range(Q),"\n")
    cat("Qdet:",det(Q),"\n")
    #diag(Q)=diag(Q)
    return (Q)
  }

  log.norm.const = function(){ #INLA computes this automatically
    return(numeric(0))
  }

  log.prior = function(){
    if(!is.null(envir)){
      #tslutt=get("tslutt",envir)
      #tstart=get("tstart",envir)
      #ystart=get("ystart",envir)
      #log.theta.prior=get("log.theta.prior",envir)
      # rescale.y=get("rescale.y",envir)
    }

    lprior =  dgamma(exp(theta[1]),2,rate = 0.15,log=TRUE) + theta[1] #sigma/kappa
    cat("lprior:",lprior,"\n")

    return (lprior)
  }

  initial = function(){
    ini = c(0)
    return (ini)
  }

  quit = function(){
    return ()
  }
  if (!length(theta)) theta = initial()

  cmd = match.arg(cmd)
  val = do.call(cmd, args = list())
  return (val)
}
