#' Linear ramp fit
#'
#' Fits the linear ramp described by Erhardt et al. (2019) to proxy values using INLA.
#'
#' @param object List object which is the output of function \code{\link{bremla_chronology_simulation}}
#' @param control.linramp Includes the data and specifications for the linear ramp
#' model fitting procedure. Must include \code{control.linramp\$proxy} and \code{control.linramp\$interval}
#' (preferably also \code{control.linramp\$interval.unit}),
#' the rest can be set to default values. See \code{\link{control.linramp.default}} for details.
#' @param print.progress Boolean. If \code{TRUE} progress will be printed to screen.
#'
#' @return Returns the same \code{object} list from the input, but appends another list of results from the linear ramp fit, including posterior marginal distributions of the hyperparameters, summary statistics, and simulations of linear ramp and transition end point (if any).
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @references Erhardt et al. (2019)
#' @seealso \code{\link{bremla_chronology_simulation},\link{events_depth_to_age}}
#' @keywords bremla linear_ramp
#'
#' @examples
#'
#'\donttest{
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
#' object = bremla_prepare(age,depth,proxy,events=eventdepths,nsims=0)
#'
#' lowerints = which.index(event_intervals$depth_int_lower.m, depth[2:length(depth)])
#' upperints = which.index(event_intervals$depth_int_upper.m, depth[2:length(depth)])
#'
#' eventnumber=13 #number between 1 and 29. specifies which transition to consider
#' depth.reference = event_intervals$NGRIP_depth_m[eventnumber]
#' interval = lowerints[eventnumber]:upperints[eventnumber]
#'
#' object = linrampfitter(object,interval,label="GI-11",depth.reference=depth.reference)
#'
#' plot(object,plot.proxydata=list(age=TRUE,depth=FALSE,xrev=FALSE,label=NULL),
#'     plot.ls = NULL,
#'     plot.inla.posterior = NULL,
#'     plot.inlasims = NULL,
#'     plot.syncsims = NULL,
#'     plot.tiepoints = NULL,
#'     plot.bias = NULL,
#'     plot.linramp = list(depth.reference=NULL,show.t0=TRUE,show.t1=TRUE,xrev=TRUE,label=NULL))
#'
#'}
#'
#' @export
#' @import INLA numDeriv
#' @importFrom INLA inla.rgeneric.define inla.emarginal inla.smarginal
#' @importFrom stats optim
#' @importFrom numDeriv grad
#'
linrampfitter = function(object,control.linramp,print.progress=FALSE){
  time.start = Sys.time()


  if(missing(control.linramp)){

    if(!is.null(object$.args$control.linramp)){
      if(print.progress){
        cat("'control.linramp' missing. Importing information from 'object'.",sep="")
      }
      control.linramp = object$.args$control.linramp
    }else{
      stop("Could not find 'control.linramp'. Stopping.")
    }
  }
  #if(!is.null(control.linramp))
  control.linramp = set.options(control.linramp,control.linramp.default())

  object$.args$control.linramp = control.linramp

  ## sample hyperparameters
  if(is.null(object$simulation$age) && is.null(object$simulation$age_sync))
    stop("Chronologies not found. Run 'bremla_synchronized_simulation' or 'bremla_chronology_simulation' first!")
  interval = control.linramp$interval



  if(control.linramp$interval.what == "depth"){
    intervalrange = which.index(range(interval),object$data$depth)
    interval = intervalrange[1]:intervalrange[2]
  }else if(control.linramp$interval.what == "age"){
    intervalrange = which.index(range(interval),object$data$age)
    interval = intervalrange[1]:intervalrange[2]
  }else if(control.linramp$interval.what == "index"){
    interval = range(interval)[1]:range(interval)[2]
  }

  if(print.progress) cat("Initializing linear ramp fit using INLA.\n",sep="")
  df_event = data.frame(xx=rev(object$data$depth[interval]),
                        yy=rev(control.linramp$proxy[interval])) #reverse: Want depth axis representing old->new

  ## default initial values for optim function
  if(is.null(control.linramp$opt.params)){
    optparams = c(round(length(interval)/2),round(length(interval)/10),df_event$yy[1],df_event$yy[length(interval)]-df_event$yy[1])
  }else{
    optparams = control.linramp$opt.params
  }

  n=length(interval)
  df0=data.frame(y=df_event$yy,x=df_event$xx)

  t_start = df0$time[1];t_end = df0$time[n];y_start = df0$y[1]
  #library(numDeriv)
  #library(compiler)
  y = df0$y
  if(print.progress) cat("Using 'optim' to find initial positions for hyperparameters in INLA.\n",sep="")
  ## for stability the x-axis is transformed to values ranging from 1:n.

  timepoints = round((df_event$xx - df_event$xx[1])/(df_event$xx[n]-df_event$xx[1])*(n-1)+1,digits=4)
  df0$time = timepoints

  ## Finding initial values in INLA optimization procedure by first using a simple optimization
  ## gradient and cost function given here
  minfun.grad = function(param, args = NULL){
    return (grad(minfun, param, args=args, method.args = list(r=6)))
  }
  minfun = function(param, args = NULL){
    yhat = linramp(timepoints,t0=param[1],dt=param[2],y0=param[3],dy=param[4])
    mse = sum((args$y-yhat)^2)
    return(sqrt(mse))
  }

  ## perform optimization to find good starting values for INLA
  param=optparams
  args=list(y=y)
  fit = optim(param,
              fn = minfun,
              gr = minfun.grad,
              method = "BFGS",
              control = list(
                abstol = 0,
                maxit = 100000,
                reltol = 1e-11),
              args = args)



  ### use least squares estimates for fixed effects as initial values in inla

  muvek = linramp(timepoints,t0=fit$par[1],dt=fit$par[2],y0=fit$par[3],dy=fit$par[4])

  init = c(fit$par[1],log(fit$par[2]),fit$par[3],fit$par[4],0,0)

  if(print.progress) cat("Fitting linear ramp model in INLA using rgeneric model specification...\n",sep="")
  ## creating linear ramp INLA model using rgeneric framework. Requires further specification, see "rgeneric.uneven" function
  time.startinla = Sys.time()
  model.rgeneric = inla.rgeneric.define(rgeneric.uneven.AR1,n=n,tstart=timepoints[1],tslutt=timepoints[n],ystart=y[1],timepoints = timepoints)
  formula = y ~ -1+ f(idx, model=model.rgeneric)

  r = inla(formula,family="gaussian", data=data.frame(y=df0$y,idx=as.integer(df0$time)),
           control.mode=list(theta=init,
                             restart=TRUE),
           num.threads = control.linramp$ncores,
           verbose=control.linramp$verbose,
           control.inla=list(h=control.linramp$h),
           control.family = list(hyper = list(prec = list(initial = 12, fixed=TRUE))) )#, num.threads = 1)
  time.endinla = Sys.time()
  elapsedinla = difftime(time.endinla,time.startinla,units="secs")[[1]]
  if(print.progress) cat("Completed in ",elapsedinla," seconds.\n",sep="")
  if(print.progress) cat("Gathering results...\n",sep="")
  object$linramp = list(timepoints=timepoints,data=df0,inlafit=r)

  ## compute posterior marginals and posterior marginal means of z^* = t0 (transition onset), dt (transition duration), y0 (initial ramp level), dy (change in ramp level) sigma (amplitude of AR(1) noise) and tau (parameter for correlation of AR(1) noise)
  t0=inla.emarginal(function(x)x,r$marginals.hyperpar$`Theta1 for idx`); dt=inla.emarginal(function(x)exp(x),r$marginals.hyperpar$`Theta2 for idx`);y0=inla.emarginal(function(x)x,r$marginals.hyperpar$`Theta3 for idx`); dy=inla.emarginal(function(x)x,r$marginals.hyperpar$`Theta4 for idx`);rho = inla.emarginal(function(x)2/(1+exp(-x))-1,r$marginals.hyperpar$`Theta5 for idx`); sigma = inla.emarginal(function(x)1/sqrt(exp(x)),r$marginals.hyperpar$`Theta6 for idx`)
  muvek = linramp(timepoints,t0=t0,dt=dt,y0=y0,dy=dy)

  t0mean = r$summary.hyperpar$mean[1]; t0lower = r$summary.hyperpar$`0.025quant`[1];t0upper = r$summary.hyperpar$`0.975quant`[1]
  #margt0 = inla.tmarginal(function(x)df$age[1]+x/(n-1)*(df$age[n]-df$age[1]),r$marginals.hyperpar$`Theta1 for idx`);
  margt0 = inla.tmarginal(function(x)df0$x[1]+x/(n-1)*(df0$x[n]-df0$x[1]),r$marginals.hyperpar$`Theta1 for idx`);
  z.t0 = inla.zmarginal(margt0,silent=TRUE)

  object$linramp$param$t0 = list(marg.t0=margt0,mean=z.t0$mean,sd=z.t0$sd,q0.025=z.t0$quant0.025,q0.5=z.t0$quant0.5,q0.975=z.t0$quant0.975)
  if((abs(r$summary.hyperpar$mean[2])>1000) || (r$summary.hyperpar$sd[2]>1000)){
    margdt=NA
    margdtpos=NA
  }else{
    margdt = inla.tmarginal(function(x)exp(x)/(n-1)*(df0$x[n]-df0$x[1]),inla.smarginal(r$marginals.hyperpar$`Theta2 for idx`))
    margdtpos = data.frame(x=-margdt$x,y=margdt$y)
    z.dt = inla.zmarginal(margdt,silent=TRUE)
    z.dtpos = inla.zmarginal(margdtpos,silent=TRUE)
    object$linramp$param$dt = list(marg.dt=margdt,mean=z.dt$mean,sd=z.dt$sd,q0.025=z.dt$quant0.025,q0.5=z.dt$quant0.5,q0.975=z.dt$quant0.975)
    object$linramp$param$dtpos = list(marg.dtpos=margdtpos,mean=z.dtpos$mean,sd=z.dtpos$sd,q0.025=z.dtpos$quant0.025,q0.5=z.dtpos$quant0.5,q0.975=z.dtpos$quant0.975)
  }
  margy0 = inla.tmarginal(function(x)x,inla.smarginal(r$marginals.hyperpar$`Theta3 for idx`))
  margdy = inla.tmarginal(function(x)x,inla.smarginal(r$marginals.hyperpar$`Theta4 for idx`))
  z.y0 = inla.zmarginal(margy0,silent=TRUE)
  z.dy = inla.zmarginal(margdy,silent=TRUE)
  object$linramp$param$y0 = list(marg.y0=margy0,mean=z.y0$mean,sd=z.y0$sd,q0.025=z.y0$quant0.025,q0.5=z.y0$quant0.5,q0.975=z.y0$quant0.975)
  object$linramp$param$dy = list(marg.dy=margdy,mean=z.dy$mean,sd=z.dy$sd,q0.025=z.dy$quant0.025,q0.5=z.dy$quant0.5,q0.975=z.dy$quant0.975)

  margsigma = inla.tmarginal(function(x)1/sqrt(exp(x)),inla.smarginal(r$marginals.hyperpar$`Theta5 for idx`))
  margtau = inla.tmarginal(function(x)exp(x),inla.smarginal(r$marginals.hyperpar$`Theta6 for idx`))
  z.sigma = inla.zmarginal(margsigma,silent=TRUE)
  z.tau = inla.zmarginal(margtau,silent=TRUE)
  object$linramp$param$sigma = list(marg.sigma=margsigma,mean=z.sigma$mean,sd=z.sigma$sd,q0.025=z.sigma$quant0.025,q0.5=z.sigma$quant0.5,q0.975=z.sigma$quant0.975)
  object$linramp$param$tau = list(marg.sigma=margtau,mean=z.tau$mean,sd=z.tau$sd,q0.025=z.tau$quant0.025,q0.5=z.tau$quant0.5,q0.975=z.tau$quant0.975)


  if(control.linramp$nsims>0 || control.linramp$nsims>0) time.startbonussample = Sys.time()

  if(control.linramp$nsims>0 && print.progress && control.linramp$nsims==0) cat("Simulating ensemble of ", control.linramp$nsims, " samples for t1 = t0 + dt...","\n",sep="")
  if(control.linramp$nsims==0 && print.progress && control.linramp$nsims>0) cat("Simulating ensemble of ",control.linramp$nsims," samples for linear ramp...","\n",sep="")
  if(control.linramp$nsims>0 && print.progress && control.linramp$nsims>0) cat("Simulating ensembles of ",control.linramp$nsims," samples for t1 = t0 + dt and ",control.linramp$nsims," samples for linear ramp...","\n",sep="")

  if(control.linramp$nsims>0 || control.linramp$nsims>0){
    nsims = max(control.linramp$nsims,control.linramp$nsims)
    samps=inla.hyperpar.sample(nsims,object$linramp$inlafit)

    hpars = matrix(NA,nrow = nsims,ncol=5)
    hpars[,1:2]=samps[,3:4] #y0,dy
    n=length(timepoints)

    hpars[,3] = df0$x[1] + (samps[,1]-1)/(n-1)*(df0$x[n]-df0$x[1]) #t0
    hpars[,4] = exp(samps[,2])/(n-1)*(df0$x[n]-df0$x[1]) #Dt
  }

  if(control.linramp$nsims>0){
    t1sims=numeric(nsims)
    t1mean = mean(t1sims)
    for(i in 1:control.linramp$nsims){
      t01 = hpars[i,3]
      dt1 = hpars[i,4]
      t1sims[i] = t01+dt1

    }

    t1dens = density(t1sims)
    margt1 = cbind(t1dens$x,t1dens$y); colnames(margt1) = c("x","y")
    z.t1 = inla.zmarginal(margt1,silent=TRUE)
    object$linramp$param$t1 = list(marginal=margt1,samples = t1sims,mean=z.t1$mean,sd=z.t1$sd,q0.025=z.t1$quant0.025,q0.5=z.t1$quant0.5,q0.975=z.t1$quant0.975)
  }
  if(control.linramp$nsims>0){
    vekmat = matrix(NA,nrow=n,ncol=control.linramp$nsims)
    for(i in 1:control.linramp$nsims){
      t01 = hpars[i,3]
      dt1 = hpars[i,4]
      vekmat[,i] = bremla::linramprev(object$linramp$data$x,t0=t01,dt=dt1,y0=hpars[i,1],dy=hpars[i,2])
    }
    vek.quant0.025 = numeric(n)
    vek.quant0.5 = numeric(n)
    vek.quant0.975 = numeric(n)
    vek.mean = numeric(n)
    for(iter in 1:n){
      dens = density(vekmat[iter,])
      vek.quant0.025[iter]=INLA::inla.qmarginal(0.05,dens)
      vek.quant0.5[iter]=INLA::inla.qmarginal(0.5,dens)
      vek.mean[iter] = mean(vekmat[iter,])
      vek.quant0.975[iter]=INLA::inla.qmarginal(0.95,dens)
    }

    object$linramp$linrampfit = list(mean = vek.mean,q0.025=vek.quant0.025,q0.5=vek.quant0.5,q0.975=vek.quant0.975)
  }
  object$time$linramp = list(inla=elapsedinla)
  if(control.linramp$nsims>0 || control.linramp$nsims>0) object$time$t1_and_ramp = difftime(Sys.time(),time.startbonussample,units="secs")[[1]]
  if((control.linramp$nsims>0 || control.linramp$nsims>0) && print.progress) cat(" completed in ",object$time$t1_and_ramp," seconds!\n",sep="")

  time.total = difftime(Sys.time(),time.start,units="secs")[[1]]

  object$linramp$.args = list(nsims=control.linramp$nsims,
                              depth.reference=control.linramp$depth.ref,
                              label=control.linramp$depth.label)
  object$time$linramp=time.total

  return(object)

}
