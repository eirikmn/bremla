#' Plot bremla model
#'
#' Plots results from bremla S3 class.
#'
#' @param x \code{bremla} S3 class. Output of \code{\link{bremla}} function
#' @param plot.proxydata List specifying how the proxies should be plotted. Either as a function of age or depth. Also includes title (\code{label}) and a boolean describing whether the x-axis should be reversed to give a chronological representation (\code{xrev})
#' @param plot.ls List specifying how the least square fit should be illustrated. \code{fitted=TRUE} plots the fitted values with label \code{label.fit}, \code{legend} specifies the legend, \code{residuals=TRUE} means the residuals are plotted with title \code{label.res}.
#' if \code{histogram=TRUE} the residuals are represented in a histogram with label \code{label.hist}, if \code{qqplot=TRUE} a quantile-quantile plot of the residuals are plotted with title \code{label.qq} and if
#' \code{acf=TRUE} the empirical autocorrelations are plotted with title \code{label.acf}.
#' @param plot.inla.posterior list specifying how the results from the inla regression fit should be plotted. If \code{posteriors=TRUE} then the posterior marginal distributions of the hyperparameters are plotted with title \code{label}.
#' @param plot.inlasims list specifying how the simulated chronologies from the INLA posterior should be plotted. \code{nsims} gives how many simulated chronologies should be included in the plot (with title \code{label}), \code{legend} specifies the legend, if \code{xrev=TRUE} the x-axis is reversed to give a chronological ordering.
#' @param plot.bias list specifying how the simulations under the assumptions of unknown counting bias should be represented. If \code{MCE} is given as a numeric vector it will be included in the plot (with title \code{label}). \code{legend} specifies the legend, if \code{xrev=TRUE} the x-axis is reversed to give a chronological ordering.
#' @param plot.linramp list specifying how the linear ramp fit should be plotted. If \code{depth.reference} is given, it will be represented by a vertical dotted line. If \code{show.t0=TRUE} the posterior marginal distribution of the onset is included (non-normalized). If \code{show.t1=TRUE} the posterior marginal distribution of the end point of the transition is included (non-normalized). If \code{xrev=TRUE} the x-axis is reversed to give a chronological ordering. \code{label} gives the title of the plot.
#' @param plot.DO_depth list specifying how the posterior distribution of the onset depth should be plotted. If \code{depth.reference} is given, it will be represented by a vertical dotted line. If \code{xrev=TRUE} the x-axis is reversed to give a chronological ordering. \code{label} gives the title of the plot.
#' @param plot.DO_age list specifying how the histogram of the simulated onset ages should be plotted. If \code{age.reference} is given, it will be represented by a vertical dotted line. If \code{xrev=TRUE} the x-axis is reversed to give a chronological ordering. \code{label} gives the title of the plot.
#' @param postscript Boolean variable indicating if postscript files should be produced instead.
#' @param pdf Boolean variable indicating if pdf files should be produced instead.
#' @param prefix The prefix for created files. Additional numbering is added.
#' @param ... Additional arguments to \code{postscripts()}, \code{pdf()} or \code{dev.new()}.
#'
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{bremla},\link{bremla_chronology_simulation}}
#' @keywords bremla plot
#'
#' @examples
#'
#' @export
#' @importFrom INLA inla inla.tmarginal inla.zmarginal inla.ar.pacf2phi
#' @importFrom grDevices dev.cur dev.new dev.off
#' @importFrom graphics abline hist legend lines par
#' @importFrom stats qqline qqnorm
#'
summary.bremla = function(object,digits=4L,...){

ut=list()

maxlength=2048L
if(sum(nchar(object$.args$call)) > maxlength){
  ut=c(ut, list(call=paste0( substr(deparse(object$.args$call),1L,maxlength),"...") ) )
}else{
  ut=c(ut, list(call=object$.args$call))
}

cpu = round(object$time$inla,digits=digits)
cpu.navn="Running INLA"


if(!is.null(object$TCR)){
  cpu=c(cpu,round(object$time$TCR,digits=digits))
  cpu.navn=c(cpu.navn,"Sampling TCR")
}
if(!is.null(object$mu)){
  cpu=c(cpu,round(object$time$mu,digits=digits))
  cpu.navn=c(cpu.navn,"Sampling mu")
}


cpu=c(cpu,round(object$time$Total,digits=digits))
cpu.navn=c(cpu.navn,"Total")
names(cpu)=cpu.navn
ut=c(ut, list(cpu.used=cpu))

if(is.lrd){
  hypers=matrix(round(as.numeric(c(object$hyperparam$means,object$hyperparam$sd,object$hyperparam$quant0.025,
                                   object$hyperparam$quant0.5,object$hyperparam$quant0.975)),digits=digits),
                nrow=length(object$inla.result$summary.hyperpar$mean),ncol=5)
}else if(object$misc$model=="ar1"){
  hypers=matrix(round(as.numeric(c(object$hyperparam$means,object$hyperparam$sd,object$hyperparam$quant0.025,
                                   object$hyperparam$quant0.5,object$hyperparam$quant0.975)),digits=digits),
                nrow=(length(object$inla.result$summary.hyperpar$mean)+(object$misc$m>1)),ncol=5)
}

hypers=as.data.frame(hypers)
colnames(hypers)=c("mean","sd","0.025quant","0.5quant","0.975quant")
hyper.navn=c()
if(!object$misc$INLA.options$control.family$hyper$prec$fixed){
  hyper.navn=c(hyper.navn,"Precision for the Gaussian observations")
}
hyper.navn=c(hyper.navn,var.name,"Sigmax","Sigmaf","F0")
if(object$misc$model == "ar1"){
  if(object$misc$m == 1){
    hyper.navn = c(hyper.navn, "p")
  }else{
    for(k in 1:object$misc$m){
      hyper.navn = c(hyper.navn, paste0("w",k))
    }
    for(k in 1:object$misc$m){
      hyper.navn = c(hyper.navn, paste0("p",k))
    }
  }
}
rownames(hypers)=hyper.navn

ut=c(ut, list(hyperpar=hypers))

if(!is.null(object$TCR)){
  tcr=matrix(round(c(object$TCR$mean,object$TCR$sd,object$TCR$quant0.025,object$TCR$quant0.5,object$TCR$quant0.975),digits=digits),nrow=1)
  tcr=as.data.frame(tcr)
  colnames(tcr)=c("mean","sd","0.025quant","0.5quant","0.975quant")
  rownames(tcr)="TCR"
  ut=c(ut, list(TCR=tcr,tcr.nsamples=object$misc$TCR.options$nsamples,Qco2=object$misc$TCR.options$Qco2))
}

if(!is.null(object$mu)){
  mu=matrix(round(c(object$mu$mean,object$mu$sd,object$mu$quant0.025,object$mu$quant0.5,object$mu$quant0.975),digits=digits),nrow=1)
  mu=as.data.frame(mu)

  #### FIKS DETTE
  ut=c(ut, list(mu = mu,mu.full.Bayesian= object$misc$mu.options$compute.mu %in% c(2,"full","complete"),mu.nsamples=object$misc$mu.options$nsamples))
}

neffp = object$inla.result$neffp
ut = c(ut, list(neffp = round(neffp, digits)))

if(!is.null(object$inla.result$dic)){
  ut=c(ut,list(dic=object$inla.result$dic))
}
if(!is.null(object$inla.result$waic)){
  ut=c(ut,list(waic=object$inla.result$waic))
}
if(!is.null(object$inla.result$mlik)){
  ut=c(ut,list(mlik=object$inla.result$mlik))
}
if(!is.null(object$inla.result$cpo$cpo) && length(object$inla.result$cpo$cpo)>0L){
  ut=c(ut, list(cpo=lapply(object$inla.result$cpo,round,digits=digits)))
}
if(!is.null(object$inla.result$summary.linear.predictor)){
  ut=c(ut, list(linear.predictor=round(object$inla.result$summary.linear.predictor),digits=digits))
}
if(!is.null(object$inla.result$summary.random) && length(object$inla.result$summary.random)>0L){
  ut=c(ut, list(random.names=names(object$inla.result$summary.random),random.model=object$inla.result$model.random))
}
if(!is.null(object$inla.result$summary.linear.predictor)){
  ut=c(ut, list(linear.predictor=round(object$inla.result$summary.linear.predictor),digits=digits))
}
ut = c(ut, list(model = object$misc$model))
ut=c(ut, list(family=object$inla.result$family))
class(ut) = "summary.inla.climate"

return(ut)
}

print.summary.bremla = function(x,digits=4L,...){
  cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
  cat("Time used:\n")
  print(x$cpu.used)
  cat("\n",sep="")

  if(!is.null(x$random.names)){
    cat("Random effects:\n",sep="")
    cat("Name\t ","Model\n",sep="")
    for(i in 1:length(x$random.names)){
      cat(paste0(x$random.names[i]," ",x$random.model[i],"\n"))
    }
    cat("\n")
  }else{
    cat("The model has no random effects\n\n")
  }

  if(!is.null(x$hyperpar)){
    cat("Model hyperparameters:\n")
    print(format(x$hyperpar,digits=digits,nsmall=2),quote=FALSE)
    cat("\n",sep="")
  }else{
    cat("This model has no hyperparameters\n\n")
  }

  if(!is.null(x$TCR)){
    cat("Transient climate response computed from ",format(x$tcr.nsamples,digits=digits,scientific=FALSE)," samples using CO2 coefficient ",format(x$Qco2,digits=digits,scientific=FALSE),":\n",sep="")
    print(format(x$TCR,digits=digits))
    cat("\n")
  }

  if(!is.null(x$mu)){
    if(x$mu.full.Bayesian){
      cat("Full Bayesian analysis of forcing response computed from ",format(x$mu.nsamples,digits=digits,scientific=FALSE)," samples.\n\n",sep="")
    }else{
      cat("Quick computation of forcing response computed from ",format(x$mu.nsamples,digits=digits,scientific=FALSE)," samples.\n\n",sep="")
    }
  }



  if(!is.null(x$neffp)){
    cat("Expected number of effective parameters (std dev): ", format(x$neffp[1], digits=digits, nsmall=2)," (",
        format(x$neffp[2], digits=digits, nsmall=2),")\n", sep="")
    cat("Number of equivalent replicates: ", format(x$neffp[3], digits=digits, nsmall=2),"\n\n",sep="")
  } else {
    cat("Expected number of effective parameters and Number of equivalent replicates not computed\n\n")
  }

  if(!is.null(x$dic)){
    cat(paste0("Deviance Information Criterion (DIC) ...: ",
               format(x$dic$dic, digits=digits, nsmall=2), "\n",
               "Effective number of parameters .........: ",
               format(x$dic$p.eff, digits=digits, nsmall=2), "\n\n"))
  }

  if (!is.null(x$waic)){
    cat(paste0("Watanabe-Akaike information criterion (WAIC) ...: ",
               format(x$waic$waic, digits=digits, nsmall=2), "\n",
               "Effective number of parameters .................: ",
               format(x$waic$p.eff, digits=digits, nsmall=2), "\n\n"))
  }

  if(!is.null(x$mlik)){
    cat(paste("Marginal log-Likelihood: ", format(x$mlik[2], digits=digits, nsmall=2),"\n",sep=""))
  }
  if(!is.null(x$cpo)){
    cat("CPO and PIT are computed\n\n")
  }
  if(!is.null(x$linear.predictor)){
    cat("Posterior marginals for linear predictor and fitted values computed\n\n")
  }

}
