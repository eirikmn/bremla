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



  if(!is.null(object$fitting)){
    cpu = as.numeric(round(object$time$fit$total,digits=digits))
    cpu.navn="Model fitting"
  }

  if(!is.null(object$simulation)){
    cpu=as.numeric(c(cpu,round(object$time$simulation$total,digits=digits)))
    cpu.navn=c(cpu.navn,"Chronology sampling")
  }
  if(!is.null(object$linramp) && !is.null(object$DO_dating)){
    cpu=as.numeric(c(cpu,round(object$time$t1_and_ramp + object$time$DO_age$total,digits=digits)))
    cpu.navn=c(cpu.navn,"DO dating")
  }
  if(!is.null(object$biases)){
    cpu=as.numeric(c(cpu,round(object$time$biases,digits=digits)))
    cpu.navn=c(cpu.navn,"Bias sampling")
  }

  cpu=as.numeric(c(cpu,round(object$time$total,digits=digits)))
  cpu.navn=c(cpu.navn,"Total")
  names(cpu)=cpu.navn
  ut=c(ut, list(cpu.used=cpu))

  if(tolower(object$.args$noise) %in% c(0,"iid","independent")){
    noise = "iid"
  }else if(tolower(object$.args$noise) %in% c(1,"ar1","ar(1)")){
    noise = "ar1"
  }else if(tolower(object$.args$noise) %in% c(2,"ar2","ar(2)")){
    noise = "ar2"
  }

  if(!is.null(object$fitting)){
    if(noise == "iid"){

      hypers = matrix(round(as.numeric(object$fitting$hyperparameters$results$sigma_epsilon),digits=digits),nrow=1)
      #hypers = matrix( round(as.numeric(object$fitting$hyperparameters$results$sigma_epsilon),digits=digits),ncol=mm )
      hypers = as.data.frame(hypers)
      colnames(hypers) = c("mean","sd","quant0.025","quant0.25","quant0.5","quant0.75","quant0.975")
      rownames(hypers) = c("sigma_epsilon")
    }else if(noise == "ar1"){
      hypers = rbind(round(as.numeric(object$fitting$hyperparameters$results$sigma_epsilon),digits=digits),
                     round(as.numeric(object$fitting$hyperparameters$results$phi),digits=digits))
      hypers = as.data.frame(hypers)
      colnames(hypers) = c("mean","sd","quant0.025","quant0.25","quant0.5","quant0.75","quant0.975")
      rownames(hypers) = c("sigma_epsilon","phi")
    }else if(noise == "ar2"){
      hypers = rbind(round(as.numeric(object$fitting$hyperparameters$results$sigma_epsilon),digits=digits),
                     round(as.numeric(object$fitting$hyperparameters$results$phi1),digits=digits),
                     round(as.numeric(object$fitting$hyperparameters$results$phi2),digits=digits))
      hypers = as.data.frame(hypers)
      colnames(hypers) = c("mean","sd","quant0.025","quant0.25","quant0.5","quant0.75","quant0.975")
      rownames(hypers) = c("sigma_epsilon","phi1","phi2")
    }
    maxlength=2048L
    formulastring = object$.args$ls.formulastring
    if(sum(nchar(formulastring)) > maxlength){
      formulastring = paste0( substr(deparse(object$.args$formulastring),1L,maxlength),"...")
    }
    fit.arg = list(formula = formulastring,noise=object$.args$noise, nevents=object$.args$nevents, method = object$.args$method)
  }

  ut=c(ut, list(hyperpar=hypers,fit.arg=fit.arg))

  if(!is.null(object$simulation)){
    sim = list(nsims = dim(object$simulation$age)[2],n = dim(object$simulation$age)[1], store.means = !is.null(object$simulation$dmean) )
  }
  ut = c(ut,sim)

  if(!is.null(object$linramp)){
    hyperramp = rbind(c(object$linramp$param$t0$mean,object$linramp$param$t0$sd,object$linramp$param$t0$q0.025,object$linramp$param$t0$q0.5,object$linramp$param$t0$q0.975),
                      c(object$linramp$param$dtpos$mean,object$linramp$param$dtpos$sd,object$linramp$param$dtpos$q0.025,object$linramp$param$dtpos$q0.5,object$linramp$param$dtpos$q0.975),
                      c(object$linramp$param$y0$mean,object$linramp$param$y0$sd,object$linramp$param$y0$q0.025,object$linramp$param$y0$q0.5,object$linramp$param$y0$q0.975),
                      c(object$linramp$param$dy$mean,object$linramp$param$dy$sd,object$linramp$param$dy$q0.025,object$linramp$param$dy$q0.5,object$linramp$param$dy$q0.975)
                      )
    colnames(hyperramp) = c("mean","sd","quant0.025","quant0.5","quant0.975")
    if(object$linramp$.args$t1.sims>0){
      hyperramp = rbind(hyperramp, c(object$linramp$param$t1$mean,object$linramp$param$t1$sd,object$linramp$param$t1$q0.025,object$linramp$param$t1$q0.5,object$linramp$param$t1$q0.975))
      rownames(hyperramp) = c("t0", "dt", "y0","dy","t1")
    }else{
      rownames(hyperramp) = c("t0", "dt", "y0","dy")
    }
    hyperramp = as.data.frame(round(hyperramp,digits=digits))
    linramplist = list(hyperramp = hyperramp,t1sims=object$linramp$.args$t1.sims,rampsims = object$linramp$.args$rampsims,label=object$linramp$.args$label,depth.reference=object$linramp$.args$depth.reference)
    ut = c(ut,linramplist)
  }

  if(!is.null(object$DO_dating)){
    DO_age = matrix(c(object$DO_dating$mean,object$DO_dating$sd,object$DO_dating$q0.025,object$DO_dating$q0.5,object$DO_dating$q0.975),nrow=1)
    colnames(DO_age) = c("mean","sd","quant0.025","quant0.5","quant0.975")
    rownames(DO_age) = "Onset age"
    DO_age = as.data.frame(round(DO_age,digits=digits))
    DOlist = list(DO_age=DO_age,datingsims = object$DO_dating$.args$nsims,label=object$DO_dating$.args$label,age.reference=object$DO_dating$.args$age.reference)
    ut = c(ut,DOlist)
  }

  if(!is.null(object$biases)){
    nbiases = object$biases$.args$nbiases

    if(nbiases == 1){
      biasparam = object$biases$.args$biasparam
    }else{
      biasparam = object$biases$.args$biasparam
    }

    biaslist = list(nbiases = nbiases, bias.model = object$biases$.args$bias.model, biasparam = biasparam, store.samples = object$biases$.args$store.samples,nsims=object$biases$.args$nsims)
    ut = c(ut,biaslist)
  }


  class(ut) = "summary.bremla"

return(ut)
}

print.summary.bremla = function(x,digits=4L,...){
  cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
  cat("Time used:\n")
  print(x$cpu.used)
  cat("\n",sep="")

  cat("The fixed component is explained by linear predictor: \n",x$fit.arg$formula,"\n\nThe noise component is explained by an ",x$fit.arg$noise," process.\n",sep="")

  if(tolower(x$fit.arg$method) %in% c("inla")){
    cat("\nThe model is fitted using INLA, with following estimates for the hyperparameters:\n")
    print(format(x$hyperpar,digits=digits,nsmall=2),quote=FALSE)
  }else{
    cat("\nThe model is fitted using least squares, with following estimates for the model parameters:\n")
    print(format(x$hyperpar,digits=digits,nsmall=2),quote=FALSE)
  }

  if(!is.null(x$nsims)){
    cat("\nSimulating ",x$nsims," chronologies.\n",sep="")
    if(x$store.means){
      cat(" Storing means")
    }
    cat("\n")
  }

  if(!is.null(x$hyperramp)){
    if(is.null(x$label)){
      lab="unlabeled"
    }else{
      lab=x$label
    }

    if(!is.null(x$depth.reference)){
      cat("\nOnset detection for ",x$label," event (reference onset depth is ",x$depth.reference,"):\n",sep="")
      #cat("The reference onset depth is ",x$depth.reference,"\n",sep="")
    }else{
      cat("\nOnset detection for ",x$label," event:\n",sep="")
    }

    print(x$hyperramp)

    if(x$rampsims>0) cat("\n",x$rampsims, " samples of linear ramp function produced.\n",sep="")

    if(!is.null(x$age.reference)){
      if(!is.null(x$DO_age)) cat("\nGenerated ",x$datingsims, " samples of onset ages (reference onset age is ",x$age.reference,").\n",sep="")
    }else{
      if(!is.null(x$DO_age)) cat("\nGenerated ",x$datingsims, " samples of onset ages.\n",sep="")
    }

    if(!is.null(x$DO_age)) print(x$DO_age)
    #if(!is.null(x$age.reference)) cat("The reference onset age is ",x$age.reference,"\n",sep="")
  }

  if(!is.null(x$biasparam)){
    cat("\nGenerated ",x$nbiases, " sets of (",x$bias.model,") biased chronologies, with parameters:\n",sep="")
    print(x$biasparam)
  }



}
