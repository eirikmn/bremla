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
plot.bremla = function(x,
                       plot.proxydata=list(age=TRUE,depth=FALSE,xrev=FALSE,label=NULL),
                       plot.ls = list(fitted=TRUE,legend=NULL,residuals=TRUE,histogram=TRUE,qqplot=TRUE,acf=TRUE,xrev=FALSE,
                                      label.fit=NULL,label.res=NULL,label.hist=NULL,label.qq=NULL,label.acf=NULL),
                       plot.inla.posterior = list(posteriors=TRUE,label=NULL),
                       plot.inlasims = list(nsims=30,legend=NULL,xrev=FALSE,label=NULL),
                       plot.bias = list(MCE=NULL,legend=NULL,xrev=FALSE,label=NULL),
                       plot.linramp = list(depth.reference=NULL,show.t0=TRUE,show.t1=TRUE,xrev=TRUE,label=NULL),
                       plot.DO_depth = list(depth.reference=NULL,xrev=TRUE,label=NULL),
                       plot.DO_age = list(age.reference=NULL,xrev=TRUE,label=NULL),
                       postscript=FALSE,
                       pdf=FALSE,
                       prefix = "bemla.plots/figure-",
                       ...){
  if(!postscript && !pdf){
    dev=getOption("device")
    if(!is.character(dev) || dev != "RStudioGD"){
      dev.new(...)
    }
  }else{
    dir = dirname(prefix)
    if (!file.exists(dir) && nchar(dir) > 0L) {
      dir.create(dir, recursive=TRUE)
    } else {
      stopifnot(file.info(dir)$isdir)
    }
  }

  figure.count = 1L
  gicc05 = x$data$y; z = x$data$z; xx = x$data$x; n=length(gicc05)
  eventindexes = x$.args$eventindexes; nevents = length(eventindexes)
  oldpar = par()

  if(!is.null(plot.proxydata)){
    figure.count <- new.plot(postscript,pdf,prefix,figure.count,...) +1
    par(mfrow=c(1,1),mar=(c(5,4.5,4,2)+0.1))
    if(plot.proxydata$age){
      xdata = gicc05; xlab="GICC05 (yr b2k)"
    }else if(plot.proxydata$depth){
      xdata = z; xlab = "Depth (m)"
    }
    xlim = c(xdata[1],xdata[n])
    if(plot.proxydata$xrev) xlim=rev(xlim)
    plot(xdata,xx,type="l",xlab=xlab,ylab=expression(paste(delta^18,"O (permil)")),main=plot.proxydata$label,xlim=xlim)
    if(nevents>0) abline(v=xdata[eventindexes],lwd=0.6,col="gray")#rgb(red=0.5,green=0.5,blue=0.5,alpha=1),lwd=0.8)

    if(postscript || pdf){
      if (names(dev.cur()) != "null device") {
        dev.off()
      }
    }
  }
  dy=x$data$dy

  if(!is.null(plot.ls) && !is.null(x$LS.fitting$fit)){
    figure.count <- new.plot(postscript,pdf,prefix,figure.count,...) +1
    par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)
    if(plot.ls$fitted){
      xlim = range(z)
      if(plot.ls$xrev) xlim=rev(xlim)
      if(is.null(plot.ls$label.fit)){
        plot.label="Least squares fit"
      }else{
        plot.label=plot.ls$label.fit
      }
      plot(z,dy,type="l",xlab=paste0("Depth (m)"),xlim=xlim,ylab=("Layers per 5 cm"),main=plot.label)
      lines(z,x$LS.fitting$fit$fitted.values,col="red")
      abline(v=z[eventindexes],col="gray",lwd=0.8)
    }
    if(!is.null(plot.ls$legend)) legend(x=leg$x,y=leg$y,legend=leg$legend,col=leg$col,lty=leg$lty,cex=leg$cex,pch=leg$pch,lwd=leg$lwd,pt.cex=leg$pt.cex,bty=leg$bty)

    if(plot.ls$residuals){
      if(is.null(plot.ls$label.res)){
        plot.label="Least squares residuals"
      }else{
        plot.label=plot.ls$label.res
      }
      plot(z,x$LS.fitting$fit$residuals,type="l",xlab="Depth (m)",ylab="Residual errors per 5 cm",main=plot.label,xlim=xlim); abline(h=0,lty=3,col="gray")
    }
    if(plot.ls$histogram){
      if(is.null(plot.ls$label.hist)){
        plot.label="Histogram"
      }else{
        plot.label=plot.ls$label.hist
      }
      hist(x$LS.fitting$fit$residuals,freq=0,col="orange",breaks=20,xlab="Residual errors per 5 cm", main=plot.label)
    }
    if(plot.ls$qqplot){
      if(is.null(plot.ls$label.qq)){
        plot.label="Q-Q Plot"
      }else{
        plot.label=plot.ls$label.qq
      }
      qqnorm(x$LS.fitting$fit$residuals,main=plot.label); qqline(x$LS.fitting$fit$residuals)
    }
    if(plot.ls$acf){
      if(is.null(plot.ls$label.acf)){
        plot.label="Autocorrelation function"
      }else{
        plot.label=plot.ls$label.acf
      }
      acf(x$LS.fitting$fit$residuals,lag.max = 30,main=plot.label,lwd=2)
    }
    if(postscript || pdf){
      if (names(dev.cur()) != "null device") {
        dev.off()
      }
    }
  }

  if(!is.null(plot.inla.posterior) && !is.null(x$fitting$fit)){
    figure.count <- new.plot(postscript,pdf,prefix,figure.count,...) +1
    if(plot.inla.posterior$posteriors){
      if(length(plot.inla.posterior$label)<=1){
        plot.label1 = plot.inla.posterior$label
        plot.label2 = plot.inla.posterior$label
        plot.label3 = plot.inla.posterior$label
      }else{
        plot.label1 = plot.inla.posterior$label[1]
      }
      plot(x$fitting$hyperparameters$posteriors$sigma_epsilon,type="l",xlab=expression(paste(sigma[epsilon])),ylab="Density",lwd=2,main=plot.label1)
      abline(v=x$fitting$hyperparameters$results$sigma_epsilon$mean)
      abline(v=c(x$fitting$hyperparameters$results$sigma_epsilon$quant0.025,x$fitting$hyperparameters$results$sigma_epsilon$quant0.975),col="gray")

      if(tolower(x$.args$noise) %in% c(1,"ar1","ar(1)")){
        if(length(plot.inla.posterior$label)>=2){
          plot.label2 = plot.inla.posterior$label[2]
        }
        plot(x$fitting$hyperparameters$posteriors$phi,xlab=expression(paste(phi)),ylab="Density",lwd=2,type="l",main=plot.label2)
        abline(v=x$fitting$hyperparameters$results$phi$mean)
        abline(v=c(x$fitting$hyperparameters$results$phi$quant0.025,x$fitting$hyperparameters$results$phi$quant0.975),col="gray")
      }else if(tolower(x$.args$noise) %in% c(2,"ar2","ar(2)")){
        if(length(plot.inla.posterior$label)>=3){
          plot.label2 = plot.inla.posterior$label[2]
          plot.label3 = plot.inla.posterior$label[3]
        }
        plot(x$fitting$hyperparameters$posteriors$phi1,xlab=expression(paste(phi[1])),ylab="Density",lwd=2,type="l",main=plot.label2)
        abline(v=x$fitting$hyperparameters$results$phi1$mean)
        abline(v=c(x$fitting$hyperparameters$results$phi1$quant0.025,x$fitting$hyperparameters$results$phi1$quant0.975),col="gray")

        plot(x$fitting$hyperparameters$posteriors$phi2,xlab=expression(paste(phi[2])),ylab="Density",lwd=2,type="l",main=plot.label3)
        abline(v=x$fitting$hyperparameters$results$phi2$mean)
        abline(v=c(x$fitting$hyperparameters$results$phi2$quant0.025,x$fitting$hyperparameters$results$phi2$quant0.975),col="gray")
      }
    }
    if(postscript || pdf){
      if (names(dev.cur()) != "null device") {
        dev.off()
      }
    }
  }

  if(!is.null(plot.inlasims) && !is.null(x$simulation)){
    figure.count <- new.plot(postscript,pdf,prefix,figure.count,...) +1
    nsims = plot.inlasims$nsims
    xlim=range(z)
    if(plot.inlasims$xrev) xlim=rev(xlim)
    ylim=range(x$simulation$summary$hpd0.025-gicc05,x$simulation$summary$hpd0.975-gicc05)*1.1
    plot(NA,xlim=xlim,ylim=ylim,xlab="Depth (m)",ylab="Simulated time scale - GICC05 (years)",main=plot.inlasims$label)
    for(iter in 1:nsims){
      lines(z,x$simulation$age[,iter]-gicc05,col="gray",lwd=0.8)
    }
    lines(z,x$simulation$summary$hpd0.025-gicc05,col="red",lwd=2)
    lines(z,x$simulation$summary$hpd0.975-gicc05,col="red",lwd=2)
    abline(h=0,lty=3)
    lines(z,x$simulation$summary$mean-gicc05,col="blue",lwd=2)
    if(postscript || pdf){
      if (names(dev.cur()) != "null device") {
        dev.off()
      }
    }
  }

  if(!is.null(plot.bias) && !is.null(x$biases)){
    figure.count <- new.plot(postscript,pdf,prefix,figure.count,...) +1
    nbiases = length(x$biases)

    yrange = c(0,0)
    for( iter in 1:nbiases){
      yrange = range(yrange,x$biases[[paste0("bias",iter)]]$quant0.975-gicc05,x$biases[[paste0("bias",iter)]]$quant0.025-gicc05)
    }
    if(!is.null(plot.bias$MCE)){
      yrange = range(yrange,-plot.bias$MCE,plot.bias$MCE)
    }
    xlim=range(z)
    if(!is.null(plot.bias$legend)){
      yrange = range(yrange,plot.bias$legend$y)
      xlim=range(xlim,plot.bias$legend$y)
    }

    if(plot.bias$xrev) xlim=rev(xlim)
    plot(NA, xlim=xlim,xlab="Depth (m)", ylab="Simulated timescale - GICC05 (years)",type="l",col="blue",ylim=yrange,main=plot.bias$label)
    abline(h=0,lty=3,lwd=1,col="gray")
    for(iter in 1:nbiases){
      lines(z,x$biases[[paste0("bias",iter)]]$quant0.975-gicc05,col="blue",lty=iter)
      lines(z,x$biases[[paste0("bias",iter)]]$quant0.025-gicc05,col="blue",lty=iter)
    }
    if(!is.null(plot.bias$MCE)){
      lines(z,plot.bias$MCE,lwd=2)
      lines(z,-plot.bias$MCE,lwd=2)
    }
    leg=plot.bias$legend
    legend(leg$x,leg$y,leg$legend,col=leg$col,lty=leg$lty,cex=leg$cex,pch=leg$pch,lwd=leg$lwd,pt.cex=leg$pt.cex)
    if(postscript || pdf){
      if (names(dev.cur()) != "null device") {
        dev.off()
      }
    }
  }

  if(!is.null(plot.linramp) && !is.null(x$linramp)){
    figure.count <- new.plot(postscript,pdf,prefix,figure.count,...) +1
    par(mfrow=c(1,1),mar=(c(5,4.5,4,2)+0.1))
    yrange = range(x$linramp$data$y,x$linramp$linrampfit$q0.025,x$linramp$linrampfit$q0.975)
    xlim = range(x$linramp$data$x)
    if(plot.linramp$xrev) xlim=rev(xlim)
    xval = x$linramp$data$x
    if(!is.null(plot.linramp$label)){
      plot.label=plot.linramp$label
    }else{
      plot.label=x$linramp$.args$label
    }
    plot(xval,x$linramp$data$y,type="l",lwd=1.25,col="gray",xlim=xlim,ylim=yrange,xlab="Depth (m)",ylab=expression(paste(delta^18, "O (permil)")),main=plot.label)
    lines(xval,x$linramp$linrampfit$q0.025,col="red",lwd=2)
    lines(xval,x$linramp$linrampfit$q0.975,col="red",lwd=2)
    lines(xval,x$linramp$linrampfit$mean,col="black",lwd=2)

    if(plot.linramp$show.t0){
      ybottom = min(x$linramp$data$y)
      ytop = min(x$linramp$linrampfit$q0.025)
      margt0=x$linramp$param$t0$marg.t0
      normt0.y = margt0[,2]/diff(range(margt0[,2]))*(ytop-ybottom)-min(margt0[,2])+ybottom
      lines(x=margt0[,1],y=normt0.y,col="blue",lwd=2)
    }
    if(plot.linramp$show.t0 && !is.null(x$linramp$param$t1)){
      ybottom = min(x$linramp$data$y)
      ytop = min(x$linramp$linrampfit$q0.025)
      margt1 = x$linramp$param$t1$marginal
      normt1.y = margt1[,2]/diff(range(margt1[,2]))*(ytop-ybottom)-min(margt1[,2])+ybottom
      lines(x=margt1[,1],y=normt1.y,col="blue",lwd=2,lty=3)
    }

    if(!is.null(plot.DO_depth$depth.reference)){
      abline(v=plot.DO_age$age.reference,lwd=2,lty=3)
    }else{
      if(!is.null(x$linramp$.args$depth.reference)){
        abline(v= x$linramp$.args$depth.reference,lwd=2,lty=3)
      }
    }
  }

  if(!is.null(plot.DO_depth) && !is.null(plot.linramp)){
    figure.count <- new.plot(postscript,pdf,prefix,figure.count,...) +1
    par(mfrow=c(1,1),mar=(c(5,4,4,2)+0.1))
    if(is.null(plot.DO_depth$label)){
      plot.label = x$linramp$.args$label
    }else{
      plot.label = plot.DO_depth$label
    }
    xlim = range(x$linramp$param$t0$marg.t0[,1])
    if(!is.null(plot.DO_depth$depth.reference)){
      xlim=range(xlim,plot.DO_depth$depth.reference)
    }else{
      if(!is.null(x$linramp$.args$depth.reference)){
        xlim=range(xlim,x$linramp$.args$depth.reference)
      }
    }
    if(plot.DO_age$xrev) xlim=rev(xlim)
    plot(x$linramp$param$t0$marg.t0,type="l",lwd=2,xlab="Onset depth (m)",ylab="Density",main=plot.label,xlim=xlim)
    abline(v=x$linramp$param$t0$mean,lwd=2)
    abline(v=c(x$linramp$param$t0$q0.025,x$linramp$param$t0$q0.975),lwd=2,col="gray")

    if(!is.null(plot.DO_depth$depth.reference)){
      abline(v=plot.DO_age$age.reference,lwd=2,lty=3)
    }else{
      if(!is.null(x$linramp$.args$depth.reference)){
        abline(v= x$linramp$.args$depth.reference,lwd=2,lty=3)
      }
    }
  }

  if(!is.null(plot.DO_age) && !is.null(x$DO_dating)){
    xlim = range(x$DO_dating$samples)
    if(plot.DO_age$xrev) xlim=rev(xlim)
    hist(x$DO_dating$samples,xlab="Onset age (y b2k)",freq=FALSE,col="orange",breaks=30,main=x$DO_dating$.args$label,xlim=xlim)
    abline(v=x$DO_dating$mean,lwd=2)
    abline(v=c(x$DO_dating$q0.025,x$DO_dating$q0.975),col="gray",lwd=2)
    if(!is.null(plot.DO_age$age.reference)){
      abline(v=plot.DO_age$age.reference,lwd=2,lty=3)
    }else{
      if(!is.null(x$DO_dating$.args$age.reference)){
        abline(v=x$DO_dating$.args$age.reference,lwd=2,lty=3)
      }
    }
  }

  return(invisible(x))
}


new.plot = function(postscript,pdf,prefix,figure.count,...)
{

  #dev = getOption("device")
  if(postscript && pdf){
    stop("Multiple file types have been seleced.")
  }
  if(postscript) {
    ending=".eps"
  }else if(pdf){
    ending=".pdf"
  }
  if(!postscript && !pdf){
    dev=getOption("device")
    if(!is.character(dev) || dev != "RStudioGD"){
      dev.new(...)
    }
  }else{

    file.found=FALSE
    while(!file.found){
      filename=paste(prefix,figure.count,ending,sep="")

      if(file.exists(filename)){
        figure.count <- figure.count +1L
      }else{
        file.found=TRUE
      }
    }
    if(postscript){
      postscript(file=filename,...)
    }else if(pdf){
      pdf(file=filename,...)
    }
  }
  return (invisible(figure.count))
}
