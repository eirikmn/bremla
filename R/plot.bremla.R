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
#' @param plot.syncsims list specifying how the simulated synchronized chronologies from the INLA posterior should be plotted. \code{nsims} gives how many simulated chronologies should be included in the plot (with title \code{label}), \code{legend} specifies the legend, if \code{xrev=TRUE} the x-axis is reversed to give a chronological ordering.
#' @param plot.tiepoints list containing specifications for plotting histogram of sampled tie-points.
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
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon theme_bw xlab ylab geom_segment geom_point
plot.bremla = function(x,
                       plot.proxydata=list(age=TRUE,depth=FALSE,xrev=FALSE,label=NULL),
                       plot.ls = list(fitted=TRUE,legend=NULL,residuals=TRUE,histogram=TRUE,qqplot=TRUE,acf=TRUE,xrev=FALSE,
                                      label.fit=NULL,label.res=NULL,label.hist=NULL,label.qq=NULL,label.acf=NULL),
                       plot.inla.posterior = list(posteriors=TRUE,label=NULL),
                       plot.inlasims = list(nsims=0,legend=NULL,xrev=FALSE,label=NULL),
                       plot.syncsims = list(nsims=0,legend=NULL,xrev=FALSE,label=NULL),
                       plot.tiepoints = list(col.hist="orange",breaks.hist=50,label.x="Age (yb2k)",
                                             label.main=NULL),
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
  ageref = x$data$y; z = x$data$z; xx = x$data$x; n=length(ageref)
  dy=x$data$dy
  reference.label = x$.args$reference.label
  if(is.null(reference.label)){
    reference.label = "reference"
  }
  eventindexes = x$.args$eventindexes; nevents = length(eventindexes)
  oldpar = par()

  if(!is.null(plot.proxydata)){
    figure.count <- new.plot(postscript,pdf,prefix,figure.count,...) +1
    par(mfrow=c(1,1),mar=(c(5,4.5,4,2)+0.1))
    if(plot.proxydata$age){
      xdata = ageref; xlab=paste0(reference.label," (yr b2k)")
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

  if(!is.null(plot.inlasims) && !is.null(x$simulation$age)){

    xlim=range(z)
    if(plot.inlasims$xrev) xlim=rev(xlim)

    ylim=range(x$simulation$summary$lower-ageref,x$simulation$summary$upper-ageref)*1.1
    if(sum(is.na(ylim))==0){
      figure.count <- new.plot(postscript,pdf,prefix,figure.count,...) +1


      nsims = plot.inlasims$nsims

      gicc = x$data$y

      fullpd = data.frame(depth = x$data$z, medians=x$simulation$summary$median-gicc,
                          lower=x$simulation$summary$lower-gicc,
                          upper=x$simulation$summary$upper-gicc)
      gg2 = ggplot(data=fullpd,aes(x=depth)) +
        geom_line(aes(y=0),color="blue",linetype="dotted",size=0.2)+
        geom_line(aes(y=medians))+
        geom_ribbon(aes(ymin=lower,ymax=upper),color="red",fill="red",alpha=0.3)+
        theme_bw()+ylab("Estimated age - GICC05 (years)")+
        xlab("NGRIP depth (m)")

      if(!is.null(nsims) && nsims>0){
        for(i in 1:min(nsims,50)){
          gg2 = gg2 + geom_line(data=data.frame(depth=x$data$z,
                                                sim = x$simulation$age[,i]-x$data$y),
                                aes(depth,sim),color="gray",alpha=0.5)
        }
      }
      print(gg2)

      if(postscript || pdf){
        if (names(dev.cur()) != "null device") {
          dev.off()
        }
      }
    }

  }

  if(!is.null(plot.syncsims) && !is.null(x$simulation$age_sync)){

    xlim=range(z)
    if(plot.inlasims$xrev) xlim=rev(xlim)

    ylim=range(x$simulation$summary$lower-ageref,x$simulation$summary$upper-ageref)*1.1
    if(sum(is.na(ylim))==0){
      figure.count <- new.plot(postscript,pdf,prefix,figure.count,...) +1


      nsims = plot.syncsims$nsims

      free_indexes = x$tie_points$free_indexes
      tie_indexes = x$tie_points$tie_indexes
      gicc_free = x$data$y[free_indexes]
      gicc_tie = x$data$y[tie_indexes]
      fullpd = data.frame(depth = x$data$z[free_indexes], medians=x$simulation$summary_sync$median[free_indexes]-gicc_free,
                          lower=x$simulation$summary_sync$lower[free_indexes]-gicc_free,
                          upper=x$simulation$summary_sync$upper[free_indexes]-gicc_free)
      gg2 = ggplot(data=fullpd,aes(x=depth)) +
        geom_line(aes(y=0),color="blue",linetype="dotted",size=0.2)+
        geom_line(aes(y=medians))+
        geom_ribbon(aes(ymin=lower,ymax=upper),color="red",fill="red",alpha=0.3)+
        theme_bw()+ylab("Estimated age - GICC05 (years)")+
        xlab("NGRIP depth (m)")
      if(x$simulation$summary_sync$lower[1] != x$simulation$summary_sync$upper[1]){
        gg2 = gg2 + geom_segment(data=data.frame(depth=x$data$z[tie_indexes],
                                                 median=x$simulation$summary_sync$median[tie_indexes]-gicc_tie,
                                                 lower=x$simulation$summary_sync$lower[tie_indexes]-gicc_tie,
                                                 upper=x$simulation$summary_sync$upper[tie_indexes]-gicc_tie),
                                 aes(x=depth,y=lower,xend=depth,yend=upper),
                                 col="magenta")
      }else{
        gg2 = gg2 + geom_point(data=data.frame(tiedepths = x$data$z[tie_indexes],
                                               tiemid = x$simulation$summary_sync$median),
                               aes(x=tiedepths,y=tiemid),
                               color="magenta")
      }
      if(!is.null(nsims) && nsims>0){
        for(i in 1:min(nsims,50)){
          gg2 = gg2 + geom_line(data=data.frame(depth=x$data$z,
                                                sim = x$simulation$age_sync[,i]-x$data$y),
                                aes(depth,sim),color="gray",alpha=0.5)
        }
      }
      print(gg2)

      if(postscript || pdf){
        if (names(dev.cur()) != "null device") {
          dev.off()
        }
      }
    }

  }

  if(!is.null(plot.tiepoints) && !is.null(x$tie_points)){
    figure.count <- new.plot(postscript,pdf,prefix,figure.count,...) +1



    if(x$tie_points$tie_n == 5){
      la=layout(mat=matrix(c(1,4,1,4,2,4,2,5,3,5,3,5) ,nrow=2))


    }else if(x$tie_points$tie_n==1){
      par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)
    }else if(x$tie_points$tie_n==2){
      par(mfrow=c(1,2),mar=c(5,4,4,2)+0.1)
    }else if(x$tie_points$tie_n==3){
      la=layout(mat=matrix(c(1,3,1,3,2,3,2,3) ,nrow=2))
    }else if(x$tie_points$tie_n==4){
      par(mfrow=c(2,2),mar=c(3,2,2,1)+0.1)
    }else if(x$tie_points$tie_n==6){
      par(mfrow=c(2,3),mar=c(5,4,4,2)+0.1)
    }

    if(is.null(plot.tiepoints$label.main)){
      labs =c()
      for(i in 1:x$tie_points$tie_n){
        labs = c(labs,paste0("Tie-point #",i))
      }
    }else{
      labs = plot.tiepoints$label.main
    }
    for(i in 1:min(x$tie_points$tie_n,6)){
      if(!is.null(x$tie_points$x.ref[i])){
        xlim=range(x$tie_points$samples[,i],x$tie_points$x.ref[i])
      }else{
        xlim=range(x$tie_points$samples[,i])
      }
      hist(x$tie_points$samples[,i],freq=0,col=plot.tiepoints$col.hist,
           breaks=plot.tiepoints$breaks.hist,xlab=plot.tiepoints$label.x,
           main=labs[i],xlim=xlim)
      if(!is.null(x$tie_points$x.ref[i])){
        abline(v=x$tie_points$x.ref[i],col="blue",lwd=2,lty=3)
      }
    }

    par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)



    if(postscript || pdf){
      if (names(dev.cur()) != "null device") {
        dev.off()
      }
    }
  }



  if(!is.null(plot.bias) && !is.null(x$biases)){
    figure.count <- new.plot(postscript,pdf,prefix,figure.count,...) +1
    nbiases = x$biases$.args$nbiases

    yrange = c(0,0)
    for( iter in 1:nbiases){
      yrange = range(yrange,x$biases[[paste0("bias",iter)]]$quant0.975-ageref,x$biases[[paste0("bias",iter)]]$quant0.025-ageref)
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
    plot(NA, xlim=xlim,xlab="Depth (m)", ylab=paste0("Simulated timescale - ",reference.label," (years)"),type="l",col="blue",ylim=yrange,main=plot.bias$label)
    abline(h=0,lty=3,lwd=1,col="gray")
    for(iter in 1:nbiases){
      lines(z,x$biases[[paste0("bias",iter)]]$quant0.975-ageref,col="blue",lty=iter)
      lines(z,x$biases[[paste0("bias",iter)]]$quant0.025-ageref,col="blue",lty=iter)
    }
    if(!is.null(plot.bias$MCE)){
      lines(z,plot.bias$MCE,lwd=2)
      lines(z,-plot.bias$MCE,lwd=2)
    }
    if(!is.null(plot.bias$legend)){
      leg=plot.bias$legend
      legend(leg$x,leg$y,leg$legend,col=leg$col,lty=leg$lty,cex=leg$cex,pch=leg$pch,lwd=leg$lwd,pt.cex=leg$pt.cex)
    }
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
