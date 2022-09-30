#' Plot bremla model
#'
#' Plots results from bremla S3 class.
#'
#' @param x \code{bremla} S3 class. Output of \code{\link{bremla}} function
#' @param plot.ls List specifying how the least square fit should be illustrated. \code{fitted=TRUE} plots the fitted values with label \code{label.fit}, \code{legend} specifies the legend, \code{residuals=TRUE} means the residuals are plotted with title \code{label.res}.
#' if \code{histogram=TRUE} the residuals are represented in a histogram with label \code{label.hist}, if \code{qqplot=TRUE} a quantile-quantile plot of the residuals are plotted with title \code{label.qq} and if
#' \code{acf=TRUE} the empirical autocorrelations are plotted with title \code{label.acf}.
#' @param plot.inla.posterior list specifying how the results from the inla regression fit should be plotted. If \code{posteriors=TRUE} then the posterior marginal distributions of the hyperparameters are plotted with title \code{label}.
#' @param plot.inlasims list specifying how the simulated chronologies from the INLA posterior should be plotted. \code{plotsims} gives how many simulated chronologies should be included in the plot (with title \code{label}), \code{legend} specifies the legend, if \code{xrev=TRUE} the x-axis is reversed to give a chronological ordering.
#' @param plot.syncsims list specifying how the simulated synchronized chronologies from the INLA posterior should be plotted. \code{plotsims} gives how many simulated chronologies should be included in the plot (with title \code{label}), \code{legend} specifies the legend, if \code{xrev=TRUE} the x-axis is reversed to give a chronological ordering.
#' @param plot.tiepoints list containing specifications for plotting histogram of sampled tie-points.
#' @param plot.bias list specifying how the simulations under the assumptions of unknown counting bias should be represented. If \code{MCE} is given as a numeric vector it will be included in the plot (with title \code{label}). \code{legend} specifies the legend, if \code{xrev=TRUE} the x-axis is reversed to give a chronological ordering.
#' @param plot.linramp list specifying how the linear ramp fit should be plotted. If \code{depth.reference} is given, it will be represented by a vertical dotted line. If \code{show.t0=TRUE} the posterior marginal distribution of the onset is included (non-normalized). If \code{show.t1=TRUE} the posterior marginal distribution of the end point of the transition is included (non-normalized). If \code{xrev=TRUE} the x-axis is reversed to give a chronological ordering. \code{label} gives the title of the plot.
#' @param plot.event_depth list specifying how the posterior distribution of the onset depth should be plotted. If \code{depth.reference} is given, it will be represented by a vertical dotted line. If \code{xrev=TRUE} the x-axis is reversed to give a chronological ordering. \code{label} gives the title of the plot.
#' @param plot.event_age list specifying how the histogram of the simulated onset ages should be plotted. If \code{age.reference} is given, it will be represented by a vertical dotted line. If \code{xrev=TRUE} the x-axis is reversed to give a chronological ordering. \code{label} gives the title of the plot.
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
#' \donttest{
#' require(stats)
#' set.seed(1)
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
#' object = bremla_chronology_simulation(object, print.progress=TRUE)
#' summary(object)
#' plot(object)
#'
#' }
#'
#' @export
#' @importFrom INLA inla inla.tmarginal inla.zmarginal inla.ar.pacf2phi
#' @importFrom grDevices dev.cur dev.new dev.off
#' @importFrom graphics abline hist legend lines par
#' @importFrom stats qqline qqnorm
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon theme_bw xlab ylab geom_segment geom_point
plot.bremla = function(x,
                       plot.ls = list(fitted=TRUE,legend=NULL,residuals=TRUE,histogram=TRUE,qqplot=FALSE,acf=TRUE,xrev=FALSE,
                                      label.fit=NULL,label.res=NULL,label.hist=NULL,label.qq=NULL,label.acf=NULL,merge=TRUE),
                       plot.inla.posterior = list(posteriors=TRUE,label=NULL),
                       plot.inlasims = list(plotsims=0,legend=NULL,xrev=FALSE,label=NULL),
                       plot.syncsims = list(plotsims=0,legend=NULL,xrev=FALSE,label=NULL),
                       plot.tiepoints = list(col.hist="orange",breaks.hist=50,label.x="Age (yb2k)",
                                             label.main=NULL),
                       plot.bias = list(MCE=NULL,legend=NULL,xrev=FALSE,label=NULL),
                       plot.linramp = list(depth.reference=NULL,show.t0=TRUE,show.t1=TRUE,xrev=TRUE,label=NULL),
                       plot.event_depth = list(depth.reference=NULL,xrev=TRUE,label=NULL),
                       plot.event_age = list(age.reference=NULL,xrev=TRUE,label=NULL),
                       postscript=FALSE,
                       pdf=FALSE,
                       prefix = "bremla.plots/figure-",
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
  ageref = x$data$age; z = x$data$depth; n=length(ageref)
  response = x$fitting$LS$fit$model[,1]
  #dy=x$data$dy
  # reference.label = x$.args$reference.label
  # if(is.null(reference.label)){
  #   reference.label = "reference"
  # }
  eventindexes = x$.args$events$eventindexes
  nevents = length(eventindexes)
  oldpar = par()




  if(!is.null(plot.ls) && !is.null(x$fitting$LS$fit)){
    figure.count <- new.plot(postscript,pdf,prefix,figure.count,...) +1
    if(plot.ls$merge){
      par(mfrow=c(2,2),mar=c(5,4,4,2)+0.1)
    }else{
      par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)
    }
    plotnum=0


    if(plot.ls$fitted){
      plotnum=plotnum+1
      xlim = range(z)
      if(plot.ls$xrev) xlim=rev(xlim)
      if(is.null(plot.ls$label.fit)){
        plot.label="Least squares fit"
      }else{
        plot.label=plot.ls$label.fit
      }
      responsename = x$.args$responsename
      response = x$data[[responsename]]
      plot(z,response,type="l",xlab=paste0("Depth (m)"),xlim=xlim,ylab=("Layers per 5 cm"),main=plot.label)
      lines(z,x$fitting$LS$fit$fitted.values,col="red")
      abline(v=z[eventindexes],col="gray",lwd=0.8)
    }
    if(!is.null(plot.ls$legend)) legend(x=leg$x,y=leg$y,legend=leg$legend,col=leg$col,lty=leg$lty,cex=leg$cex,pch=leg$pch,lwd=leg$lwd,pt.cex=leg$pt.cex,bty=leg$bty)

    if(plot.ls$residuals){
      plotnum=plotnum+1
      if(is.null(plot.ls$label.res)){
        plot.label="Least squares residuals"
      }else{
        plot.label=plot.ls$label.res
      }
      plot(z,x$fitting$LS$fit$residuals,type="l",xlab="Depth (m)",ylab="Residual errors per 5 cm",main=plot.label,xlim=xlim); abline(h=0,lty=3,col="gray")
    }
    if(plot.ls$histogram){
      plotnum=plotnum+1
      if(is.null(plot.ls$label.hist)){
        plot.label="Histogram: Residuals"
      }else{
        plot.label=plot.ls$label.hist
      }
      hist(x$fitting$LS$fit$residuals,freq=0,col="orange",breaks=20,xlab="Residual errors per 5 cm", main=plot.label)
    }
    if(plot.ls$qqplot){
      plotnum=plotnum+1
      if(is.null(plot.ls$label.qq)){
        plot.label="Q-Q Plot"
      }else{
        plot.label=plot.ls$label.qq
      }
      qqnorm(x$fitting$LS$fit$residuals,main=plot.label); qqline(x$fitting$LS$fit$residuals)
    }
    if(plotnum>4){
      par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)
    }
    if(plot.ls$acf){
      plotnum=plotnum+1
      if(is.null(plot.ls$label.acf)){
        plot.label="Autocorrelation function"
      }else{
        plot.label=plot.ls$label.acf
      }
      acf(x$fitting$LS$fit$residuals,lag.max = 30,main=plot.label,lwd=2)
    }
    if(plot.ls$merge){
      par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)
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
      plot(x$fitting$inla$hyperparameters$posteriors$sigma_epsilon,type="l",xlab=expression(paste(sigma[epsilon])),ylab="Density",lwd=2,main=plot.label1)
      abline(v=x$fitting$inla$hyperparameters$results$sigma_epsilon$mean)
      abline(v=c(x$fitting$inla$hyperparameters$results$sigma_epsilon$quant0.025,x$fitting$inla$hyperparameters$results$sigma_epsilon$quant0.975),col="gray")

      if(tolower(x$.args$control.fit$noise) %in% c(1,"ar1","ar(1)")){
        if(length(plot.inla.posterior$label)>=2){
          plot.label2 = plot.inla.posterior$label[2]
        }
        plot(x$fitting$inla$hyperparameters$posteriors$phi,xlab=expression(paste(phi)),ylab="Density",lwd=2,type="l",main=plot.label2)
        abline(v=x$fitting$inla$hyperparameters$results$phi$mean)
        abline(v=c(x$fitting$inla$hyperparameters$results$phi$quant0.025,x$fitting$inla$hyperparameters$results$phi$quant0.975),col="gray")
      }else if(tolower(x$.args$control.fit$noise) %in% c(2,"ar2","ar(2)")){
        if(length(plot.inla.posterior$label)>=3){
          plot.label2 = plot.inla.posterior$label[2]
          plot.label3 = plot.inla.posterior$label[3]
        }
        plot(x$fitting$inla$hyperparameters$posteriors$phi1,xlab=expression(paste(phi[1])),ylab="Density",lwd=2,type="l",main=plot.label2)
        abline(v=x$fitting$inla$hyperparameters$results$phi1$mean)
        abline(v=c(x$fitting$inla$hyperparameters$results$phi1$quant0.025,x$fitting$inla$hyperparameters$results$phi1$quant0.975),col="gray")

        plot(x$fitting$inla$hyperparameters$posteriors$phi2,xlab=expression(paste(phi[2])),ylab="Density",lwd=2,type="l",main=plot.label3)
        abline(v=x$fitting$inla$hyperparameters$results$phi2$mean)
        abline(v=c(x$fitting$inla$hyperparameters$results$phi2$quant0.025,x$fitting$inla$hyperparameters$results$phi2$quant0.975),col="gray")
      }
    }
    if(postscript || pdf){
      if (names(dev.cur()) != "null device") {
        dev.off()
      }
    }
  }

  if(!is.null(plot.inlasims) && !is.null(x$.args$control.sim$synchronized)){


  if(x$.args$control.sim$synchronized %in% c(FALSE,2)){

    xlim=range(z)
    if(plot.inlasims$xrev) xlim=rev(xlim)


      figure.count <- new.plot(postscript,pdf,prefix,figure.count,...) +1


      plotsims = plot.inlasims$plotsims

      ageref = x$data$age

      if(!is.null(x$simulation$summary)){

        fullpd = data.frame(depth = x$data$depth,
                            medians=x$simulation$summary$median-ageref,
                            lower=x$simulation$summary$lower-ageref,
                            upper=x$simulation$summary$upper-ageref)
      }else{
        fullpd = data.frame(depth = x$data$depth)
      }

      gg2 = ggplot(data=fullpd,aes(x=.data$depth)) +
        geom_line(aes(y=0),color="blue",linetype="dotted",size=0.2)+
        theme_bw()+ylab("Estimated age - reference (years)")+
        xlab("NGRIP depth (m)")

      if(!is.null(x$simulation$summary)){
        gg2 = gg2 +
          geom_line(aes(y=.data$medians))+
          geom_ribbon(aes(ymin=.data$lower,ymax=.data$upper),color="red",fill="red",alpha=0.3)
      }



      if(!is.null(plotsims) && plotsims>0){
        for(i in 1:min(plotsims,50)){
          gg2 = gg2 + geom_line(data=data.frame(depth=x$data$depth,
                                                sim = x$simulation$age[,i]-ageref),
                                aes(.data$depth,.data$sim),color="gray",alpha=0.5)
        }
      }
      if((!is.null(plotsims) && plotsims>0) || !is.null(x$simulation$summary) ){
        print(gg2)
      }

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
        abline(v=x$tie_points$x.ref[i],col="blue",lwd=3,lty=3)
      }
    }

    par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)



    if(postscript || pdf){
      if (names(dev.cur()) != "null device") {
        dev.off()
      }
    }
  }



  if(!is.null(plot.syncsims) && !is.null(x$simulation$age_sync)){

    xlim=range(z)
    if(plot.inlasims$xrev) xlim=rev(xlim)

    ylim=range(x$simulation$summary$lower-ageref,x$simulation$summary$upper-ageref)*1.1
    if(sum(is.na(ylim))==0){
      figure.count <- new.plot(postscript,pdf,prefix,figure.count,...) +1


      plotsims = plot.syncsims$plotsims

      free_indexes = x$tie_points$free_indexes
      tie_indexes = x$tie_points$tie_indexes
      ageref_free = x$data$age[free_indexes]
      ageref_tie = x$data$age[tie_indexes]

      if(!is.null(x$simulation$summary_sync)){
        fullpd = data.frame(depth = x$data$depth[free_indexes],
                            medians=x$simulation$summary_sync$median[free_indexes]-ageref_free,
                            lower=x$simulation$summary_sync$lower[free_indexes]-ageref_free,
                            upper=x$simulation$summary_sync$upper[free_indexes]-ageref_free)
      }else{
        fullpd = data.frame(depth = x$data$depth[free_indexes])
      }

      gg2 = ggplot(data=fullpd,aes(x=.data$depth)) +
        geom_line(aes(y=0),color="blue",linetype="dotted",size=0.2)+
        theme_bw()+ylab("Estimated age - reference (years)")+
        xlab("NGRIP depth (m)")

      if(!is.null(x$simulation$summary_sync)){
        gg2 = gg2+geom_line(aes(y=.data$medians))+
          geom_ribbon(aes(ymin=.data$lower,ymax=.data$upper),color="red",fill="red",alpha=0.3)
      }

      if(!is.null(x$simulation$summary_sync)){
        if(x$simulation$summary_sync$lower[tie_indexes[1]] != x$simulation$summary_sync$upper[tie_indexes[1]]){
          gg2 = gg2 + geom_segment(data=data.frame(depth=x$data$depth[tie_indexes],
                                                   median=x$simulation$summary_sync$median[tie_indexes]-ageref_tie,
                                                   lower=x$simulation$summary_sync$lower[tie_indexes]-ageref_tie,
                                                   upper=x$simulation$summary_sync$upper[tie_indexes]-ageref_tie),
                                   aes(x=.data$depth,y=.data$lower,xend=.data$depth,yend=.data$upper),
                                   col="magenta")
        }else{
          gg2 = gg2 + geom_point(data=data.frame(tiedepths = x$data$depth[tie_indexes],
                                                 tiemid = x$simulation$summary_sync$median),
                                 aes(x=.data$tiedepths,y=.data$tiemid),
                                 color="magenta")
        }
      }

      if(!is.null(plotsims) && plotsims>0){
        for(i in 1:min(plotsims,50)){
          gg2 = gg2 + geom_line(data=data.frame(depth=x$data$depth,
                                                sim = x$simulation$age_sync[,i]-ageref),
                                aes(.data$depth,.data$sim),color="gray",alpha=0.5)
        }
      }
      if((!is.null(plotsims) && plotsims>0) || !is.null(x$simulation$summary_sync) ){
        print(gg2)
      }


      if(postscript || pdf){
        if (names(dev.cur()) != "null device") {
          dev.off()
        }
      }
    }

  }


  if(!is.null(plot.bias) && !is.null(x$biases)){
    figure.count <- new.plot(postscript,pdf,prefix,figure.count,...) +1
    nbiases = x$biases$.args$nbiases

    yrange = c(0,0)
    for( iter in 1:nbiases){
      yrange = range(yrange,x$biases[[paste0("bias",iter)]]$quant0.975-ageref,
                     x$biases[[paste0("bias",iter)]]$quant0.025-ageref)
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
    if(is.null(x$.args$reference.label)){
      plot(NA, xlim=xlim,xlab="Depth (m)", ylab=paste0("Simulated timescale - reference timescale (years)"),type="l",col="blue",ylim=yrange,main=plot.bias$label)
    }else{
      plot(NA, xlim=xlim,xlab="Depth (m)", ylab=paste0("Simulated timescale - ",x$.args$reference.label," (years)"),type="l",col="blue",ylim=yrange,main=plot.bias$label)
    }

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

    ybottom = min(x$linramp$data$y)
    ytop = min(x$linramp$linrampfit$q0.025)
    epsilon = (ytop-ybottom)*0.3
    if(plot.linramp$show.t0 && !is.null(x$linramp$param)){

      margt0=x$linramp$param$t0$marg.t0
      #normt0.y = margt0[,2]/diff(range(margt0[,2]))*(ytop-ybottom)-min(margt0[,2])+ybottom
      normt0.y = (margt0[,2]-min(margt0[,2]))/diff(range(margt0[,2]))*(ytop-ybottom)+ybottom-epsilon
    }
    if(plot.linramp$show.t0 && !is.null(x$linramp$param$t1)){
      # ybottom = min(x$linramp$data$y)
      # ytop = min(x$linramp$linrampfit$q0.025)
      margt1 = x$linramp$param$t1$marginal
      normt1.y = margt1[,2]/diff(range(margt1[,2]))*(ytop-ybottom)-min(margt1[,2])+ybottom-epsilon
    }
    yrange2=range(yrange,normt0.y,normt1.y)
    plot(xval,x$linramp$data$y,type="l",lwd=1.25,col="gray",xlim=xlim,ylim=yrange2,xlab="Depth (m)",
         ylab=expression(paste(delta^18, "O (permil)")),main=plot.label)
    lines(xval,x$linramp$linrampfit$q0.025,col="red",lwd=2)
    lines(xval,x$linramp$linrampfit$q0.975,col="red",lwd=2)
    lines(xval,x$linramp$linrampfit$mean,col="black",lwd=2)

    lines(x=margt0[,1],y=normt0.y,col="blue",lwd=2)
    lines(x=margt1[,1],y=normt1.y,col="blue",lwd=2,lty=3)

    if(!is.null(plot.event_depth$depth.reference)){
      abline(v=plot.event_age$age.reference,lwd=2,lty=3)
    }else{
      if(!is.null(x$linramp$.args$depth.reference)){
        abline(v= x$linramp$.args$depth.reference,lwd=2,lty=3)
      }
    }
  }

  if(!is.null(plot.event_depth) && !is.null(plot.linramp) && !is.null(x$event_dating)){
    figure.count <- new.plot(postscript,pdf,prefix,figure.count,...) +1
    par(mfrow=c(1,1),mar=(c(5,4,4,2)+0.1))
    if(is.null(plot.event_depth$label)){
      plot.label = x$linramp$.args$label
    }else{
      plot.label = plot.event_depth$label
    }
    xlim = range(x$linramp$param$t0$marg.t0[,1])
    if(!is.null(plot.event_depth$depth.reference)){
      xlim=range(xlim,plot.event_depth$depth.reference)
    }else{
      if(!is.null(x$linramp$.args$depth.reference)){
        xlim=range(xlim,x$linramp$.args$depth.reference)
      }
    }
    if(plot.event_age$xrev) xlim=rev(xlim)
    plot(x$linramp$param$t0$marg.t0,type="l",lwd=2,xlab="Onset depth (m)",ylab="Density",main=plot.label,xlim=xlim)
    abline(v=x$linramp$param$t0$mean,lwd=2)
    abline(v=c(x$linramp$param$t0$q0.025,x$linramp$param$t0$q0.975),lwd=2,col="gray")

    if(!is.null(plot.event_depth$depth.reference)){
      abline(v=plot.event_age$age.reference,lwd=2,lty=3)
    }else{
      if(!is.null(x$linramp$.args$depth.reference)){
        abline(v= x$linramp$.args$depth.reference,lwd=2,lty=3)
      }
    }
  }

  if(!is.null(plot.event_age) && !is.null(x$event_dating)){
    xlim = range(x$event_dating$samples)
    if(plot.event_age$xrev) xlim=rev(xlim)
    hist(x$event_dating$samples,xlab="Onset age (y b2k)",freq=FALSE,col="orange",breaks=30,main=x$event_dating$.args$label,xlim=xlim)
    abline(v=x$event_dating$mean,lwd=2)
    abline(v=c(x$event_dating$q0.025,x$event_dating$q0.975),col="gray",lwd=2)
    if(!is.null(plot.event_age$age.reference)){
      abline(v=plot.event_age$age.reference,lwd=2,lty=3)
    }else{
      if(!is.null(x$event_dating$.args$age.reference)){
        abline(v=x$event_dating$.args$age.reference,lwd=2,lty=3)
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
