#' Summary bremla model
#'
#' Summarizes results from bremla S3 class.
#'
#' @param object \code{bremla} S3 class. Output of \code{\link{bremla}} function
#' @param digits Number of digits displayed.
#' @param ... Other arguments
#'
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{bremla},\link{bremla_chronology_simulation}}
#' @keywords bremla summary
#'
#' @examples
#' \donttest{
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
#' object = bremla(age,depth,proxy,events=eventdepths,nsims=10000,
#'   synchronization = list(locations=c(11050,12050,13050,22050,42050),
#'                           locations.type="age",method="adolphi",
#'                           samples=NULL),
#'   print.progress=TRUE
#'   )
#'   summary(object)
#' }
#'
#' @export
#' @method summary bremla
summary.bremla = function(object,
                          digits=4L,
                          ...){

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
  if(!is.null(object$linramp) && !is.null(object$event_dating)){
    cpu=as.numeric(c(cpu,round(object$time$t1_and_ramp + object$time$event_age$total,digits=digits)))
    cpu.navn=c(cpu.navn,"event dating")
  }
  if(!is.null(object$biases)){
    cpu=as.numeric(c(cpu,round(object$time$biases,digits=digits)))
    cpu.navn=c(cpu.navn,"Bias sampling")
  }
  if(!is.null(object$fitting)){
    if(!is.null(object$time$total)){
      cpu=as.numeric(c(cpu,round(object$time$total,digits=digits)))
      cpu.navn=c(cpu.navn,"Total")
    }

    names(cpu)=cpu.navn
    ut=c(ut, list(cpu.used=cpu))
  }


  if(tolower(object$.args$control.fit$noise) %in% c(0,"iid","independent")){
    noise = "iid"
  }else if(tolower(object$.args$control.fit$noise) %in% c(1,"ar1","ar(1)")){
    noise = "ar1"
  }else if(tolower(object$.args$control.fit$noise) %in% c(2,"ar2","ar(2)")){
    noise = "ar2"
  }

  if(!is.null(object$fitting)){
    if(noise == "iid"){

      hypers = matrix(round(as.numeric(object$fitting$inla$hyperparameters$results$sigma_epsilon),digits=digits),nrow=1)
      #hypers = matrix( round(as.numeric(object$fitting$hyperparameters$results$sigma_epsilon),digits=digits),ncol=mm )
      hypers = as.data.frame(hypers)
      colnames(hypers) = c("mean","sd","quant0.025","quant0.25","quant0.5","quant0.75","quant0.975")
      rownames(hypers) = c("sigma_epsilon")
    }else if(noise == "ar1"){
      hypers = rbind(round(as.numeric(object$fitting$inla$hyperparameters$results$sigma_epsilon),digits=digits),
                     round(as.numeric(object$fitting$inla$hyperparameters$results$phi),digits=digits))
      hypers = as.data.frame(hypers)
      colnames(hypers) = c("mean","sd","quant0.025","quant0.25","quant0.5","quant0.75","quant0.975")
      rownames(hypers) = c("sigma_epsilon","phi")
    }else if(noise == "ar2"){
      hypers = rbind(round(as.numeric(object$fitting$inla$hyperparameters$results$sigma_epsilon),digits=digits),
                     round(as.numeric(object$fitting$inla$hyperparameters$results$phi1),digits=digits),
                     round(as.numeric(object$fitting$inla$hyperparameters$results$phi2),digits=digits))
      hypers = as.data.frame(hypers)
      colnames(hypers) = c("mean","sd","quant0.025","quant0.25","quant0.5","quant0.75","quant0.975")
      rownames(hypers) = c("sigma_epsilon","phi1","phi2")
    }
    maxlength=2048L
    formulastring = format(object$.args$formula.ls)
    if(sum(nchar(formulastring)) > maxlength){
      formulastring = paste0( substr(deparse(object$.args$formula.ls),1L,maxlength),"...")
    }
    fit.arg = list(formula = formulastring,noise=object$.args$control.fit$noise,
                   nevents=object$.args$events$nevents,
                   method = object$.args$control.fit$method)

    ut=c(ut, list(hyperpar=hypers,fit.arg=fit.arg))
  }



  if(!is.null(object$simulation)){
    sim = list(nsims = object$.args$control.sim$nsims,
               n = nrow(object$data),
               syncstatus = object$.args$control.sim$synchronized,
               store.means = !is.null(object$simulation$dmean) )
    ut = c(ut,sim)
  }


  if(!is.null(object$tie_points)){
    tiepoints = list(tie_n=object$tie_points$tie_n,
                     free_n=object$tie_points$free_n,
                     method=object$tie_points$method,
                     nsims=object$tie_points$nsims,
                     locations=object$tie_points$locations,
                     locations.type=object$tie_points$unit)
    ut = c(ut,list(tiepoints=tiepoints))
  }

  if(!is.null(object$linramp)){
    hyperramp = rbind(round(c(object$linramp$param$t0$mean,
                              object$linramp$param$t0$sd,
                              object$linramp$param$t0$q0.025,
                              object$linramp$param$t0$q0.5,
                              object$linramp$param$t0$q0.975),
                            digits=digits),
                      round(c(object$linramp$param$dtpos$mean,
                              object$linramp$param$dtpos$sd,
                              object$linramp$param$dtpos$q0.025,
                              object$linramp$param$dtpos$q0.5,
                              object$linramp$param$dtpos$q0.975),
                            digits=digits),
                      round(c(object$linramp$param$y0$mean,
                              object$linramp$param$y0$sd,
                              object$linramp$param$y0$q0.025,
                              object$linramp$param$y0$q0.5,
                              object$linramp$param$y0$q0.975),
                            digits=digits),
                      round(c(object$linramp$param$dy$mean,
                              object$linramp$param$dy$sd,
                              object$linramp$param$dy$q0.025,
                              object$linramp$param$dy$q0.5,
                              object$linramp$param$dy$q0.975),
                            digits=digits)
                      )
    colnames(hyperramp) = c("mean","sd","quant0.025","quant0.5","quant0.975")
    if(object$linramp$.args$nsims>0){
      hyperramp = rbind(hyperramp,
                        c(object$linramp$param$t1$mean,
                          object$linramp$param$t1$sd,
                          object$linramp$param$t1$q0.025,
                          object$linramp$param$t1$q0.5,
                          object$linramp$param$t1$q0.975))
      rownames(hyperramp) = c("t0", "dt", "y0","dy","t1")
    }else{
      rownames(hyperramp) = c("t0", "dt", "y0","dy")
    }
    hyperramp = as.data.frame(round(hyperramp,digits=digits))
    linramplist = list(hyperramp = hyperramp,t1sims=object$linramp$.args$nsims,
                       rampsims = object$linramp$.args$nsims,
                       label=object$linramp$.args$label,
                       depth.reference=object$linramp$.args$depth.reference)
    ut = c(ut,linramplist)
  }

  if(!is.null(object$event_dating)){
    event_age = matrix(round(c(object$event_dating$mean,
                            object$event_dating$sd,
                            object$event_dating$q0.025,
                            object$event_dating$q0.5,
                            object$event_dating$q0.975),
                          digits=digits),
                    nrow=1)
    colnames(event_age) = c("mean","sd","quant0.025","quant0.5","quant0.975")
    rownames(event_age) = "Onset age"
    event_age = as.data.frame(round(event_age,digits=digits))
    eventlist = list(event_age=event_age,datingsims = object$event_dating$.args$nsims,label=object$event_dating$.args$label,age.reference=object$event_dating$.args$age.reference)
    ut = c(ut,eventlist)
  }

  if(!is.null(object$biases)){
    nbiases = object$biases$.args$nbiases

    # if(nbiases == 1){
    #   #biasparam = object$biases$.args$biasparam
    #   biasparam = matrix(object$biases$.args$biasparam,nrow=1)
    #   colnames(biasparam) = c("param1","param2")
    #   rownames(biasparam) = ""
    #   biasparam = as.data.frame(biasparam)
    # }else{
      biasparam = t(matrix(object$biases$.args$biasparam,ncol=nbiases))
      colnames(biasparam) = c("param1","param2")
      rownames(biasparam) = 1:nbiases
      biasparam = as.data.frame(biasparam)
    #}

    biaslist = list(nbiases = nbiases, bias.model = object$biases$.args$bias.model, biasparam = biasparam, store.samples = object$biases$.args$store.samples,biasnsims=object$biases$.args$nsims)
    ut = c(ut,biaslist)
  }
    ut = c(ut,reference.label=list(reference.label=object$.args$reference.label))

  class(ut) = "summary.bremla"

return(ut)
}

#' Print summary.bremla class
#'
#' Prints summary.bremla S3 class.
#'
#' @param x \code{summary.bremla} S3 class. Output of \code{\link{summary.bremla}} function
#' @param digits Number of digits displayed.
#' @param ... Other arguments
#'
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{bremla},\link{bremla_chronology_simulation}}
#' @keywords bremla summary.print
#'
#' @examples
#' \donttest{
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

#' object = bremla(age,depth,proxy,events=eventdepths,nsims=100,
#'   synchronization = list(locations=c(11050,12050,13050,22050,42050),
#'                           locations.type="age",method="adolphi",
#'                           samples=NULL),
#'   print.progress=TRUE
#'   )
#' ss = summary(object)
#' print(ss)
#' }
#'
#' @export
#' @method print summary.bremla
print.summary.bremla = function(x,
                                digits=4L,
                                ...){
  if(!is.null(x$call)) cat("\nCall:\n",deparse(x$call),"\n\n",sep="")

  if(!is.null(x$cpu.used)){
    cat("Time used:\n")
    print(x$cpu.used)
    cat("\n",sep="")
  }

  if(!is.null(x$fit.arg$formula)){
    cat("The fixed component is explained by linear predictor: \n",x$fit.arg$formula,"\n\nThe noise component is explained by an ",x$fit.arg$noise," process.\n",sep="")
  }else{
    cat("bremla object is prepared, but no analysis has been performed yet.")
  }

  if(!is.null(x$fit.arg$method)){
    if(tolower(x$fit.arg$method) %in% c("inla")){
      cat("\nThe model is fitted using INLA, with following estimates for the hyperparameters:\n")
      print(format(x$hyperpar,digits=digits,nsmall=2),quote=FALSE)
    }else{
      cat("\nThe model is fitted using least squares, with following estimates for the model parameters:\n")
      print(format(x$hyperpar,digits=digits,nsmall=2),quote=FALSE)
    }
  }


  if(!is.null(x$nsims)){
    if(!is.null(x$reference.label)){
      cat("\nSimulating ",x$nsims," chronologies, using ",x$reference.label," as reference.\n",sep="")
    }else{
      cat("\nSimulating ",x$nsims," chronologies.\n",sep="")
    }

    if(x$store.means){
      cat(" Storing means\n")
    }
    #cat("\n")
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

    if(!is.null(x$tiepoints)){
      cat("\n",x$tiepoints$nsims, " chronologies sampled using ", x$tiepoints$tie_n ,
          " tie-point distributions",sep="")
      if(tolower(x$tiepoints$method) %in% c("adolphi")){
        cat(" (Adolphi).",sep="")
      }else if(tolower(x$tiepoints$method) %in% c("precomputed", "given")){
        cat(" (precomputed).",sep="")
      }else if(tolower(x$tiepoints$method) %in% c("normal", "gauss","gaussian")){
        cat(" (Gaussian).",sep="")
      }else if(tolower(x$tiepoints$method) %in% c("semigauss","semi-gauss","quasigauss","quasi-gauss",
                                                  "skewered","skewered-gauss")){
        cat(" (merged Gaussians).",sep="")
      }else{cat(".")}

      if(x$tiepoints$tie_n <= 6){
        cat("\nTie-points are fixed at ")
        if(tolower(x$tiepoints$locations.type) %in% c("depth","z")){
          cat("NGRIP depths (m) ",sep="")
          cat(x$tiepoints$locations, sep=", ")
          cat(".")
        }else if(tolower(x$tiepoints$locations.type) %in% c("age","y")){
          cat(x$reference.label," ages (yb2k) ",sep="")
          cat(x$tiepoints$locations, sep=", ")
          cat(".")
        }else if(tolower(x$tiepoints$locations.type) %in% c("index","ind")){
          cat("indices ",sep="")
          cat(x$tiepoints$locations, sep=", ")
          cat(".")
        }
        cat("\n",sep="")
      }

    }

    if(x$rampsims>0) cat("\n",x$rampsims, " samples of linear ramp function produced.\n",sep="")

    if(!is.null(x$age.reference)){
      if(!is.null(x$event_age)) cat("\nGenerated ",x$datingsims, " samples of onset ages (reference onset age is ",x$age.reference,").\n",sep="")
    }else{
      if(!is.null(x$event_age)) cat("\nGenerated ",x$datingsims, " samples of onset ages.\n",sep="")
    }

    if(!is.null(x$event_age)) print(x$event_age)
    #if(!is.null(x$age.reference)) cat("The reference onset age is ",x$age.reference,"\n",sep="")
  }

  if(!is.null(x$biasparam)){
    cat("\nGenerated ",x$biasnsims," samples for ",x$nbiases, " sets of (",x$bias.model,") biased chronologies, with parameters:\n",sep="")
    print(x$biasparam)
  }



}
