#' Print bremla model
#'
#' Prints essentials from bremla S3 class.
#'
#' @param x \code{bremla} S3 class. Output of \code{\link{bremla}} function
#' @param digits Number of digits displayed.
#' @param ... Other arguments
#'
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{bremla},\link{bremla_chronology_simulation}}
#' @keywords bremla sumary
#'
#' @examples
#' \donttest{
#' require(stats)
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
#' }
#' @export
#' @method print bremla
print.bremla = function(x,
                        digits=4L,
                        ...){

  if(!is.null(x$.args$call)){
    cat("Call:\n")
    cat(deparse(x$.args$call),"\n\n",sep="")
    cat("Time used:\n",sep="")
  }


  if(!is.null(x$fitting)){
    cpu = as.numeric(round(x$time$fit$total,digits=digits))
    cpu.navn="Model fitting"
  }else{
    cat("bremla object has been prepared, but no analysis has been performed yet.")
  }

  if(!is.null(x$simulation)){
    cpu=as.numeric(c(cpu,round(x$time$simulation$total,digits=digits)))
    cpu.navn=c(cpu.navn,"Chron. sampling")
  }
  if(!is.null(x$linramp) && !is.null(x$event_dating)){
    cpu=as.numeric(c(cpu,round(x$time$t1_and_ramp + x$time$event_age$total,digits=digits)))
    cpu.navn=c(cpu.navn,"Event dating")
  }
  if(!is.null(x$biases)){
    cpu=as.numeric(c(cpu,round(x$time$biases,digits=digits)))
    cpu.navn=c(cpu.navn,"Bias sampling")
  }

  if(!is.null(x$fitting)){
    if(!is.null(x$time$total)){
      cpu=as.numeric(c(cpu,round(x$time$total,digits=digits)))
      cpu.navn=c(cpu.navn,"Total")
    }

    names(cpu)=cpu.navn

    print(cpu)
  }




  if(!is.null(x$.args$noise)){
    if(tolower(x$.args$noise) %in% c(0,"iid","independent")){
      noise = "iid"
    }else if(tolower(x$.args$noise) %in% c(1,"ar1","ar(1)")){
      noise = "ar1"
    }else if(tolower(x$.args$noise) %in% c(2,"ar2","ar(2)")){
      noise = "ar2"
    }
  }



    if(!is.null(x$fitting)){
      if(noise == "iid"){

        hypers = matrix(round(as.numeric(x$fitting$inla$hyperparameters$results$sigma_epsilon),digits=digits),nrow=1)
        #hypers = matrix( round(as.numeric(x$fitting$hyperparameters$results$sigma_epsilon),digits=digits),ncol=mm )
        hypers = as.data.frame(hypers)
        colnames(hypers) = c("mean","sd","quant0.025","quant0.25","quant0.5","quant0.75","quant0.975")
        rownames(hypers) = c("sigma_epsilon")
      }else if(noise == "ar1"){
        hypers = rbind(round(as.numeric(x$fitting$inla$hyperparameters$results$sigma_epsilon),digits=digits),
                       round(as.numeric(x$fitting$inla$hyperparameters$results$phi),digits=digits))
        hypers = as.data.frame(hypers)
        colnames(hypers) = c("mean","sd","quant0.025","quant0.25","quant0.5","quant0.75","quant0.975")
        rownames(hypers) = c("sigma_epsilon","phi")
      }else if(noise == "ar2"){
        hypers = rbind(round(as.numeric(x$fitting$inla$hyperparameters$results$sigma_epsilon),digits=digits),
                       round(as.numeric(x$fitting$inla$hyperparameters$results$phi1),digits=digits),
                       round(as.numeric(x$fitting$inla$hyperparameters$results$phi2),digits=digits))
        hypers = as.data.frame(hypers)
        colnames(hypers) = c("mean","sd","quant0.025","quant0.25","quant0.5","quant0.75","quant0.975")
        rownames(hypers) = c("sigma_epsilon","phi1","phi2")
      }
      formulastring = format(x$.args$formula.input)

      if(!is.null(x$.args$events)){
        if(x$.args$events$fill_formula){
          formulastring = paste0(formulastring," + psi_fill(degree=",x$.args$events$degree,
                                 ", n_events=",x$.args$events$nevents,")")
        }
      }
      # maxlength=2048L
      # formulastring = format(x$.args$formula.ls)
      # if(sum(nchar(formulastring)) > maxlength){
      #   formulastring = paste0( substr(deparse(x$.args$formula.ls),1L,maxlength),"...")
      # }
      cat("The fixed component is explained by linear predictor: \n",formulastring,
          "\n\nThe noise component is explained by an ",noise," process.\n\n",sep="")

      print(hypers)
    }

  if(!is.null(x$simulation)){
    nsims = dim(x$simulation$age)[2]
    cat("\nGenerated ",nsims," chronologies.\n",sep="")
  }

    if(!is.null(x$tie_points)){
      nsims = dim(x$tie_points$nsims)
      if(tolower(x$tie_points$method) %in% c("adolphi")){
        cat("\nGenerated ",nsims," samples from ",x$tie_points$tie_n," Adolphi tie-point distributions.\n",sep="")
      }else if(tolower(x$tie_points$method) %in% c("gauss","gaussian","normal")){
        cat("\nGenerated ",nsims," samples from ",x$tie_points$tie_n," normal tie-point distributions.\n",sep="")
      }else if(tolower(x$tie_points$method) %in% c("precomputed","given")){
        cat("\nUsing ",nsims," precomputed tie-point samples.\n",sep="")
      }else if(tolower(x$tie_points$method) %in% c("semigauss","semigaussian","skewered","skewered-gauss",
                                                   "semi-gauss","semi-gaussian","merged","merged-normal")){
        cat("\nGenerated ",nsims," samples from ",x$tie_points$tie_n," merged-normal tie-point distributions.\n",sep="")
      }

    }
    if(!is.null(x$simulation$age_sync)){
      nsims = dim(x$simulation$age_sync)[2]
      cat("\nGenerated ",nsims," synchronized chronologies.\n",sep="")
    }

  if(!is.null(x$linramp)){
    hyperramp = matrix(round(c(x$linramp$param$t0$mean,x$linramp$param$t0$sd,x$linramp$param$t0$q0.025,x$linramp$param$t0$q0.5,x$linramp$param$t0$q0.975),digits=digits),nrow=1)

    rownames(hyperramp) = "onset depth"
    hyperramp = as.data.frame(hyperramp)
    lab = x$linramp$.args$label
    if(is.null(lab)) lab="unlabeled"
    cat("\nEstimated onset depth and age for event ",lab,": \n",sep="")
    if(!is.null(x$event_dating)){
      event_age = matrix(round(c(x$event_dating$mean,x$event_dating$sd,x$event_dating$q0.025,x$event_dating$q0.5,x$event_dating$q0.975),digits=digits),nrow=1)
      rownames(event_age) = "onset age"
      hyperramp = rbind(hyperramp,event_age)

    }
    colnames(hyperramp) = c("mean","sd","quant0.025","quant0.5","quant0.975")
    print(hyperramp)
    }

    if(!is.null(x$biases)){
      nbiases = x$biases$.args$nbiases
        biasparam = t(matrix(x$biases$.args$biasparam,ncol=nbiases))
        colnames(biasparam) = c("param1","param2")
        rownames(biasparam) = 1:nbiases
        biasparam = as.data.frame(biasparam)

      cat("\nGenerated ",x$biases$.args$nsims," samples for ",nbiases, " sets of (",x$biases$.args$bias.model,") biased chronologies, with parameters:\n",sep="")
      print(biasparam)
    }


  return(invisible(x))
}
