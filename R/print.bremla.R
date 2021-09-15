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
#'
#' @export
#'
print.bremla = function(x,digits=4L,...){

  cat("Call:\n")
  cat(deparse(x$.args$call),"\n\n",sep="")
  cat("Time used:\n",sep="")

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

    print(cpu)




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
      cat("The fixed component is explained by linear predictor: \n",formulastring,"\n\nThe noise component is explained by an ",noise," process.\n\n",sep="")

      print(hypers)
    }

  if(!is.null(object$simulation)){
    nsims = dim(object$simulation$age)[2]
    cat("\nGenerated ",nsims," chronologies.\n",sep="")
  }

  if(!is.null(object$linramp)){
    hyperramp = matrix(round(c(object$linramp$param$t0$mean,object$linramp$param$t0$sd,object$linramp$param$t0$q0.025,object$linramp$param$t0$q0.5,object$linramp$param$t0$q0.975),digits=digits),nrow=1)

    rownames(hyperramp) = "onset depth"
    hyperramp = as.data.frame(hyperramp)
    lab = object$linramp$.args$label
    if(is.null(lab)) lab="unlabeled"
    cat("\nEstimated onset depth and age for event ",lab,": \n",sep="")
    if(!is.null(object$DO_dating)){
      DO_age = matrix(round(c(object$DO_dating$mean,object$DO_dating$sd,object$DO_dating$q0.025,object$DO_dating$q0.5,object$DO_dating$q0.975),digits=digits),nrow=1)
      rownames(DO_age) = "onset age"
      hyperramp = rbind(hyperramp,DO_age)

    }
    colnames(hyperramp) = c("mean","sd","quant0.025","quant0.5","quant0.975")
    print(hyperramp)
    }

    if(!is.null(object$biases)){
        biasparam = t(matrix(object$biases$.args$biasparam,ncol=nbiases))
        colnames(biasparam) = c("param1","param2")
        rownames(biasparam) = 1:nbiases
        biasparam = as.data.frame(biasparam)

      cat("\nGenerated ",object$biases$.args$nsims," samples for ",object$biases$.args$nbiases, " sets of (",object$biases$.args$bias.model,") biased chronologies, with parameters:\n",sep="")
      print(biasparam)
    }


  return(invisible(x))
}
