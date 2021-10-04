#' Bremla preparation function
#'
#' Prepares input and formats data for bremla and functions therein.
#'
#' @param age Vector of observed ages y_0,...,y_n
#' @param depth Vector of depths z_0,...,z_n corresponding to \code{age}
#' @param proxy Vector of observed proxy (e.g. del-18O) at depths \code{depth}
#' @param events Vector describing locations (in dimension specified by \code{eventmeasure}) of climate transitions for use in linear regression model
#' @param nsims Number of chronologies to be simulated
#' @param eventmeasure Character describing in which dimension the climate events are located ("age", "depth", "index")
#' @param reg.model List of booleans specifies which effects to be included in the linear regression model: intercept (\code{const}), linear response to depth (\code{depth1}), quadratic response to depth (\code{depth2}), linear response to proxy (\code{proxy}), individual intercepts for each climate period (\code{psi0}), individual linear responses to depth for each climate period (\code{psi1})
#' @param noise Character specifying the noise model: independent identically distributed (\code{iid}), first order autoregressive (\code{ar1}) or second order autoregressive (\code{ar2})
#' @param method Character specifying how the model is fitted. Currently only least squares (\code{ls}) and INLA (\code{inla}) are supported, and least squares is always run.
#'
#' @return Returns a list containing data in internal formatting (\code{data}), linear regression formula (\code{ls.formula}), inla formula (\code{formula}), input settings and indices corresponding to climate transitions (\code{.args}) and z_0, y_0 and x_0 (\code{preceeding})
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{bremla},\link{bremla_chronology_simulation}}
#' @keywords bremla preparation
#'
#' @examples
#'
#' @export
bremla_prepare = function(age,depth,proxy, events=NULL,nsims=10000, eventmeasure = "depth",reg.model = list(
  const=FALSE,depth1=FALSE,depth2=TRUE,proxy=TRUE,psi0=TRUE,psi1=TRUE), noise="ar1", method="inla",reference.label=NULL){
  n = length(age)
  y = age[2:n]
  dy = diff(age)
  z = depth[2:n]

  data = data.frame(y=y,dy=dy,z=z)

  formulastring = "dy ~ "
  if(!reg.model$const) formulastring=paste0(formulastring,"-1") else formulastring=paste0(formulastring,"1")
  if(reg.model$depth1) formulastring=paste0(formulastring," + z")
  if(reg.model$depth2) {
    formulastring=paste0(formulastring," + z2")
    data[["z2"]]=z^2
  }
  if(reg.model$proxy) {
    if(missing(proxy)){
      stop("'proxy' is missing. Shutting down...")
    }
    formulastring=paste0(formulastring," + x")
    x = proxy[2:n]
    data[["x"]]=x
  }
  if(!is.null(events)){
    eventindexes = numeric(length(events))

    if(tolower(eventmeasure) %in% c("depth","z")){
      eventindexes = which.index (events, z)
    }else if(tolower(eventmeasure) %in% c("age","time","y")){
      eventindexes = which.index (events, y)
    }
    # for(i in 1:length(events)){
    #   if(tolower(eventmeasure) %in% c("depth","z")){
    #     eventindexes = which.index (events, record)
    #     if(events[i] < min(z) || events[i]>max(z)){
    #       warning(paste0("Event ",i,", located at ",events[i]," is outside the interval covered by 'depth' (",min(z),", ",max(z),"). The event will be omitted!"))
    #       eventindexes[i] = NA
    #     }else{
    #       eventindexes[i] = which(abs(events[i]-z) == min(abs(events[i]-z)))
    #     }
    #   }else if(tolower(eventmeasure) %in% c("age","time","y")){
    #     if(events[i] < min(y) || events[i]>max(y)){
    #       warning(paste0("Event ",i,", located at ",events[i]," is outside the interval covered by 'age' (",min(age),", ",max(age),"). The event will be omitted!"))
    #       eventindexes[i] = NA
    #     }else{
    #       eventindexes[i] = which(abs(events[i]-y) == min(abs(events[i]-y)))
    #     }
    #   }
    #
    # }
    eventindexes = unique(c(1,eventindexes[!is.na(eventindexes)]))
    nevents = length(eventindexes)
  }
  if(reg.model$psi0) {
    if(is.null(events)){
      stop("'events' are missing. Shutting down...")
    }


    for(i in 2:nevents){
      ev = numeric(n-1)
      ev[ eventindexes[i-1]:(eventindexes[i]-1) ] = data$z[eventindexes[i-1]:(eventindexes[i]-1)]
      data[[paste0("a",i-1)]] = ev
      konst = numeric(n-1)
      konst[eventindexes[i-1]:(eventindexes[i]-1)] = 1
      data[[paste0("c",i-1)]] = konst
      formulastring = paste0(formulastring, " + a",i-1," + c",i-1)
    }
  }
  if(tolower(method)=="inla") data$idy=data$z

  object = list(data=data)
  #object$model = list(reg.model=reg.model,noise=noise,method=method,eventmeasure=eventmeasure)
  object$.args = list(reg.model = reg.model, noise=noise,method=method, eventmeasure=eventmeasure)
  object$.args$ls.formulastring = formulastring
  object$.args$eventindexes = eventindexes
  object$ls.formula = as.formula(formulastring)
  object$.args$nevents=nevents
  object$preceeding = list(y0 = age[1],z0=depth[1],x0=proxy[1])

  if(tolower(method)=="inla"){
    if(tolower(noise) %in% c("iid","independent","ar(0)")){
      formulastring = paste0(formulastring, " + f(idy,model=\"iid\")")
    }else if(tolower(noise) %in% c("ar1","ar(1)")){
      formulastring = paste0(formulastring, " + f(idy,model=\"ar1\")")
    }else if(tolower(noise) %in% c("ar2","ar(2)")){
      formulastring = paste0(formulastring, " + f(idy,model=\"ar\",order=2)")
    }
  }
  formula = as.formula(formulastring)
  object$formula = formula
  object$.args$formulastring=formulastring
  object$.args$reference.label=reference.label

  ## move to separate function

  return(object)

}
