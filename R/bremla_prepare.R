#' Bremla preparation function
#'
#' Prepares input and formats data for bremla and functions therein.
#'
#' @param formula Formula describing the response and covariates (partial
#' covariates (phi) can be filled out using \code{events$fill_formula})
#' @param data data.frame including response, all covariates included in
#' \code{formula}, as well as 'depth' and 'age' describing both axes of reference chronology.
#' First row is used to extract initial values 'y0' and 'z0', the remaining variables
#' in the first row is typically not used and can be left as \code{NA}. Partial
#' covariates (phi) can be filled out using 'events$fill_data'
#' @param reference.label Character label of reference timescale. Used in \code{\link{plot},\link{summary}}.
#' @param nsims Number of chronologies to be simulated.
#' @param events List object describing the specifics of the climate transitions
#' used in 'formula'. Must include an item called \code{locations}. See \code{?events.default} for details.
#' @param synchronization List object describing specifics related to tie-points
#' and their distribution. If \code{synchronization$method="adolphi"} a specific
#' set of instructions are employed. If not, this list must include \code{locations}
#' argument and preferably \code{locations_unit} to determine which axis (default \code{unit="depth"}).
#' See \code{\link{synchronization.default}} for details.
#' @param control.fit List object describing specifics related to the fitting procedure.
#' See \code{\link{control.fit.default}} for details.
#' @param control.sim List object describing specifics related to simulating chronologies.
#' See \code{\link{control.sim.default}} for details.
#' @param control.linramp List containing specifications for fitting a linear ramp model to given data.
#' If used, must include \code{control.linramp\$proxy} and \code{control.linramp\$interval}. See
#' \code{\link{control.linramp.default}} for details.
#' @param control.transition_dating List object describing specifics related to
#' the estimation of a given transition. Must include \code{interval} and \code{proxy}
#' used in linear ramp estimation. See \code{\link{control.transition_dating.default}} for details.
#' @param control.bias List object describing specifics related to the simulation of
#' unknown stochastic bias. See \code{\link{control.bias.default}} for details.
#'
#' @return Returns a list containing data in internal formatting (\code{data}), linear regression formula (\code{ls.formula}), inla formula (\code{formula}), input settings and indices corresponding to climate transitions (\code{.args}) and z_0, y_0 and x_0 (\code{preceeding})
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{bremla},\link{bremla_chronology_simulation}}
#' @keywords bremla preparation
#'
#' @examples
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
#' @export
bremla_prepare = function(formula,data,reference.label=NULL,
                          nsims=NULL,events = NULL,
                          synchronization=NULL,
                          control.fit=NULL,
                          control.sim=NULL,
                          control.linramp=NULL,
                          control.transition_dating=NULL,
                          control.bias=NULL){


  ## fill out missing input arguments with defaults:
  if(!is.null(events)) events = set.options(events,events.default())
  if(!is.null(synchronization)) synchronization = set.options(synchronization,
                                                              synchronization.default())
  if(!is.null(control.fit)) control.fit = set.options(control.fit,control.fit.default())
  if(!is.null(control.sim)) control.sim = set.options(control.sim,control.sim.default())
  if(!is.null(control.linramp)) control.linramp = set.options(
                  control.linramp,control.linramp.default())
  if(!is.null(control.transition_dating)) control.transition_dating = set.options(
                  control.transition_dating,control.transition_dating.default())
  if(!is.null(control.bias)) control.bias = set.options(control.bias,
                                                        control.bias.default())

  if("age" %in% colnames(data)){
    age=data$age
  }else if("y" %in% colnames(data)){
    age=data$y
    colnames(data$y)="age"
  }else{
    stop("Could not find 'age' in data.frame. Stopping!")
  }

  if("depth" %in% colnames(data)){
    depth=data$depth
  }else if("z" %in% colnames(data)){
    depth=data$z
    colnames(data$z)="depth"
  }else{
    stop("Could not find 'depth' in data.frame. Stopping!")
  }

  n = length(age)
  initials = data[1,]
  data_obj = data[2:n,]

  z = depth[2:n]
  y = age[2:n]
  formula.input=formula

  if(!is.null(events)){
    eventindexes = numeric(length(events))
    if(tolower(events$locations_unit) %in% c("depth","z")){ #find indexes corresponding to 'events' location
      eventindexes = which.index (events$locations, z)
    }else if(tolower(events$locations_unit) %in% c("age","time","y")){
      eventindexes = which.index (events$locations, y)
    }
    eventindexes = unique(c(1,eventindexes[!is.na(eventindexes)]))
    nevents = length(eventindexes)
  }else{
    nevents=0
    eventindexes=NULL
  }
  events$nevents=nevents
  events$eventindexes=eventindexes



  formulastring = format(formula.input)

  if(!is.null(events)){
    if(events$fill_formula){
      for(i in 2:nevents){
        formulastring = paste0(formulastring, " + psi0_",i-1)
        if(events$degree>=1) formulastring=paste0(formulastring," + psi1_",i-1)
        if(events$degree==2) formulastring=paste0(formulastring," + psi2_",i-1)
      }
      formulastring = paste0(formulastring, " + psi0_",nevents)
      if(events$degree>=1) formulastring=paste0(formulastring," + psi1_",nevents)
      if(events$degree==2) formulastring=paste0(formulastring," + psi2_",nevents)
    }
    if(events$fill_data){
      for(i in 2:nevents){
        konst = numeric(n-1)
        konst[eventindexes[i-1]:(eventindexes[i]-1)] = 1
        data_obj[[paste0("psi0_",i-1)]] = konst

        if(events$degree>=1){
          ev1 = numeric(n-1)
          ev1[ eventindexes[i-1]:(eventindexes[i]-1) ] = z[eventindexes[i-1]:(eventindexes[i]-1)]
          data_obj[[paste0("psi1_",i-1)]] = ev1
        }
        if(events$degree == 2){
          ev2 = numeric(n-1)
          ev2[ eventindexes[i-1]:(eventindexes[i]-1) ] = z[eventindexes[i-1]:(eventindexes[i]-1)]^2
          data_obj[[paste0("psi2_",i-1)]] = ev2
        }

        konst = numeric(n-1)
        konst[eventindexes[nevents]:(n-1)] = 1
        data_obj[[paste0("psi0_",nevents)]] = konst

        if(events$degree>=1){
          ev1 = numeric(n-1)
          ev1[ eventindexes[nevents]:(n-1) ] = z[eventindexes[nevents]:(n-1)]
          data_obj[[paste0("psi1_",nevents)]] = ev1
        }
        if(events$degree == 2){
          ev2 = numeric(n-1)
          ev2[ eventindexes[nevents]:(n-1) ] = z[eventindexes[nevents]:(n-1)]^2
          data_obj[[paste0("psi2_",nevents)]] = ev2

      }
    }
    }
  }



  if(!is.null(nsims)){
    if(!is.null(control.sim)) control.sim$nsims = nsims
    if(!is.null(control.bias)) control.bias$nsims = nsims
    if(!is.null(synchronization)) synchronization$nsims = nsims
  }

  if(!is.null(control.fit)){ #add random effect to formula string for use in INLA
    if(tolower(control.fit$method) %in% "inla" ){
      if(tolower(control.fit$noise) %in% c("iid","independent","ar(0)")){
        formulastring_inla = paste0(formulastring, " + f(idy,model=\"iid\")")
      }else if(tolower(control.fit$noise) %in% c("ar1","ar(1)")){
        formulastring_inla = paste0(formulastring, " + f(idy,model=\"ar1\")")
      }else if(tolower(control.fit$noise) %in% c("ar2","ar(2)")){
        formulastring_inla = paste0(formulastring, " + f(idy,model=\"ar\",order=2)")
      }
    }
    formula=as.formula(formulastring_inla)
  }else{
    formula=as.formula(formulastring)
  }

  object= list(data=data_obj,formula=formula,
               initials=initials,
               original.chron = data.frame(depth=depth,age=age))
  if(!is.null(events)) object$events=events
  #store input arguments
  str = cleanstring(format(formulastring))
  tildepos = str_locate(str,"~")
  responsename = str_sub(str,start=1L,end=tildepos-1)

  lat.selection = lat.selector(format(formulastring))
  object$.internal=list(formula.ls = as.formula(formulastring),
                        lat.selection=lat.selection)
  object$.args=list(formula.ls=as.formula(formulastring),
                    formula.input=formula.input,
                    responsename=responsename,
                    data.input=data,
                    reference.label=reference.label,
                    events=events,
                    synchronization=synchronization,
                    control.fit=control.fit,
                    control.sim=control.sim,
                    control.linramps=control.linramp,
                    control.transition_dating=control.transition_dating,
                    control.bias=control.bias)


  class(object) = "bremla"

  return(object)

}

#
#   #
#
#
#
#
#
#
#
#
#   #
#
#
#   y = age[2:n]
#   dy = diff(age)
#   z = depth[2:n]
#
#   data = data.frame(y=y,dy=dy,z=z) #format dataset into data.frame object
#
#
#
#   ## create formula string to be used for least squares (and INLA)
#   if(transform %in% c("log","logarithmic")){
#     data$logdy = log(dy)
#     formulastring = "logdy ~ "
#   }else{
#     formulastring = "dy ~ "
#   }
#
#
#
#   ## model components are expressed using reg.model
#
#   if(!reg.model$const) formulastring=paste0(formulastring,"-1") else formulastring=paste0(formulastring,"1")
#   if(reg.model$depth1) formulastring=paste0(formulastring," + z")
#   if(reg.model$depth2) {
#     formulastring=paste0(formulastring," + z2")
#     data[["z2"]]=z^2
#   }
#   if(reg.model$proxy) {
#     if(missing(proxy)){
#       stop("'proxy' is missing. Shutting down...")
#     }
#     formulastring=paste0(formulastring," + x")
#     x = proxy[2:n]
#     data[["x"]]=x
#   }
#   if(!is.null(events)){
#     eventindexes = numeric(length(events))
#
#     if(tolower(eventmeasure) %in% c("depth","z")){ #find indexes corresponding to 'events' location
#       eventindexes = which.index (events, z)
#     }else if(tolower(eventmeasure) %in% c("age","time","y")){
#       eventindexes = which.index (events, y)
#     }
#
#     eventindexes = unique(c(1,eventindexes[!is.na(eventindexes)]))
#     nevents = length(eventindexes)
#   }
#
#   if(reg.model$psi0 && !is.null(events)) {
#       for(i in 2:nevents){ #express covariates for psi-functions (for each climate period)
#         ev = numeric(n-1)
#         ev[ eventindexes[i-1]:(eventindexes[i]-1) ] = data$z[eventindexes[i-1]:(eventindexes[i]-1)]
#         data[[paste0("a",i-1)]] = ev
#         konst = numeric(n-1)
#         konst[eventindexes[i-1]:(eventindexes[i]-1)] = 1
#         data[[paste0("c",i-1)]] = konst
#         formulastring = paste0(formulastring, " + a",i-1," + c",i-1)
#       }
#   }
#
#   if(tolower(method)=="inla") data$idy=data$z
#
#   ## create main object which will include all input and output from function
#   object = list(data=data)
#   object$.args = list(reg.model = reg.model, noise=noise,method=method)
#   object$.args$ls.formulastring = formulastring
#
#   if(!is.null(events)){
#     object$.args$eventindexes=eventindexes
#     object$.args$events=events
#     object$.args$eventmeasure=eventmeasure
#     object$.args$nevents=nevents
#   }else{
#     object$.args$eventindexes=NULL
#     object$.args$events=NULL
#     object$.args$eventmeasure="missing"
#     object$.args$nevents=0
#   }
#
#   object$ls.formula = as.formula(formulastring)
#
#   object$preceeding = list(y0 = age[1],z0=depth[1],x0=proxy[1])
#
#
#
#
#   if(tolower(method)=="inla"){ #add random effect to formula string for use in INLA
#     if(tolower(noise) %in% c("iid","independent","ar(0)")){
#       formulastring = paste0(formulastring, " + f(idy,model=\"iid\")")
#     }else if(tolower(noise) %in% c("ar1","ar(1)")){
#       formulastring = paste0(formulastring, " + f(idy,model=\"ar1\")")
#     }else if(tolower(noise) %in% c("ar2","ar(2)")){
#       formulastring = paste0(formulastring, " + f(idy,model=\"ar\",order=2)")
#     }
#   }
#
#
#
#
#   formula = as.formula(formulastring)
#   object$formula = formula
#   object$.args$formulastring=formulastring
#   object$.args$reference.label=reference.label
#   object$.args$transform = transform
#   object$.args$proxy.type = proxy.type

