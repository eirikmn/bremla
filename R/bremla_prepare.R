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
#' covariates (phi) can be filled out using 'events$fill_data'.
#' @param nsims Number of chronologies to be simulated.
#' @param reference.label Character label of reference timescale. Used in \code{\link{plot},\link{summary}}.
#' @param x.label Character label for the x-axis (depth).
#' @param y.label Character label for the y-axis (age).
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
#' depth2 = depth^2/depth[1]^2 #normalize for stability
#'
#' data = data.frame(age=age,dy=dy,proxy=proxy,depth=depth,depth2=depth2)
#' data = rbind(c(y0,NA,NA,z0,NA),data) #First row is only used to extract y0 and z0.
#'
#' events=list(locations=c(1210,1220,1240))
#' control.fit = list(ncores=2,noise="ar1")
#' synchronization=list(method="gauss")
#' control.sim=list(synchronized=2,
#'                  summary=list(compute=TRUE))
#'
#' control.bias=NULL
#' object = bremla_prepare(formula,data,nsims=5000,reference.label="simulated timescale",
#'                         events = events,
#'                         synchronization=synchronization,
#'                         control.fit=control.fit,
#'                         control.sim=control.sim,
#'                         control.bias=control.bias)
#' summary(object)
#' @export
bremla_prepare = function(formula,data,nsims=NULL,reference.label=NULL,
                          x.label=NULL,y.label=NULL,
                          events = NULL,
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
  if(is.null(control.fit)) {
    control.fit=list(noise="ar1")

  }

  control.fit = set.options(control.fit,control.fit.default())

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
    #colnames(data$y)="age"
    colnames(data)[which("y"==colnames(data))]="age"
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
    }else{
      eventindexes = events$locations
    }
    #eventindexes = unique(c(1,eventindexes[!is.na(eventindexes)]))
    eventindexes = unique(c(1,eventindexes[!is.na(eventindexes)],length(y)))
    nevents = length(eventindexes)-1

    events$nevents=nevents
    events$eventindexes=eventindexes
  }else{
    nevents=0
    eventindexes=NULL
  }




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



  if(!is.null(nsims)){
    if(!is.null(control.sim)) control.sim$nsims = nsims
    if(!is.null(control.bias)) control.bias$nsims = nsims
    if(!is.null(synchronization)) synchronization$nsims = nsims
  }

  if(!is.null(control.fit)){ #add random effect to formula string for use in INLA
    if(tolower(control.fit$method) %in% "inla" ){
      if(!is.null(control.fit$rgeneric)){ #rgeneric model is specified

        control.fit$noise = "rgeneric"

        ntheta = length(control.fit$rgeneric$from.theta)

        model.rgeneric = INLA::inla.rgeneric.define(rgeneric.fitting, n=n-1,
                                              ntheta=ntheta,
                                              Q = control.fit$rgeneric$Q,
                                              log.prior=control.fit$rgeneric$log.prior)

        formulastring_inla = paste0(formulastring, " + f(idy,model=model.rgeneric)")

        # formula=dy ~ -1 + depth + f(idy,model=model.rgeneric)


      }else{

        hyperprior = control.fit$hyperprior
        if(!is.null(hyperprior)){

            hyperstring = paste0(", hyper = list(")
            hypernames = names(hyperprior)
            n_hyper = length(hypernames)
            liststring = paste(hyperprior)

            for( iter in 1:n_hyper){
              if(iter > 1){
                hyperstring = paste0(hyperstring, ", ")
              }
              itername = hypernames[iter]
              hyperstring = paste0(hyperstring,itername," = ")
              hyperstring = paste0(hyperstring,liststring[iter])
            }
            hyperstring = paste0(hyperstring,")")


          if(tolower(control.fit$noise) %in% c("iid","independent","ar(0)")){
            formulastring_inla = paste0(formulastring, " + f(idy,model=\"iid\"", hyperstring,")")
          }else if(tolower(control.fit$noise) %in% c("ar1","ar(1)")){
            formulastring_inla = paste0(formulastring, " + f(idy,model=\"ar1\"", hyperstring,")")
          }else if(tolower(control.fit$noise) %in% c("ar2","ar(2)")){
            formulastring_inla = paste0(formulastring, " + f(idy,model=\"ar\",order=2", hyperstring,")")
          }

        }else{
          if(tolower(control.fit$noise) %in% c("iid","independent","ar(0)")){
            formulastring_inla = paste0(formulastring, " + f(idy,model=\"iid\")")
          }else if(tolower(control.fit$noise) %in% c("ar1","ar(1)")){
            formulastring_inla = paste0(formulastring, " + f(idy,model=\"ar1\")")
          }else if(tolower(control.fit$noise) %in% c("ar2","ar(2)")){
            formulastring_inla = paste0(formulastring, " + f(idy,model=\"ar\",order=2)")
          }
        }



        model.rgeneric=NULL
      }
    }
    formula_inla=as.formula(formulastring_inla)
  }else{
    formula_inla=as.formula(formulastring)
  }


  object= list(data=data_obj,formula=formula_inla,
               initials=initials,
               original.chron = data.frame(depth=depth,age=age))
  if(!is.null(events)) object$events=events
  #store input arguments
  str = cleanstring(format(formulastring))

  tildepos = str_locate(str,"~")[1]
  responsename = str_sub(str,start=1L,end=tildepos-1)

  lat.selection = lat.selector(format(formulastring))

  object$.internal=list(formula.ls = as.formula(formulastring),
                        lat.selection=lat.selection)
  object$.args=list(formula.ls=as.formula(formulastring),
                    formula.input=formula.input,
                    responsename=responsename,
                    data.input=data,
                    model.rgeneric=model.rgeneric,
                    reference.label=reference.label,
                    x.label=x.label,
                    y.label=y.label,
                    events=events,
                    synchronization=synchronization,
                    control.fit=control.fit,
                    control.sim=control.sim,
                    control.linramp=control.linramp,
                    control.transition_dating=control.transition_dating,
                    control.bias=control.bias)

if(tolower(control.fit$noise) %in% c("rgeneric","custom")){
  object$.args$model.rgeneric
}
  class(object) = "bremla"

  return(object)

}
