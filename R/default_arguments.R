#' Default variables in events
#'
#' Sets the default variables in the list \code{events} used to specify the climatic
#' periods and separating events used in the linear predictor. The list contains the following arguments:
#' \itemize{
#'   \item{\code{locations} }{Numeric describing the location of the transitions between climatic periods.
#'    Default value is \code{NULL}, and must be provided by user.}
#'   \item{\code{locations_unit} }{Character describing on which axis the locations are given. Can be \code{age},
#'   \code{depth} or \code{index}. Default value is \code{depth}.}
#'   \item{\code{degree} }{Integer giving the polynomial degree the individual period-specific
#'   covariates should be in. Set to \code{0} to only include constants between transitions,
#'   \code{1} to also include a linear function between covariates and \code{2} to include a square
#'   function. Default value is \code{1}.}
#'   \item{\code{fill_data} }{Boolean telling whether or not the \code{data} object should be filled to include
#'   the climate period covariates. Default value is \code{TRUE}.}
#'   \item{\code{fill_formula} }{Boolean telling whether or not the \code{formula} object should be
#'   filled to include description of the period-specific covariates (with degree
#'   as given in to \code{events\$degree}.). Default value is \code{TRUE}.}
#' }
#' @return Returns a list including default values for all variables in \code{events}.
#'
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{bremla},\link{set.options}}
#' @keywords bremla events default
#'
events.default <- function(){
  return(list(
       locations=NULL,
       locations_unit="depth",
       degree=1, #can be 0, 1 or 2
       fill_data=TRUE,
       fill_formula=TRUE
       )
  )
}

#' Default variables in control.fit
#'
#' Sets the default variables in the list \code{control.fit} used to specify the
#' fitting procedure. The list contains the following arguments:
#' \itemize{
#'   \item{\code{noise} }{Character which noise process is assumed for the residuals.
#'   \code{"iid"}, \code{"ar1"} and \code{"ar2"} are currently. However, if synchronization is to
#'   be performed only \code{"ar1"} is currently supported. Default input is \code{"ar1"}.}
#'   \item{\code{method} }{Character describing which method should be used to fit the process.
#'   Only \code{inla} is currently supported, but more will hopefully come in the future.}
#'   \item{\code{verbose} }{Boolean telling whether or not the \code{inla}-program should run in
#'   verbose mode. Default \code{FALSE}.}
#'   \item{\code{ncores} }{Integer used to set the number of cores to be used in the \code{inla} fitting
#'   procedure. Default value is \code{1}.}
#'   \item{\code{transform} }{Character indicating which transformation should be applied to
#'   the observations. Can be \code{"log"} or \code{"identity"}. If synchronization is
#'   intended then only \code{"identity"} is currently supported, which is also the
#'   default value here.}
#' }
#' @return Returns a list including default values for all variables in \code{control.fit}.
#'
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{bremla},\link{bremla_modelfitter},\link{set.options}}
#' @keywords bremla fitting default
#'
control.fit.default <- function(){
  return(list(
    noise="ar1",
    method="inla",
    verbose=FALSE,
    log.theta.prior=NULL,
    ncores=1,
    transform="identity"
    )
  )
}
#' Default variables in control.sim
#'
#' Sets the default variables in the list \code{control.sim} used to specify the
#' fitting procedure. The list contains the following arguments:
#' \itemize{
#'   \item{\code{synchronized} }{Boolean indicating whether or not synchronization is to be performed.
#'   Default value is \code{TRUE}, but this is set to \code{FALSE} if supplied tie-points are not found.}
#'   \item{\code{nsims} }{Integer giving how many simulations should be generated. Default value is \code{10000}.}
#'   \item{\code{ncores} }{Integer which sets the number of cores to be used to sample from
#'   the joint posterior obtained by INLA. Default value is \code{2}.}
#'   \item{\code{store.everything} }{Boolean indicating whether or not to store superfluous information
#'   such as the simulated mean vectors (from fixed effects). Default \code{FALSE}.}
#'   \item{\code{summary} }{List object with two items
#'     \itemize{
#'     \item{\code{compute} }{Boolean indicating if summary statistics (posterior marginal means median
#'     and credible intervals) should be computed. Default \code{TRUE}.}
#'     \item{\code{CI.type} }{character giving which type of credible intervals to be computed.
#'     \code{"quantiles"} for quantiles based intervals (Default) or \code{"hpd"} for
#'     highest posterior density intervals.}}
#'     }
#' }
#' @return Returns a list including default values for all variables in \code{control.sim}.
#'
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{bremla},\link{bremla_chronology_simulation},
#' \link{bremla_synchronized_simulation},\link{set.options}}
#' @keywords bremla simulation default
#'
control.sim.default <- function(){
  return(list(
    synchronized = TRUE,
    nsims=10000,
    ncores=2,
    store.everything=FALSE,
    summary=list(
      compute=TRUE,
      CI.type="quantiles"
      )
    )
  )
}

#' Default variables in synchronization
#'
#' Sets the default variables in the list \code{synchronization} used to specify the
#' tie-points used for synchronization. The list contains the following arguments:
#' \itemize{
#'   \item \code{locations} Numeric describing the locations of the tie-points. Default value
#'   is here \code{NULL} meaning that this has to be given by the user.
#'   \item \code{locations_unit} Character giving which axis the \code{locations} are given.
#'   Can be \code{"depth"} (default), \code{"age"} or \code{"index"}.
#'   \item \code{samples} data.frame giving pre-computed samples. This should be used if
#'   you want to provide the sampled tie-points yourself. \code{nrow} must correspond to
#'   the number of samples and \code{ncol} to the number of tie-points.
#'   Default value is \code{NULL}.
#'   \item \code{method} Character giving which method should be used to generate tie-points.
#'   Currently supported: \code{"adolphi"} (default) uses tie-point distributions given by
#'   Adolphi et al. (2018) and Muscheler et al. (2020), \code{"normal"} uses a normal
#'   distribution with mean and standard deviations supplied by user (TODO).
#'   If \code{"adolphi"} is used, then \code{locations}, \code{locations_unit} is
#'   filled out accordingly.
#'   \item \code{x.ref} Numeric giving a reference for the x-axis (given by \code{locations_unit})
#'   that will be included in future \code{plot}.
#'   \item \code{params} List object containing parameters for parametric tie-point distributions.
#'   for \code{method="gauss"} the list has two numeric objects called \code{mean} and \code{sd}
#'   representing the mean and standard deviation for each tie point, respectively.
#' }
#' @return Returns a list including default values for all variables in
#' \code{synchronization}.
#'
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{bremla}, \link{tiepointsimmer},
#' \link{bremla_synchronized_simulation},\link{set.options}}
#' @keywords bremla simulation default
#'
synchronization.default <- function(){
    return(list(
      locations=NULL,
      locations_unit="depth",
      samples=NULL,
      nsims=10000,
      method="adolphis",
      x.ref = NULL,
      params=NULL
    )
  )
}

#' Default variables in control.linramp
#'
#' Sets the default variables in the list \code{control.linramp} used to specify the
#' linear ramp model fitting procedure for a given transition. The list contains the following arguments:
#' \itemize{
#'   \item \code{label} Character used to describe the name of the transition in future \code{plot}. Default \code{NULL}.
#'    \item \code{proxy} Numeric giving the proxy data to perform the linear ramp model fit.
#'    Default value is \code{NULL} meaning that it must be provided by the user.
#'    \item \code{interval} Interval describing the window of \code{proxy} to perform the
#'    linear ramp model fit. Default value is \code{NULL} so it must be given by the user.
#'    \item \code{interval.unit} Character giving which axis is used for interval. Can be
#'    \code{"depth"} (default), \code{"age"} or \code{"index"}.
#'    \item \code{depth.label} Character describing what label to place on x-axis (depth) in a
#'    future \code{plot}. Default value is \code{""}.
#'    \item \code{depth.ref} Numeric giving a reference value for the depth which will be
#'    included in a future \code{plot}.
#'    \item \code{opt.params} Numeric describing initial values for \code{\link{optim}} function
#'    used to find good initial values for the \code{inla} fitting procedure. If
#'    this is set to \code{NULL} (default) the code will provide its own initial values.
#'    \item \code{compute.t1} Boolean indicating whether or not transition end point \code{t1}
#'    should be computed. Default \code{TRUE}.
#'    \item \code{nsims} Integer giving how many Monte Carlo samples should be used to
#'    compute transition end point \code{t1} if \code{compute.t1} is \code{TRUE}.
#'    \item \code{rescale.y.factor} If proxies (y-axis) should be rescaled, this can be done using
#'    this command. \code{y_new = y_old*rescale.y.factor}.
#'    \item \code{imp.fit} Boolean indicating whether initial values for \code{inla} be computed using \code{optim}.
#'   \item{\code{log.theta.prior}}{Function argument which returns the logarithmic value of the joint
#'   prior of the rescaled hyperparameters \code{theta} used in the linear ramp model.
#'   If \code{NULL} (default) then default priors  will be assumed (see code in
#'   \code{rgeneric.r} file).}
#'    \item \code{h} Step length of optimization procedure in the \code{inla}-program.
#'    Default value is \code{0.01}.
#'    \item \code{ncores} Integer giving the number of cores used in the \code{inla}-program.
#'    For \code{rgeneric} models which are used here, stability can sometimes be
#'    improved if the number of cores is set to \code{1}, which is why it is the
#'    default here.
#'    \item \code{silent} If equal to \code{TRUE} or \code{1L}, then \code{inla} will run in "silent" mode,
#'    If equal to \code{2L} then error messages will also be suppressed.
#'    \item \code{verbose} Boolean stating whether or not the \code{inla}-program should
#'    run in verbose mode. Default \code{FALSE}.
#' }
#' @return Returns a list including default values for all variables in \code{control.linramp}.
#'
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{bremla}, \link{linrampfitter},
#' \link{events_depth_to_age},\link{set.options}}
#' @keywords bremla transition onset default
control.linramp.default <- function(){
  return(list(
      label=NULL,
      proxy=NULL,
      interval=NULL,
      interval.unit="depth",
      depth.label="",
      depth.ref = NULL, #function will allow this
      opt.params=NULL, #function will recognize this
      compute.t1=TRUE,
      nsims=30000,
      rescale.y.factor=1,
      imp.fit=TRUE,
      log.theta.prior=NULL,
      h=0.01,
      ncores=1,
      silent=FALSE,
      verbose=FALSE
      )
    )
}

#' Default variables in control.transition_dating
#'
#' Sets the default variables in the list \code{control.transition_dating} used to
#' specify the procedure for estimating the dating uncertainty of a given transition
#' following a linear ramp model fit.
#' The list contains the following arguments:
#'    \itemize{
#'      \item \code{label} Character used to describe the name of the transition in
#'      future \code{plot}. Default \code{NULL}.
#'      \item \code{nsims} Integer giving the number of samples to be generated for
#'      onset transition age. Default value is \code{10000}.
#'      \item \code{sync} Boolean indicating if the synchronized age-depth model should
#'      be used. Default value is \code{TRUE}, but this is set to \code{fALSE} if
#'      synchronized samples are not yet computed.
#'      \item \code{age.label} Character describing what label to place on x-axis (age)
#'      in a future \code{plot}.
#'      Default value is \code{""}.
#'      \item \code{age.ref} Numeric giving a reference value for the age which will be
#'      included in a future \code{plot}.
#'      }
#' @return Returns a list including default values for all variables in
#' \code{control.transition_dating}.
#'
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{bremla}, \link{linrampfitter},
#' \link{events_depth_to_age},\link{set.options}}
#' @keywords bremla transition onset default
control.transition_dating.default <- function(){
  return(list(
      label=NULL,
      nsims=10000,
      sync=TRUE,
      age.label="",
      age.ref=NULL #function will allow this
    )
    )
}


#' Default variables in control.bias
#'
#' Sets the default variables in the list \code{control.bias} used to specify the
#' how an unknown stochastic bias should be included. The list contains the following arguments:
#' \itemize{
#'   \item{\code{bias.model} }{Character describing which model should be used for the unknown bias.
#'   Currently only \code{"uniform"} is supported.}
#'   \item{\code{biasparams} }{Gives the parameters for the \code{bias.model}. If \code{ncol}>1,
#'   the function will be applied to each column consecutively.}
#'   \item{\code{nsims} }{Integer giving the number of samples to be produced for the model. Default value
#'   is \code{10000}.}
#'   \item{\code{store.samples} }{Boolean indicating if samples should be stored. Default value is \code{FALSE}.}
#' }
#' @return Returns a list including default values for all variables in \code{control.bias}.
#'
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{bremla}, \link{bremla_biased_chronologies},
#' \link{set.options}}
#' @keywords bremla bias default
#'
control.bias.default <- function(){
  return(list(
    bias.model="uniform",
    biasparams = c(0.99,1.01),
    nsims=10000,
    store.samples = FALSE
  )
)
}

#' Import default arguments
#'
#' Fills out missing arguments in list \code{opt} with default arguments in list
#' \code{default.opt}.
#' @param opt List object with different specifications.
#' @param default.opt List of default variables corresponding to \code{opt}.
#'
#' @return Returns the \code{opt} list, but with values from \code{default.opt} inserted
#' in missing values.
#'
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{bremla},\link{bremla_prepare}}
#' @keywords bremla default
set.options <- function(opt,default.opt){
  temp = default.opt

  if(length(opt)>0){
    for(i in 1:length(opt)){
      if(names(opt)[i] %in% names(default.opt)){
        if(!is.list(opt[[i]])){
          temp[[ names(opt)[i] ]] <- opt[[i]]
        }else{
          for(j in 1:length(opt[[i]])){
            temp[[ names(opt)[i] ]][[names(opt[[i]])[j]]] <- opt[[i]][[j]]
          }
        }
      }else{
        temp[[names(opt)[i]]] <- opt[[i]]
      }
    }
  }
  return(temp)
}


